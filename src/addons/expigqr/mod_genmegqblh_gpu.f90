MODULE mod_genmegqblh_gpu

  USE ISO_C_BINDING, ONLY: C_LOC, C_PTR
  USE mod_prec
  USE mod_gpu
  USE modmain, ONLY: zzero, zone, natmtot, natmcls, nstsv, nkptnr, &
                     ngrid, ngrtot, ngkmax, ivg, &
                     lmmaxapw, nufrmax, nspinor, spinpol, igfft, ivgig, cfunir
  USE mod_mpi_grid
  USE mod_timer
  USE mod_papi
  USE mod_prec
  USE mod_addons, ONLY: debug_level, dim_k, ias2ic, &
                        pt_megqblh, pt_megqblh_mt, pt_megqblh_it
  USE mod_addons_q, ONLY: ngq, igqig, sfacgq
  USE mod_nrkp, ONLY: spinor_ud

#ifdef _PACK_gntuju_
  USE mod_expigqr, ONLY: expigqr22, gntuju_packed, megqblh, bmegqblh, &
                         nmegqblh, idxkq, nbandblhloc, &
                         ltranblhloc, ntranblhloc, idxtranblhloc, &
                         nrownz, ncolnz, npackdim, irownz, icolnz, irowmap_wf1
#else
  USE mod_expigqr, ONLY: expigqr22, gntuju, megqblh, bmegqblh, &
                         nmegqblh, idxkq, nbandblhloc, &
                         ltranblhloc, ntranblhloc, idxtranblhloc
#endif /* _PACK_gntuju_ */
                         
#ifdef _USE_3M_
  USE mod_lapack, ONLY: ZGEMM3M
#else
  USE mod_lapack, ONLY: ZGEMM
#endif /* _USE_3M_ */

  IMPLICIT NONE
  
  ! Table of band indices per spin projection (replaces l1 check)
  ! First dimension is nstsv, but only the first nj elements will be used
  ! Second dimension is nspinor ( = 1 when spinpol == .FALSE, = 2 when .TRUE. )
  ! Third dimension is nkptnrloc
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: jbandidx

  ! Number of 2nd-variational states per spin projection
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: njbands
  INTEGER :: nj, njmax

  ! Number of muffin-tin elements
#ifdef _PACK_gntuju_
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nmt
#endif
  INTEGER :: nmtmax

  ! Block size for batched ZGEMM
  !INTEGER, PARAMETER :: nb = 64

  ! Number of blocks
  INTEGER :: nblock1, nblock2
  
  ! Number of batches
  INTEGER :: nbatch1, nbatch2

  ! Number of bands associated with the bra state vectors
  INTEGER :: nband1

  ! Number of G+q vectors for a particular value of q-vector
  INTEGER :: ngqiq

  ! Translation table for each batch index
  ! Dimensions are natmtot, ngqiq, nblock, respectively
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: batchidx

  ! Matrices for muffin-tin calculation (batched ZGEMM)
  COMPLEX(KIND=dz), DIMENSION(:,:,:), ALLOCATABLE :: bgntuju, b1, b2

  ! Device pointers related to the muffin-tin calculation
  ! (Only used for CUDA version)
  TYPE(C_PTR), DIMENSION(:,:,:), ALLOCATABLE :: d_b1, d_b2
  TYPE(C_PTR), DIMENSION(:,:,:), ALLOCATABLE :: d_gntuju, d_sfacgq, d_wfsvmt1

  ! Array of device pointers for batching the muffin-tin calculation
  TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: dptr_b1, dptr_b2, dptr_gntuju

  ! Temporary array to hold results for muffin-tin calculation (batched ZGEMM)
  ! (will be removed after everything is ported to GPU)
  COMPLEX(KIND=dz), DIMENSION(:,:,:,:), ALLOCATABLE :: wftmp1mt

  ! Temporary arrays to hold results for final calculation (ZGEMM3M)
  COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: wftmp1, wftmp2

  ! Device pointers related to the final calculation (ZGEMM3M)
  ! (Only used for CUDA version)
  TYPE(C_PTR), DIMENSION(:,:,:), ALLOCATABLE :: d_wfsvmt2, d_wfir1

CONTAINS

!==============================================================================
! Sets constant module variables on CPU and copies/allocates them on device

  SUBROUTINE genmegqblh_allocmodvar_const( ikloc, iq )
    USE modmain, ONLY: natmtot, nspinor, nstsv, lmmaxapw, nufrmax
    USE mod_addons, ONLY: nkptnrloc
    USE mod_addons_q, ONLY: ngq
    USE mod_expigqr, ONLY: nbandblhloc
    IMPLICIT NONE

    ! Input arguments
    INTEGER, INTENT(IN) :: ikloc, iq

    ! Number of G+q vectors for a particular value of q-vector
    ngqiq = ngq(iq)

    ! Convert to the corresponding ist1 loop in getmeidx()
    ! as stored into nbandblhloc
    nband1 = nbandblhloc(ikloc)

    ! Number of blocks and batches, blocked version
    !nblock1 = CEILING( REAL(nband1)/REAL(nb) )
    !nbatch1 = ngqiq * natmtot * nblock1

    ! Number of blocks and batches, unblocked version
    nblock1 = 1
    nbatch1 = ngqiq * natmtot

    ! Note: Number of muffin-tin elements nmt is an array when _PACK_gntuju_
    !       is enabled. The maximum value of this array is called nmtmax
    !       and is set in gengntuju()

    IF( iq == 1 .AND. ikloc == 1 ) THEN
#ifdef _OPENACC
       ! Copy/allocate constants to device
       !$ACC ENTER DATA CREATE( nj, njmax, nblock1, nbatch1, nband1, ngqiq ) &
       !$ACC            COPYIN( natmtot, nspinor, nstsv, lmmaxapw, nufrmax, &
       !$ACC                    nkptnrloc, ngq, nmtmax )

#elif defined(_CUDA_)

       ! Allocate constants on device
       !CALL cudaMalloc( d_natmtot, ... )
       !CALL cudaMalloc( d_nstsv, ... )
       !CALL cudaMalloc( d_nspinor, ... )
       !CALL cudaMalloc( d_lmmaxapw, ... )
       !CALL cudaMalloc( d_nufrmax, ... )
       !CALL cudaMalloc( d_nj, ... )
       !CALL cudaMalloc( d_nblock1, ... )
       !CALL cudaMalloc( d_nbatch1, ... )
       !CALL cudaMalloc( d_nband1, ... )
       !CALL cudaMalloc( d_ngqiq, ... )

       ! Copy constants H->D
       !CALL cudaMemcpy( ... )

#endif /* _OPENACC || _CUDA_ */
    END IF ! iq == 1 && ikloc == 1

    ! Copy variables that are dependent on outer loop vars { ikloc, iq }
#ifdef _OPENACC
    !$ACC UPDATE DEVICE( nblock1, nbatch1, nband1, ngqiq )
#elif defined(_CUDA_)
    !CALL cudaMemcpy( ... )
#endif /* _OPENACC || _CUDA_ */

    RETURN
  END SUBROUTINE genmegqblh_allocmodvar_const

!==============================================================================
! Cleans up constant module variables on device

  SUBROUTINE genmegqblh_freemodvar_const( ikloc, iq )
    USE mod_addons, ONLY: nkptnrloc
    USE mod_addons_q, ONLY: nvqloc
    IMPLICIT NONE

    ! Input arguments
    INTEGER, INTENT(IN) :: ikloc, iq

    IF( iq == nvqloc .AND. ikloc == nkptnrloc ) THEN
#ifdef _OPENACC

       !$ACC EXIT DATA DELETE ( natmtot, nspinor, nstsv, lmmaxapw, nufrmax, &
       !$ACC                    nkptnrloc, ngq, nmtmax, &
       !$ACC                    nj, njmax, nblock1, nbatch1, nband1, ngqiq )

#elif defined(_CUDA_)

       !CALL cudaFree( ... )

#endif /* _OPENACC || _CUDA_ */
    END IF ! iq == nvqloc && ikloc == nkptnrloc

    RETURN
  END SUBROUTINE genmegqblh_freemodvar_const

!==============================================================================
! Allocates module variables on CPU and device (GPU)
! that are spin-dependent, i.e., related to genmegqblh_countbands() kernel

  SUBROUTINE genmegqblh_allocmodvar_spin( ikloc, iq, ispn )
    USE modmain, ONLY: nspinor, nstsv
    USE mod_addons, ONLY: nkptnrloc
    IMPLICIT NONE

    ! Input arguments
    INTEGER, INTENT(IN) :: ikloc, iq, ispn

    IF( iq == 1 .AND. ikloc == 1 .AND. ispn == 1 ) THEN

       ! Allocate array for number of states per spin projection on CPU memory
       ALLOCATE( njbands(nspinor,nkptnrloc) )

       ! Allocate array for table of states per spin projection on CPU memory
       ALLOCATE( jbandidx(nstsv,nspinor,nkptnrloc) )

#ifdef _OPENACC

       ! Allocate arrays on device
       !$ACC ENTER DATA CREATE( jbandidx, njbands )

#elif defined(_CUDA_)

       ! Allocate arrays on device
       !CALL cudaMalloc( jbandidx, ... )
       !CALL cudaMalloc( njbands, ... )

       ! Zero the array
       !CALL cudaMemset( ... )

#endif /* _OPENACC || _CUDA_ */

    END IF ! iq == 1 && ikloc == 1 && ispn = 1

    RETURN
  END SUBROUTINE genmegqblh_allocmodvar_spin

!==============================================================================
! Cleans up module variables on CPU and device (GPU)
! that are spin-dependent, i.e., related to genmegqblh_countbands() kernel

  SUBROUTINE genmegqblh_freemodvar_spin( ikloc, iq, ispn )
    USE modmain, ONLY: nspinor
    USE mod_addons, ONLY: nkptnrloc
    USE mod_addons_q, ONLY: nvqloc
    IMPLICIT NONE

    ! Input arguments
    INTEGER, INTENT(IN) :: ikloc, iq, ispn

    IF( iq == nvqloc .AND. ikloc == nkptnrloc .AND. ispn == nspinor ) THEN

#ifdef _OPENACC

       ! Clean up device
       !$ACC EXIT DATA DELETE( jbandidx, njbands )

#elif defined(_CUDA_)

       ! Clean up device
       !CALL cudaFree( jbandidx, ... )
       !CALL cudaFree( njbands, ... )

#endif /* _OPENACC || _CUDA_ */

       ! Clean up CPU memory
       DEALLOCATE( jbandidx )
       DEALLOCATE( njbands )

    END IF ! iq == nvqloc && ikloc == nkptnrloc && ispn == nspinor )

    RETURN
  END SUBROUTINE genmegqblh_freemodvar_spin

!==============================================================================
! Allocates module variables on CPU and device (GPU)
! related to the muffin-tin part calculation (batched ZGEMM)

  SUBROUTINE genmegqblh_allocmodvar_mt()

    USE modmain, ONLY: natmtot, lmmaxapw, nufrmax

#ifdef _PACK_gntuju_
    USE mod_expigqr, ONLY: npackdim
#endif /* _PACK_gntuju_ */

    IMPLICIT NONE

    ! Allocate batching index on CPU
    ALLOCATE( batchidx( natmtot, ngqiq, nblock1 ) )

    ! Allocate temporary array to store results on CPU
    ALLOCATE( wftmp1mt( nmtmax, njmax, natmtot, ngqiq ) )

    ! Allocate batch arrays for the temporary matrices on CPU
    ! Note: we don't use bgntuju in the OpenACC implementation
#ifdef _PACK_gntuju_
       ALLOCATE( b1( npackdim, njmax, nbatch1 ) )
       ALLOCATE( b2( npackdim, njmax, nbatch1 ) )
#else
       ALLOCATE( b1( nmtmax, njmax, nbatch1 ) )
       ALLOCATE( b2( nmtmax, njmax, nbatch1 ) )
#endif /* _PACK_gntuju_ */

#ifdef _OPENACC

    ! Allocate array of device pointers on CPU
    ALLOCATE( dptr_gntuju( nbatch1 ) )
    ALLOCATE( dptr_b1( nbatch1 ) )
    ALLOCATE( dptr_b2( nbatch1 ) )

!--DEBUG
!    ALLOCATE( bgntuju( nmtmax, nmtmax, nbatch1 ))
!--DEBUG

    ! Allocate arrays on device
    !$ACC ENTER DATA CREATE( batchidx, wftmp1mt, &
    !$ACC                    b1, b2, dptr_gntuju, dptr_b1, dptr_b2 )

    sz_wftmp1mt = sz_z * INT(lmmaxapw,KIND=dl) * INT(nufrmax,KIND=dl) &
                       * INT(njmax,KIND=dl) &
                       * INT(natmtot,KIND=dl) * INT(ngqiq,KIND=dl)

#ifdef _PACK_gntuju_
    sz_b1 = sz_z * INT(npackdim,KIND=dl) &
                 * INT(njmax,KIND=dl) * INT(nbatch1,KIND=dl)
    sz_b2 = sz_z * INT(npackdim,KIND=dl) &
                 * INT(njmax,KIND=dl) * INT(nbatch1,KIND=dl)
#else
    sz_b1 = sz_z * INT(nmtmax,KIND=dl) &
                 * INT(njmax,KIND=dl) * INT(nbatch1,KIND=dl)
    sz_b2 = sz_z * INT(nmtmax,KIND=dl) &
                 * INT(njmax,KIND=dl) * INT(nbatch1,KIND=dl)
#endif /* _PACK_gntuju_ */

#elif defined(_CUDA_)

    ! Allocate arrays for "input" data on device
    !CALL cudaMalloc( d_gntuju, ... )
    !CALL cudaMalloc( d_sfacgq, ... )
    !CALL cudaMalloc( d_wfsvmt1, ... )

    ! Transfer data H->D
    !CALL cudaMemcpy( d_wfsvmt1, wfsvmt1, cudaMemcpyHostToDevice )
    !CALL cudaMemcpy( d_sfacgq,  sfacgq,  cudaMemcpyHostToDevice )
    !CALL cudaMemcpy( d_gntuju,  gntuju,  cudaMemcpyHostToDevice )

    ! Allocate batch arrays for the temporary matrices on device
    !CALL cudaMalloc( d_b1, ... )
    !CALL cudaMalloc( d_b2, ... )

    ! Allocate array of device pointers on device
    !CALL cudaMalloc( dptr_gntuju, ... )
    !CALL cudaMalloc( dptr_b1, ... )
    !CALL cudaMalloc( dptr_b2, ... )

#else

    ! OpenMP version needs bgntuju
    ! Allocate batch array for gntuju
#ifdef _PACK_gntuju_
    ALLOCATE( bgntuju( npackdim, npackdim, nbatch1 ))
#else
    ALLOCATE( bgntuju( nmtmax, nmtmax, nbatch1 ))
#endif /* _PACK_gntuju_ */

#endif /* _OPENACC || _CUDA_ */

    RETURN
  END SUBROUTINE genmegqblh_allocmodvar_mt

!==============================================================================
! Cleans up module variables on CPU and device (GPU)
! related to the muffin-tin part calculation (batched ZGEMM)

  SUBROUTINE genmegqblh_freemodvar_mt()
    IMPLICIT NONE

#ifdef _OPENACC

    ! Clean up device
    !$ACC EXIT DATA DELETE( batchidx, wftmp1mt, &
    !$ACC                   b1, b2, dptr_gntuju, dptr_b1, dptr_b2 )

#elif defined(_CUDA_)

    ! Clean up device
    !CALL cudaFree( d_gntuju )
    !CALL cudaFree( d_sfacgq )
    !CALL cudaFree( d_wfsvmt1 )
    !CALL cudaFree( d_b1 )
    !CALL cudaFree( d_b2 )
    !CALL cudaFree( dptr_gntuju )
    !CALL cudaFree( dptr_b1 )
    !CALL cudaFree( dptr_b2 )

#endif /* _OPENACC || _CUDA_ */

    ! Clean up CPU memory
    IF( ALLOCATED(wftmp1mt)    ) DEALLOCATE( wftmp1mt )
    IF( ALLOCATED(b1)          ) DEALLOCATE( b1 )
    IF( ALLOCATED(b2)          ) DEALLOCATE( b2 )
    IF( ALLOCATED(batchidx)    ) DEALLOCATE( batchidx )
#ifdef _OPENACC
    IF( ALLOCATED(dptr_gntuju) ) DEALLOCATE( dptr_gntuju )
    IF( ALLOCATED(dptr_b1)     ) DEALLOCATE( dptr_b1 )
    IF( ALLOCATED(dptr_b2)     ) DEALLOCATE( dptr_b2 )

!--DEBUG
!    IF( ALLOCATED(bgntuju)     ) DEALLOCATE( bgntuju )
!--DEBUG

#elif defined(_CUDA_)

    !CALL cudaFree( ... )

#else
    IF( ALLOCATED(bgntuju)     ) DEALLOCATE( bgntuju )
#endif /* _OPENACC */

    RETURN
  END SUBROUTINE genmegqblh_freemodvar_mt

!==============================================================================
! Kernel 0: Counts how many 2nd-variational states are spin up/down,
!           and returns a list of such states as jbandidx
!==============================================================================

  SUBROUTINE genmegqblh_countbands( ikloc, iq, ispn )

    USE modmain, ONLY: nstsv, spinpol, nspinor
    USE mod_addons, ONLY: nkptnrloc
    USE mod_nrkp, ONLY: spinor_ud
    USE mod_expigqr, ONLY: bmegqblh, nbandblhloc, &
                           idxtranblhloc, ltranblhloc, ntranblhloc
    USE mod_mpi_grid, ONLY: mpi_grid_map
    IMPLICIT NONE

    ! Input arguments
    INTEGER, INTENT(IN) :: ikloc, iq, ispn

    ! Internal variables
    INTEGER, DIMENSION(nkptnrloc) :: iktable
    INTEGER :: is, ik, ikglobal, j, spinproj
    LOGICAL :: lpaired, lspinproj

    IF( iq == 1 .AND. ikloc == 1 .AND. ispn == 1 ) THEN    

       ! Zero out jbandidx
#ifdef _OPENACC
       !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
       !$ACC   PRESENT( nstsv, nspinor, nkptnrloc, jbandidx )
#else
       !$OMP PARALLEL DO
#endif /* _OPENACC */
       DO ik = 1, nkptnrloc
          DO is = 1, nspinor
             DO j = 1, nstsv
                jbandidx(j,is,ik) = 0
             END DO ! nstsv
          END DO ! nspinor
       END DO ! nkptnrloc
#ifdef _OPENACC
       !$ACC END PARALLEL LOOP
#else
       !$OMP END PARALLEL DO
#endif /* _OPENACC */

       ! Zero out njbands
#ifdef _OPENACC
       !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
       !$ACC   PRESENT( nspinor, nkptnrloc, njbands )
#else
       !$OMP PARALLEL DO
#endif /* _OPENACC */
       DO ik = 1, nkptnrloc
          DO is = 1, nspinor
             njbands(is,ik) = 0
          END DO ! nspinor
       END DO ! nkptnrloc
#ifdef _OPENACC
       !$ACC END PARALLEL LOOP
#else
       !$OMP END PARALLEL DO
#endif /* _OPENACC */

       ! Perform ik lookup for global k-point
       !$OMP MASTER
       DO ik = 1, nkptnrloc
          iktable(ik) = mpi_grid_map(nkptnr,dim_k,loc=ik)
       END DO ! nkptnrloc
       !$OMP END MASTER

       IF( spinpol ) THEN

#ifdef _OPENACC
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
          !$ACC   COPYIN( iktable ) &
          !$ACC   PRIVATE( nj, ikglobal ) &
          !$ACC   PRESENT( jbandidx, njbands, &
          !$ACC            nkptnrloc, nspinor, nstsv, spinor_ud, ltranblhloc )
#else
          !$OMP PARALLEL DO COLLAPSE(2) &
          !$OMP   PRIVATE( nj, ikglobal, lpaired, lspinproj )
#endif /* _OPENACC */
          DO ik = 1, nkptnrloc
             DO is = 1, nspinor

                ! Initialize values
                nj = 0
                ikglobal = iktable(ik)

                ! Begin search algorithm for this spin projection
                ! Note: it is desirable to keep the jbandidx indices in order
                !       especially for collinear magnetism
                !       (where j will usually be consecutive)
#ifdef _OPENACC
                !$ACC LOOP SEQ PRIVATE( lpaired, lspinproj )
#endif /* _OPENACC */
                DO j = 1, nstsv

                   ! Test whether the band contributes to transitions
                   lpaired = ltranblhloc(j,ik)

                   ! Test whether it's the "correct" spin projection
                   lspinproj = lpaired .AND. ( spinor_ud(is,j,ikglobal) /= 0 )

                   IF( lspinproj ) THEN

                      ! Increment nj, then save band index j
                      ! Note: ATOMICs not needed since LOOP SEQ is in effect
                      nj = nj + 1
                      jbandidx(nj,is,ik) = j

                   END IF ! lspinproj

                END DO ! nstsv
#ifdef _OPENACC
                !$ACC END LOOP
#endif /* _OPENACC */

                ! Assign njbands
                njbands(is,ik) = nj

             END DO ! nspinor
          END DO ! nkptnrloc
#ifdef _OPENACC
          !$ACC END PARALLEL LOOP
#else
          !$OMP END MASTER
#endif /* _OPENACC */

       ELSE

#ifdef _OPENACC
          !$ACC PARALLEL LOOP GANG VECTOR &
          !$ACC   PRIVATE( nj ) &
          !$ACC   PRESENT( jbandidx, njbands, &
          !$ACC            nkptnrloc, nstsv, ltranblhloc )
#else
          !$OMP PARALLEL DO &
          !$OMP   PRIVATE( nj )
#endif /* _OPENACC */
          DO ik = 1, nkptnrloc

             ! Initialize values
             nj = 0
             is = 1

             ! Begin search algorithm (spinless)
             ! Note: it is desirable to keep the jbandidx indices in order
             !       (j will usually be consecutive)
#ifdef _OPENACC
             !$ACC LOOP SEQ
#endif /* _OPENACC */
             DO j = 1, nstsv

                ! Test whether the band contributes to transitions
                ! Note: When spinpol == .FALSE. the spinor_ud array doesn't exist
                IF( ltranblhloc(j,ik) ) THEN

                   ! Increment nj, then save band index j
                   ! Note: ATOMICs not needed since LOOP SEQ is in effect
                   nj = nj + 1
                   jbandidx(nj,is,ik) = j

                END IF ! ltranblhloc

             END DO ! nstsv
#ifdef _OPENACC
             !$ACC END LOOP
#endif /* _OPENACC */

             ! Assign njbands
             njbands(is,ik) = nj

          END DO ! nkptnrloc
#ifdef _OPENACC
          !$ACC END PARALLEL LOOP
#else
          !$OMP END PARALLEL DO
#endif /* _OPENACC */

       END IF ! spinpol

    END IF ! iq == 1 && ikloc == 1 && ispn = 1

    ! Find maximum value of nj and store it in njmax
#ifdef _OPENACC
    !$ACC PARALLEL PRESENT( njmax, njbands )
#endif /* _OPENACC */
    njmax = 0
#ifdef _OPENACC
    !$ACC LOOP GANG VECTOR COLLAPSE(2) REDUCTION(MAX:njmax)
#else
    !$OMP PARALLEL DO COLLAPSE(2) REDUCTION(MAX:njmax)
#endif /* _OPENACC */
    DO ik = 1, nkptnrloc
       DO is = 1, nspinor
          njmax = MAX( njmax, njbands(is,ik) )
       END DO ! nspinor
    END DO ! nkptnrloc
#ifdef _OPENACC
    !$ACC END LOOP
    !$ACC END PARALLEL
#else
    !$OMP END PARALLEL DO
#endif /* _OPENACC */

    ! Assign value for current ikloc and ispn
    !$ACC SERIAL COPYIN( ispn, ikloc )
    nj = njbands(ispn,ikloc)
    !$ACC END SERIAL

    ! Transfer result D->H
    !$ACC UPDATE HOST( jbandidx, njbands, njmax, nj )
    !$ACC WAIT

#if EBUG >= 1
    IF( ispn == 1 ) THEN
       spinproj = 1
    ELSE
       spinproj = -1
    END IF
    WRITE(*,*) 'countbands: ', spinproj, ' ikloc=', ikloc, ' ik=', iktable(ikloc), ' nj=', nj
    WRITE(*,*) jbandidx(:,ispn,ikloc)
#endif /* DEBUG */

    IF( nj > nband1 ) THEN
       WRITE(*,*) 'Warning[countbands]: nj ', nj, ' > nband1 ', nband1
    END IF

    RETURN
  END SUBROUTINE genmegqblh_countbands

!==============================================================================
! Kernel 1: Fill in bgntuju (or dptr_gntuju) and b1 arrays, and zero b2 array
!==============================================================================

  SUBROUTINE genmegqblh_fillbatch( wfsvmt1, ikloc, ispn )
    USE modmain, ONLY: zzero, natmtot, nspinor, nstsv, lmmaxapw, nufrmax
#ifdef _PACK_gntuju_
    USE mod_expigqr, ONLY: gntuju_packed, bmegqblh, idxtranblhloc, &
                           irowmap_wf1
#else
    USE mod_expigqr, ONLY: gntuju, bmegqblh, idxtranblhloc
#endif /* _PACK_gntuju_ */
    USE mod_addons, ONLY: ias2ic
    USE mod_addons_q, ONLY: sfacgq
    USE mod_nrkp, ONLY: spinor_ud

#ifdef _OPENACC
    USE openacc
#endif /* _OPENACC */

#ifdef _OPENMP
    USE omp_lib
#endif /* _OPENMP */

    USE mod_prof

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: ikloc, ispn
    COMPLEX(KIND=dz), DIMENSION(:,:,:,:,:), INTENT(IN) :: wfsvmt1

    ! Internal variables
    !INTEGER, PARAMETER :: nb = 64              ! Block size for ZGEMM batching
    INTEGER :: iblock                           ! Block index
    INTEGER(KIND=dl) :: ibatch                  ! Batch index
    INTEGER :: k1, k2, ki, nsize                ! Dummy variables for batching
    INTEGER(KIND=dl) :: i, j, ist1, ic, ig, ias ! Data access and/or loop indices
    INTEGER(KIND=dl) :: i1, i2, imt             ! Data access and/or loop indices
    INTEGER :: tid                              ! Thread ID
    INTEGER(KIND=dl) :: flop_b1

    ! Debugging variables
    LOGICAL :: li1, li2, limt, lki, list1, liasw, liass, lig, lispn, libatch
    INTEGER :: tmp, intmax
    INTEGER :: ist, iold
    INTEGER :: ltmp, lcond, lup, ldn, lrange, lpaired
    LOGICAL :: lcollapse4

#ifdef _CUDA_

    ! Call CUDA C++ kernel
    !CALL genmegqblh_fillbatch_cu_( ... )
    !CALL cudaMemset( b2, ... )

#else

    ! Batching by block size nb for idxhiband
    ! TODO: Re-enable if needed for larger systems
!  iblock = 0
!  DO k1 = 1, nj, nb
!     k2 = MIN( nj, k1+nb-1 )
!     nsize = k2 - k1 + 1
!     iblock = iblock + 1
    iblock = 1 ! Unblocked version

!--DEBUG
!    WRITE(*,*) 'entered genmegqblh_fillbatch, ikloc=', ikloc, ' ispn=',ispn
!--DEBUG
    
#ifdef _OPENACC
    ! Stub for multi-GPU support
    ! TODO: generalize for AMD GPUs
    !CALL acc_set_device_num( devnum, acc_device_nvidia )
#endif /* _OPENACC */

    CALL profstart( "Muffin-tin fill batchidx" )

    ! Fill in batchidx, the translation table for ibatch <-> {ig,ias,iblock}
    ! TODO: Change bounds for blocked algorithm
#ifdef _OPENACC
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) WAIT &
    !$ACC   COPYIN( iblock ) &
    !$ACC   PRESENT( natmtot, ngqiq, batchidx ) &
    !$ACC   PRIVATE( ibatch )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ibatch )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ! Blocked version
          ! ibatch = (iblock-1)*natmtot*ngq(iq) + (ig-1)*natmtot + ias
          ! batchidx(ias,ig,iblock) = ibatch

          ! Unblocked version (iblock = 1)
          ibatch = (ig-1)*natmtot + ias
          batchidx(ias,ig,iblock) = ibatch

       END DO ! ias
    END DO ! ig
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP
    !$ACC WAIT
    !$ACC UPDATE HOST( batchidx )
    !$ACC WAIT
#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

#if EBUG >= 2
    WRITE(*,*) 'fillbatch: ispn1=', ispn, ' nj=', nj, ' nbatch1=', nbatch1
#endif /* DEBUG */

    CALL profend( "Muffin-tin fill batchidx" )
    CALL profstart( "Muffin-tin zero b1" )

    ! Zero out b1 batch array
#if defined(_OPENACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) PRESENT( b1 )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED)
#endif /* _OPENACC || _OPENMP */
    DO ibatch = 1, SIZE(b1,3)
       DO ist1 = 1, SIZE(b1,2)
          DO imt = 1, SIZE(b1,1)
             b1(imt,ist1,ibatch) = zzero
          END DO ! imt
       END DO ! ist1
    END DO ! ibatch
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP
#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

    CALL profend( "Muffin-tin zero b1" )
    CALL profstart( "Muffin-tin fill b1" )

    ! Fill in b1 batch array
#if defined(_OPENACC) && defined(_PACK_gntuju_)
    !$ACC PARALLEL LOOP COLLAPSE(3) GANG &
    !$ACC   COPYIN( iblock, ikloc, ispn, irowmap_wf1 ) &
    !$ACC   PRIVATE( ic, i1, i2, ibatch, i, j, ist1, &
    !$ACC            li1, li2, limt, lki, list1, liasw, liass, lig, &
    !$ACC            lispn, libatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nj, nmtmax, lmmaxapw, nufrmax, &
    !$ACC            ias2ic, batchidx, jbandidx, idxtranblhloc, bmegqblh, &
    !$ACC            wfsvmt1, sfacgq, b1 )
#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
    !$ACC PARALLEL LOOP COLLAPSE(3) GANG &
    !$ACC   COPYIN( iblock, ikloc, ispn ) &
    !$ACC   PRIVATE( imt, ibatch, i, j, ist1, &
    !$ACC            li1, li2, limt, lki, list1, liasw, liass, lig, &
    !$ACC            lispn, libatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nj, nmtmax, lmmaxapw, nufrmax, &
    !$ACC            ias2ic, batchidx, jbandidx, idxtranblhloc, bmegqblh, &
    !$ACC            wfsvmt1, sfacgq, b1 )   
#elif defined(_OPENMP) && defined(_PACK_gntuju_)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) &
    !$OMP   PRIVATE( i1, i2, ic, ibatch, i, j, ist1, &
    !$OMP            li1, li2, limt, lki, list1, liasw, liass, lig, &
    !$OMP            lispn, libatch )
#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) &
    !$OMP   PRIVATE( imt, ibatch, i, j, ist1, &
    !$OMP            li1, li2, limt, lki, list1, liasw, liass, lig, &
    !$OMP            lispn, libatch )
#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot
          DO ki = 1, nj

#if defined(_OPENACC) && defined(_PACK_gntuju_)
             ic = ias2ic(ias)
             !$ACC LOOP VECTOR
             DO imt = 1, nmt(ic,ig)
                ! Use permutation vector after translation in gengntuju()
                i1 = irowmap_wf1(1,imt,ic,ig)
                i2 = irowmap_wf1(2,imt,ic,ig)
#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
             !$ACC LOOP COLLAPSE(2) VECTOR
             DO i2 = 1, nufrmax
                DO i1 = 1, lmmaxapw
                   imt = (i2-1)*lmmaxapw + i1
#elif defined(_OPENMP) && defined(_PACK_gntuju_)
             ic = ias2ic(ias)
             !$OMP SIMD
             DO imt = 1, nmt(ic,ig)
                ! Use permutation vector after translation in gengntuju()
                i1 = irowmap_wf1(1,imt,ic,ig)
                i2 = irowmap_wf1(2,imt,ic,ig)
#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
             !$OMP SIMD COLLAPSE(2)
             DO i2 = 1, nufrmax
                DO i1 = 1, lmmaxapw
                   imt = (i2-1)*lmmaxapw + i1
#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */

#if EBUG > 2
                ! Check array bounds
                ! i1, i2
                li1  = ( i1 >= LBOUND(wfsvmt1,1) ) .AND. &
                       ( i1 <= UBOUND(wfsvmt1,1) )
                li2  = ( i2 >= LBOUND(wfsvmt1,2) ) .AND. &
                       ( i2 <= UBOUND(wfsvmt1,2) )
                IF( .NOT. li1 ) THEN
                   WRITE(*,*) 'Error(fillbatch): i1 ', i1, &
                              ' reading wfsvmt1 out of bounds', &
                              LBOUND(wfsvmt1,1), UBOUND(wfsvmt1,1)
                   STOP
                END IF
                IF( .NOT. li2 ) THEN
                   WRITE(*,*) 'Error(fillbatch): i2 ', i2, &
                              ' reading wfsvmt1 out of bounds', &
                              LBOUND(wfsvmt1,2), UBOUND(wfsvmt1,2)
                   STOP
                END IF
                ! imt
                limt = ( imt >= LBOUND(b1,1) ) .AND. &
                       ( imt <= UBOUND(b1,1) )
                IF( .NOT. limt ) THEN
                   WRITE(*,*) 'Error(fillbatch): imt ', imt, &
                              ' writing b1 out of bounds', &
                              LBOUND(b1,1), UBOUND(b1,1)
                   STOP
                END IF
#endif /* DEBUG */

                ! Note: putting this here because both OpenMP and OpenACC
                !       don't like breaking up consecutive DO statements,
                !       even when commented out
                !DO j = 1, nsize ! Blocked version

                ibatch = batchidx(ias,ig,iblock)

                !j = k1 + ki - 1     ! Blocked version
                j = jbandidx( ki, ispn, ikloc ) ! Unblocked version

                i = idxtranblhloc( j, ikloc )
                ist1 = bmegqblh(1,i,ikloc)

#if EBUG > 2
                ! Check array bounds
                ! ki, ist1
                list1 = ( ist1 >= LBOUND(wfsvmt1,5) ) .AND. &
                        ( ist1 <= UBOUND(wfsvmt1,5) )
                lki   = ( ki >= LBOUND(b1,2) )        .AND. &
                        ( ki <= UBOUND(b1,2) )
                IF( .NOT. list1 ) THEN
                   WRITE(*,*) 'fillbatch: ist1 ', ist1, &
                              ' reading wfsvmt1 out of bounds', &
                              LBOUND(wfsvmt1,5), UBOUND(wfsvmt1,5)
                   STOP
                END IF
                IF( .NOT. lki ) THEN
                   WRITE(*,*) 'fillbatch: ki ', ki, &
                              ' writing b1 out of bounds', &
                              LBOUND(b1,2), UBOUND(b1,2)
                   STOP
                END IF

                ! ias
                liasw = ( ias >= LBOUND(wfsvmt1,3) ) .AND. &
                        ( ias <= UBOUND(wfsvmt1,3) )
                liass = ( ias >= LBOUND(sfacgq,2) ) .AND. &
                        ( ias <= UBOUND(sfacgq,2) )
                IF( .NOT. liasw ) THEN
                   WRITE(*,*) 'fillbatch: ias ', ias, &
                              ' reading wfsvmt1 out of bounds', &
                              LBOUND(wfsvmt1,3), UBOUND(wfsvmt1,3)
                   STOP
                END IF
                IF( .NOT. liass ) THEN
                   WRITE(*,*) 'fillbatch: ias ', ias, &
                              ' reading sfacgq out of bounds', &
                              LBOUND(sfacgq,2), UBOUND(sfacgq,2)
                   STOP
                END IF

                ! ig
                lig = ( ig >= LBOUND(sfacgq,1) ) .AND. &
                      ( ig <= UBOUND(sfacgq,1) )
                IF( .NOT. lig ) THEN
                   WRITE(*,*) 'fillbatch: ig ', ig, &
                              ' reading sfacgq out of bounds', &
                              LBOUND(sfacgq,1), UBOUND(sfacgq,1)
                   STOP
                END IF

                ! ispn
                lispn = ( ispn >= LBOUND(wfsvmt1,4) ) .AND. &
                        ( ispn <= UBOUND(wfsvmt1,4) )
                IF( .NOT. lispn ) THEN
                   WRITE(*,*) 'fillbatch: ispn ', ispn, &
                              ' reading wfsvmt1 out of bounds', &
                              LBOUND(wfsvmt1,4), UBOUND(wfsvmt1,4)
                   STOP
                END IF

                ! ibatch
                libatch = ( ibatch >= LBOUND(b1,3) ) .AND. ( ibatch <= UBOUND(b1,3) )
                IF( .NOT. libatch ) THEN
                   WRITE(*,*) 'fillbatch: ibatch ', ibatch, &
                              ' writing b1 out of bounds', &
                              LBOUND(b1,3), UBOUND(b1,3)
                   STOP
                END IF

#endif /* DEBUG */

                ! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
                b1(imt,ki,ibatch) = DCONJG( wfsvmt1(i1,i2,ias,ispn,ist1) * &
                                            sfacgq(ig,ias) )

#if defined(_OPENACC) && defined(_PACK_gntuju_)
             END DO ! imt
             !$ACC END LOOP
#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
                END DO ! i1
             END DO ! i2
             !$ACC END LOOP
#elif defined(_OPENMP) && defined(_PACK_gntuju_)
             END DO ! imt
             !$OMP END SIMD
#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
                END DO ! i1
             END DO ! i2
             !$OMP END SIMD
#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */

          END DO ! ki
       END DO ! ias
    END DO ! ig
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP

    ! fillbatch() contributes some GPU flops when filling in b1
    ! Note: flop_b1 is local to this subroutine
#ifdef _PACK_gntuju_
    flop_b1 = 0
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP   DEFAULT(SHARED) PRIVATE(ic) REDUCTION(+:flop_b1)
    DO ig = 1, ngqiq
       DO ias = 1, natmtot
          ic = ias2ic(ias)
          ! 7 FLOP (6 complex multiply + 1 CONJG) per ki and imt
          flop_b1 = flop_b1 + 7_dl * INT(nj,KIND=dl) &
                                   * INT(nmt(ic,ig),KIND=dl)
       END DO ! ias
    END DO ! ig
    !$OMP END PARALLEL DO
#else
    ! 7 FLOP (6 complex multiply + 1 CONJG) per ig, ias, ki, i1, and i2
    flop_b1 = 7_dl * INT(ngqiq,KIND=dl) * INT(natmtot,KIND=dl) &
                   * INT(nj,KIND=dl) &
                   * INT(nufrmax,KIND=dl) * INT(lmmaxapw,KIND=dl)
#endif /* _PACK_gntuju_ */
    ! Note: flop_fillbatch is declared in mod_gpu
    flop_fillbatch = flop_fillbatch + flop_b1

#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

    CALL profend( "Muffin-tin fill b1" )

    ! Note: no need to zero out b2 batch array,
    ! since both ZGEMM() and magmablas_zgemm_batched() zeroes it out anyway

    CALL profstart( "Muffin-tin fill pointers" )

#if defined(_OPENACC) && defined(_PACK_gntuju_)

    ! Fill in array of device pointers
    ! TODO: Change bounds for blocked algorithm

    !$ACC DATA PRESENT( gntuju_packed, b1, b2 )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC   COPYIN( iblock ) PRIVATE( ic, ibatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nj, nmtmax, batchidx, ias2ic, &
    !$ACC            dptr_gntuju, dptr_b1, dptr_b2 )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)
          ic = ias2ic(ias)

          dptr_gntuju(ibatch) = C_LOC( gntuju_packed(1,1,ic,ig) )
          dptr_b1(ibatch) = C_LOC( b1(1,1,ibatch) )
          dptr_b2(ibatch) = C_LOC( b2(1,1,ibatch) )

       END DO ! ias
    END DO ! ig
    !$ACC END PARALLEL LOOP

    ! gntuju_packed, b1, b2
    !$ACC END DATA

#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
    ! Fill in array of device pointers
    ! TODO: Change bounds for blocked algorithm

    !$ACC DATA PRESENT( gntuju, b1, b2 )
    !$ACC HOST_DATA USE_DEVICE( gntuju, b1, b2 )
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
    !$ACC   COPYIN( iblock ) PRIVATE( ic, ibatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nj, nmtmax, batchidx, ias2ic, &
    !$ACC            dptr_gntuju, dptr_b1, dptr_b2 )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)
          ic = ias2ic(ias)

          dptr_gntuju(ibatch) = C_LOC( gntuju(1,1,ic,ig) )
          dptr_b1(ibatch) = C_LOC( b1(1,1,ibatch) )
          dptr_b2(ibatch) = C_LOC( b2(1,1,ibatch) )

       END DO ! ias
    END DO ! ig
    !$ACC END PARALLEL LOOP

    ! gntuju, b1, b2
    !$ACC END HOST_DATA
    !$ACC END DATA

#elif defined(_OPENMP) && defined(_PACK_gntuju_)
    ! Copy gntuju_packed to bgntuju
    ! TODO: Change bounds for blocked algorithm

    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP   PRIVATE( ic, ibatch )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)
          ic = ias2ic(ias)

          !$OMP SIMD COLLAPSE(2)
          DO i2 = 1, nmtmax
             DO i1 = 1, nmtmax
                bgntuju(i1,i2,ibatch) = gntuju_packed(i1,i2,ic,ig)
             END DO ! i1
          END DO ! i2
          !$OMP END SIMD

       END DO ! ias
    END DO ! ig
    !$OMP END PARALLEL DO

#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
    ! Copy gntuju to bgntuju

    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP   PRIVATE( ic, ibatch )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)
          ic = ias2ic(ias)

          !$OMP SIMD COLLAPSE(2)
          DO i2 = 1, nmtmax
             DO i1 = 1, nmtmax
                bgntuju(i1,i2,ibatch) = gntuju(i1,i2,ic,ig)
             END DO ! i1
          END DO ! i2
          !$OMP END SIMD

       END DO ! ias
    END DO ! ig
    !$OMP END PARALLEL DO

#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */

    CALL profend( "Muffin-tin fill pointers" )

!--DEBUG
!    !$ACC PARALLEL LOOP COLLAPSE(2) &
!    !$ACC   COPY( iblock ) &
!    !$ACC   PRIVATE( ig, ias, ic, ibatch ) &
!    !$ACC   PRESENT( natmtot, ngqiq, &
!    !$ACC            batchidx, ias2ic, &
!    !$ACC            gntuju, bgntuju )
!    DO ig = 1, ngqiq
!       DO ias = 1, natmtot
!          ibatch = batchidx(ias,ig,iblock)
!          ic = ias2ic(ias)
!          bgntuju(:,:,ibatch) = gntuju(:,:,ic,ig)
!       END DO ! ias
!    END DO ! ig
!    !$ACC END PARALLEL LOOP
!--DEBUG

!  END DO ! k1

#endif /* _CUDA_ */

!--DEBUG
!    WRITE(*,*) 'exiting genmegqblh_fillbatch'
!--DEBUG
    
    RETURN
  END SUBROUTINE genmegqblh_fillbatch

!==============================================================================
! Kernel 2: Perform batched ZGEMM b2(:,:) = bgntuju(:,:) x b1(:,:)
!==============================================================================

  SUBROUTINE genmegqblh_batchzgemm(nbatch)

    !USE modmain, ONLY: zzero, zone
    USE ISO_C_BINDING, ONLY: C_INT

#ifdef _MAGMA_
    USE mod_magma
    USE mod_gpu
#elif defined(_CUBLAS_)
#else
    USE mod_gpu, ONLY: zgemm_strided_batched_omp
#endif /* _MAGMA_ */

    IMPLICIT NONE

    ! Input argument
    INTEGER, INTENT(IN) :: nbatch

    ! Internal variables
    COMPLEX(KIND=dz), PARAMETER :: alpha = (1._dd,0._dd)
    COMPLEX(KIND=dz), PARAMETER :: beta  = (0._dd,0._dd)
    INTEGER :: m, n, k, lda, ldb, ldc
    INTEGER :: ncolA, ncolB, ncolC, stA, stB, stC, ibatch
    INTEGER, DIMENSION(nbatch) :: d_m, d_n, d_k, d_lda, d_ldb, d_ldc

  !-2a-------------------------------------------------------------------------
    IF( usemagma ) THEN
  !----------------------------------------------------------------------------
  ! MAGMA version: use magmablas_zgemm_batched()

#ifdef _PACK_gntuju_
       ! gntuju is packed
       m = npackdim
       n = nj
       k = npackdim
       lda = SIZE(gntuju_packed,1)
       ldb = SIZE(b1,1)
       ldc = SIZE(b2,1)
#else
       ! gntuju isn't packed
       m = nmtmax
       n = nj
       k = nmtmax
       lda = SIZE(gntuju,1)
       ldb = SIZE(b1,1)
       ldc = SIZE(b2,1)
#endif /*_PACK_gntuju_ */
          
#if EBUG > 0
       WRITE(*,*) 'batchzgemm: using zgemm_batched'
       WRITE(*,*) 'batchzgemm: nbatch=', nbatch, &
                  ' m=', m, ' n=', n, ' k=', k, &
                  ' lda=', lda, ' ldb=', ldb, ' ldc=', ldc
#endif /* DEBUG */

       !$ACC DATA PRESENT( dptr_gntuju, dptr_b1, dptr_b2 )
       !$ACC HOST_DATA USE_DEVICE( dptr_gntuju, dptr_b1, dptr_b2 )

       ! Perform batched ZGEMM on device using MAGMA (pointer mode)
       CALL zgemm_batched_gpu_acc_magma_ptr( 'N', 'N', &
                                             m, n, k, &
                                             alpha, C_LOC(dptr_gntuju), lda, &
                                                    C_LOC(dptr_b1),     ldb, &
                                             beta,  C_LOC(dptr_b2),     ldc, &
                                             nbatch )
       ! Note: mod_gpu::zgemm_batched_gpu_acc_magma_ptr() already calls
       !       magma_queue_sync(), so no need to explicitly add it here

       !dptr_gntuju, dptr_b1, dptr_b2
       !$ACC END HOST_DATA
       !$ACC END DATA

       !$ACC WAIT

       ! FLOP formula taken from MAGMA testing/flops.h
       ! Note: flop_batchzgemm is declared in mod_gpu
       flop_batchzgemm = flop_batchzgemm + 8_dl * INT(m,KIND=dl) &
                                                * INT(n,KIND=dl) &
                                                * INT(k,KIND=dl) &
                                                * INT(nbatch,KIND=dl)

  !-2b-------------------------------------------------------------------------
    !ELSE IF( usecublas )
  !----------------------------------------------------------------------------

       ! Perform batched ZGEMM on device using cuBLAS
       !CALL cublasZgemmBatched( blashandle, CUBLAS_OP_N, CUBLAS_OP_N, ... )

  !-2c-------------------------------------------------------------------------
    ELSE ! Fall back to CPU only using OpenMP
  !----------------------------------------------------------------------------
  ! OpenMP version: batched strided ZGEMM

!--DEBUG
       !$ACC UPDATE HOST( b1, bgntuju )
!--DEBUG

#ifdef _PACK_gntuju_
       ! gntuju is packed
       m = npackdim
       n = nj
       k = npackdim
#else
       ! gntuju isn't packed
       m = nmtmax
       n = nj
       k = nmtmax
#endif /*_PACK_gntuju_ */

       lda = SIZE(bgntuju,1)
       ldb = SIZE(b1,1)
       ldc = SIZE(b2,1)

#if EBUG > 0
       WRITE(*,*) 'batchzgemm: nbatch=', nbatch, ' m=', m, ' n=', n, ' k=', k
#endif /* DEBUG */

       ncolA = SIZE(bgntuju,2)
       ncolB = SIZE(b1,2)
       ncolC = SIZE(b2,2)

       ! Set up strides
       stA = lda * ncolA
       stB = ldb * ncolB
       stC = ldc * ncolC

       ! Perform batched ZGEMM on CPU using OpenMP parallel do
       ! b2(1:nmt,1:nj) = bgntuju(1:nmt,1:nmt) x b1(1:nmt,1:nj)
       CALL zgemm_strided_batched_omp( 'N', 'N', m, n, k, &
                                       alpha, bgntuju, lda, stA, &
                                              b1,      ldb, stB, &
                                       beta,  b2,      ldc, stC, &
                                       nbatch )

  !----------------------------------------------------------------------------
    END IF ! CPU/GPU method
  !----------------------------------------------------------------------------

    RETURN
  END SUBROUTINE genmegqblh_batchzgemm

!==============================================================================
! Kernel 3: Save results to wftmp1mt and transfer back to CPU (for now)
!==============================================================================

  SUBROUTINE genmegqblh_fillresult( wftmp1mt )
    USE modmain, ONLY: natmtot, nstsv, lmmaxapw, nufrmax, zzero
#ifdef _PACK_gntuju_
    USE mod_expigqr, ONLY: irownz, irowmap_res
#endif /* _PACK_gntuju_ */

#ifdef _OPENACC
    USE openacc
#endif /* _OPENACC */

#ifdef _OPENMP
     USE omp_lib
#endif /* _OPENMP */
    
    IMPLICIT NONE

    ! Argument
    COMPLEX(KIND=dz), DIMENSION(:,:,:,:), INTENT(INOUT) :: wftmp1mt
    
    ! Internal variables
    INTEGER(KIND=dl) :: j, ist, i1, imt, iblock, ibatch, ias, ic, ig, tid

    ! Debugging variables
    LOGICAL :: li1w, li1b, li2, lj, list1, liasw, lig, libatch

#ifdef _CUDA_

    ! Zero out wftmp1mt()
    !CALL cudaMemset( ... )

    ! Call CUDA C++ kernel
    !CALL genmegqblh_fillresult_cu_( ... )

    ! Transfer data from device to host
    !CALL cudaMemcpy( wfsvmt1, d_wftmp1mt, cudaMemcpyDeviceToHost )

    ! Clean up device pointers
    !CALL cudaFree( d_wfsvmt1 )
    !CALL cudaFree( d_sfacgq )
    !CALL cudaFree( d_gntuju )
    !CALL cudaFree( d_bgntuju )
    !CALL cudaFree( d_b1 )
    !CALL cudaFree( d_b2 )

#else

!--DEBUG
!    WRITE(*,*) 'entered genmegqblh_fillresult
!--DEBUG

    ! Blocked version
!  DO iblock = 1, nblock
!     k1 = (j-1)*nb + 1
!     IF( iblock == nblock ) THEN
!        k2 = idxhiband
!     ELSE
!        k2 = j*nb
!     END IF

    iblock = 1 ! Unblocked version

    ! Stub for multi-GPU support
    ! TODO: generalize for AMD GPUs
    !CALL acc_set_device_num( devnum, acc_device_nvidia )

#if defined(_OPENACC) && defined(_PACK_gntuju_)
    ! Zero out and fill in wftmp1mt on device

    !$ACC PARALLEL LOOP COLLAPSE(2) GANG &
    !$ACC   PRESENT( ngqiq, natmtot, lmmaxapw, nufrmax, nj, ias2ic, &
    !$ACC            batchidx, b2, wftmp1mt ) &
    !$ACC   COPYIN( iblock, irowmap_res ) &
    !$ACC   PRIVATE( ibatch, ic, &
    !$ACC            li1w, li1b, lj, list1, liasw, lig, libatch )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ic = ias2ic(ias)
          ibatch = batchidx(ias,ig,iblock)

          !$ACC LOOP COLLAPSE(2) VECTOR
          DO j = 1, nj
             DO imt = 1, nmt(ic,ig)

#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
    ! Fill in wftmp1mt on device

    !$ACC PARALLEL LOOP COLLAPSE(2) GANG &
    !$ACC   PRESENT( ngqiq, natmtot, nmtmax, nj, ias2ic, &
    !$ACC            batchidx, b2, wftmp1mt ) &
    !$ACC   COPYIN( iblock ) &
    !$ACC   PRIVATE( ibatch, &
    !$ACC            li1w, li1b, lj, list1, liasw, lig, libatch )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)

          !$ACC LOOP COLLAPSE(2) VECTOR
          DO j = 1, nj
             DO imt = 1, nmtmax

#elif defined(_OPENMP) && defined(_PACK_gntuju_)
    ! Zero out wftmp1mt and copy b2 to wftmp1mt (with unpacking)

    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE( tid, ibatch, ic, &
    !$OMP            li1w, li1b, lj, list1, liasw, lig, libatch )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ic = ias2ic(ias)
          ibatch = batchidx(ias,ig,iblock)

#if EBUG > 2
                tid = omp_get_thread_num()
                WRITE(*,*) 'genmegqblh_fillresult: tid=', tid, &
                           ' ias=', ias, ' ig=', ig, &
                           ' ibatch=', ibatch
#endif /* DEBUG */

          !$OMP SIMD COLLAPSE(2)
          DO j = 1, nj
             DO imt = 1, nmt(ic,ig)

#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
    ! Copy b2 to wftmp1mt

    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE( tid, ibatch, &
    !$OMP            li1w, li1b, lj, list1, liasw, lig, libatch )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)

#if EBUG > 2
          tid = omp_get_thread_num()
          WRITE(*,*) 'genmegqblh_fillresult: tid=', tid, &
                     ' ias=', ias, ' ig=', ig, &
                     ' ibatch=', ibatch
#endif /* DEBUG */

          !$OMP SIMD COLLAPSE(2)
          DO j = 1, nj
             DO imt = 1, nmtmax

#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */

#if EBUG > 2
                ! Check array bounds
                ! i1, imt
                li1w = ( imt >= LBOUND(wftmp1mt,1) ) .AND. &
                       ( imt <= UBOUND(wftmp1mt,1) )
                li1b = ( imt >= LBOUND(b2,1) )       .AND. &
                       ( imt <= UBOUND(b2,1) )
                IF( .NOT. li1w ) THEN
                   WRITE(*,*) 'fillresult: imt ', imt, &
                              ' writing wftmp1mt out of bounds', &
                              LBOUND(wftmp1mt,1), UBOUND(wftmp1mt,1)
                   STOP
                END IF
                IF( .NOT. li1b ) THEN
                   WRITE(*,*) 'fillresult: imt ', imt, &
                              ' reading b2 out of bounds', &
                              LBOUND(b2,1), UBOUND(b2,1)
                   STOP
                END IF

                ! j, ist1
                list1 = ( j >= LBOUND(wftmp1mt,2) ) .AND. &
                        ( j <= UBOUND(wftmp1mt,2) )
                lj    = ( j >= LBOUND(b2,2) )       .AND. &
                        ( j <= UBOUND(b2,2) )
                IF( .NOT. list1 ) THEN
                   WRITE(*,*) 'fillresult: j ', j, &
                              ' writing wftmp1mt out of bounds', &
                              LBOUND(wftmp1mt,2), UBOUND(wftmp1mt,2)
                   STOP
                END IF
                IF( .NOT. lj ) THEN
                   WRITE(*,*) 'fillresult: j ', j, &
                              ' reading b2 out of bounds', &
                              LBOUND(b2,2), UBOUND(b2,2)

                   STOP
                END IF

                ! ias
                liasw = ( ias >= LBOUND(wftmp1mt,3) ) .AND. &
                        ( ias <= UBOUND(wftmp1mt,3) )
                IF( .NOT. liasw ) THEN
                   WRITE(*,*) 'fillresult: ias ', ias, &
                              ' writing wftmp1mt out of bounds', &
                              LBOUND(wftmp1mt,3), UBOUND(wftmp1mt,3)
                   STOP
                END IF

                ! ig
                lig = ( ig >= LBOUND(wftmp1mt,4) ) .AND. &
                      ( ig <= UBOUND(wftmp1mt,4) )
                IF( .NOT. lig ) THEN
                   WRITE(*,*) 'fillresult: ig ', ig, &
                              ' writing wftmp1mt out of bounds', &
                              LBOUND(wftmp1mt,4), UBOUND(wftmp1mt,4)
                   STOP
                END IF

                ! ibatch
                libatch = ( ibatch >= LBOUND(b2,3) ) .AND. &
                          ( ibatch <= UBOUND(b2,3) )
                IF( .NOT. libatch ) THEN
                   WRITE(*,*) 'fillresult: ibatch ', ibatch, &
                        ' reading b2 out of bounds', &
                        LBOUND(b2,3), UBOUND(b2,3)
                   STOP
                END IF
#endif /* DEBUG */

                wftmp1mt(imt,j,ias,ig) = b2(imt,j,ibatch)

#if defined(_OPENACC) && defined(_PACK_gntuju_)
             END DO ! imt
          END DO ! j
          !$ACC END LOOP

          ! Note: this can be a zero-trip loop
          !$ACC LOOP COLLAPSE(2) VECTOR
          DO j = 1, nj
             DO imt = nmt(ic,ig) + 1, nmtmax
                wftmp1mt(imt,j,ias,ig) = zzero
             END DO ! imt
          END DO ! j
          !$ACC END LOOP

#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
             END DO ! i1
          END DO ! j
          !$ACC END LOOP

#elif defined(_OPENMP) && defined(_PACK_gntuju_)
             END DO ! imt
          END DO ! j
          !$OMP END SIMD

          ! Note: this can be a zero-trip loop
          !$OMP SIMD COLLAPSE(2)
          DO j = 1, nj
             DO imt = nmt(ic,ig) + 1, nmtmax
                wftmp1mt(imt,j,ias,ig) = zzero
             END DO ! imt
          END DO ! j
          !$OMP END SIMD

#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
             END DO ! i1
          END DO ! j
          !$OMP END SIMD
#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */

        END DO ! ias
     END DO ! ig
#ifdef _OPENACC
     !$ACC END PARALLEL LOOP
     !$ACC WAIT
#elif defined(_OPENMP)
     !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

#endif /* _CUDA_ */

    RETURN
  END SUBROUTINE genmegqblh_fillresult

!==============================================================================

END MODULE mod_genmegqblh_gpu
