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
  
  ! Table of spin-up/dn states (replaces l1 check)
  ! Dimension is nstsv, but only the first nstspin elements will be used
  INTEGER, DIMENSION(:), ALLOCATABLE :: spinstidx

  ! Number of 2nd-variational states per spin
  INTEGER :: nstspin

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
    USE modmain, ONLY: natmtot, lmmaxapw, nufrmax
    USE mod_addons_q, ONLY: ngq
    USE mod_expigqr, ONLY: nbandblhloc
    IMPLICIT NONE

    ! Arguments
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

#ifdef _OPENACC

    ! Copy/allocate constants to device
    ! Note: natmtot, nspinor, and nstsv are declared in modmain
    !$ACC ENTER DATA CREATE( nstspin ) &
    !$ACC            COPYIN( natmtot, nstsv, nspinor, lmmaxapw, nufrmax, &
    !$ACC                    nblock1, nbatch1, nband1, nmtmax, ngqiq )

#elif defined(_CUDA_)

    ! Allocate constants on device
    !CALL cudaMalloc( d_natmtot, ... )
    !CALL cudaMalloc( d_nstsv, ... )
    !CALL cudaMalloc( d_nspinor, ... )
    !CALL cudaMalloc( d_lmmaxapw, ... )
    !CALL cudaMalloc( d_nufrmax, ... )
    !CALL cudaMalloc( d_nstspin, ... )
    !CALL cudaMalloc( d_nblock1, ... )
    !CALL cudaMalloc( d_nbatch1, ... )
    !CALL cudaMalloc( d_nband1, ... )
    !CALL cudaMalloc( d_ngqiq, ... )

    ! Copy constants H->D
    !CALL cudaMemcpy( ... )

#endif /* _OPENACC || _CUDA_ */

    RETURN
  END SUBROUTINE genmegqblh_allocmodvar_const

!==============================================================================
! Cleans up constant module variables on device

  SUBROUTINE genmegqblh_freemodvar_const
    IMPLICIT NONE

#ifdef _OPENACC

    !$ACC EXIT DATA DELETE ( natmtot, nspinor, nstsv, lmmaxapw, nufrmax, &
    !$ACC                    nstspin, &
    !$ACC                    nblock1, nbatch1, nband1, nmtmax, ngqiq )

#elif defined(_CUDA_)

    !CALL cudaFree( ... )

#endif /* _OPENACC || _CUDA_ */

    RETURN
  END SUBROUTINE genmegqblh_freemodvar_const

!==============================================================================
! Allocates module variables on CPU and device (GPU)
! that are spin-dependent, i.e., related to genmegqblh_countspin() kernel

  SUBROUTINE genmegqblh_allocmodvar_spin()
    USE modmain, ONLY: nstsv
    IMPLICIT NONE

     ! Allocate array for table of states per spin projection on CPU memory
     ALLOCATE( spinstidx(nstsv) )

#ifdef _OPENACC

     ! Allocate array for table of states per spin projection on device
     !$ACC ENTER DATA CREATE( spinstidx )

#elif defined(_CUDA_)

     ! Allocate array for table of states per spin projection on device
     !CALL cudaMalloc( spinstidx, ... )

     ! Zero the array
     !CALL cudaMemset( ... )

#endif /* _OPENACC || _CUDA_ */

     RETURN
  END SUBROUTINE genmegqblh_allocmodvar_spin

!==============================================================================
! Cleans up module variables on CPU and device (GPU)
! that are spin-dependent, i.e., related to genmegqblh_countspin() kernel

  SUBROUTINE genmegqblh_freemodvar_spin()
    IMPLICIT NONE

#ifdef _OPENACC

    ! Clean up device
    !$ACC EXIT DATA DELETE( spinstidx )

#elif defined(_CUDA_)

    ! Clean up device
    !CALL cudaFree( spinstidx, ... )

#endif /* _OPENACC || _CUDA_ */

    ! Clean up CPU memory
    DEALLOCATE( spinstidx )

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
    ALLOCATE( wftmp1mt( lmmaxapw*nufrmax, nstspin, natmtot, ngqiq ) )

    ! Allocate batch arrays for the temporary matrices on CPU
    ! Note: we don't use bgntuju in the OpenACC implementation
#ifdef _PACK_gntuju_
    ALLOCATE( b1( npackdim, nstspin, nbatch1 ) )
    ALLOCATE( b2( npackdim, nstspin, nbatch1 ) )
#else
    ALLOCATE( b1( nmtmax, nstspin, nbatch1 ) )
    ALLOCATE( b2( nmtmax, nstspin, nbatch1 ) )
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
!           and returns a list of such states as spinstidx
!==============================================================================

  SUBROUTINE genmegqblh_countspin( ispn, ikloc, ik )

    USE modmain, ONLY: nstsv, spinpol
    USE mod_nrkp, ONLY: spinor_ud
    USE mod_expigqr, ONLY: bmegqblh, nbandblhloc, &
                           idxtranblhloc, ltranblhloc, ntranblhloc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: ispn, ikloc, ik

    ! Internal variables
    INTEGER :: iband, spinproj
    LOGICAL :: lpaired, lspinproj

    ! Zero out spinstidx
#ifdef _OPENACC
    !$ACC PARALLEL LOOP PRESENT( nstsv, spinstidx )
#else
    !$OMP PARALLEL DO
#endif /* _OPENACC */
    DO iband = 1, nstsv
       spinstidx(iband) = 0
    END DO ! nstsv
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP
#else
    !$OMP END PARALLEL DO
#endif /* _OPENACC */

#ifdef _OPENACC
    !$ACC KERNELS COPYIN( ikloc, ik ) &
    !$ACC   PRESENT( spinstidx, nstspin, &
    !$ACC            nstsv, spinor_ud, ltranblhloc )
#else
    !$OMP MASTER
#endif /* _OPENACC */

    ! Initialize value
    nstspin = 0

    ! Begin search algorithm
    ! TODO: Parallelize
#ifdef _OPENACC
    !$ACC LOOP SEQ PRIVATE( iband, lpaired, lspinproj )
#endif /* _OPENACC */
    DO iband = 1, nstsv

       ! Test whether the band contributes to transitions
       lpaired = ltranblhloc(iband,ikloc)

       IF( lpaired ) THEN

          IF( spinpol ) THEN
             ! Test the condition (Are we counting spin up or spin down states?)
             lspinproj = ( spinor_ud(ispn,iband,ik) /= 0 )
          ELSE
             ! When spinpol == .FALSE. the array spinor_ud doesn't even exist
             lspinproj = .TRUE.
          END IF

          IF( lspinproj ) THEN

             ! Increment nstspin, then save band index iband
             ! Note: ATOMICs not needed since LOOP SEQ is in effect
             nstspin = nstspin + 1
             spinstidx(nstspin) = iband

          END IF ! lspinproj

       END IF ! lpaired

    END DO ! nstsv
#ifdef _OPENACC
    !$ACC END LOOP
    !$ACC END KERNELS
#else
    !$OMP END MASTER
#endif /* _OPENACC */

    ! Transfer result D->H
    !$ACC UPDATE HOST( spinstidx, nstspin )
    !$ACC WAIT

#if EBUG >= 1
    IF( ispn == 1 ) THEN
       spinproj = 1
    ELSE
       spinproj = -1
    END IF
    WRITE(*,*) 'countspin: ', spinproj, ' ikloc=', ikloc, ' ik=', ik, ' nstspin=', nstspin
    WRITE(*,*) spinstidx
#endif /* DEBUG */

    IF( nstspin > nband1 ) THEN
       WRITE(*,*) 'Warning[countspin]: nstspin ', nstspin, ' > nband1 ', nband1
    END IF

    RETURN
  END SUBROUTINE genmegqblh_countspin

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

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: ikloc, ispn
    COMPLEX(KIND=dz), DIMENSION(:,:,:,:,:), INTENT(IN) :: wfsvmt1

    ! Internal variables
    !INTEGER, PARAMETER :: nb = 64         ! Block size for ZGEMM batching
    INTEGER :: iblock                      ! Block index
    INTEGER(KIND=dl) :: ibatch                      ! Batch index
    INTEGER :: k1, k2, ki, nsize           ! Dummy variables for batching
    INTEGER(KIND=dl) :: iband, i, ist1, ic, ig, ias ! Data access and/or loop indices
    INTEGER(KIND=dl) :: i1, i2, imt                 ! Data access and/or loop indices
    INTEGER :: tid                         ! Thread ID

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
!  DO k1 = 1, nstspin, nb
!     k2 = MIN( nstspin, k1+nb-1 )
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

    ! Fill in batchidx, the translation table for ibatch <-> {ig,ias,iblock}
#ifdef _OPENACC
    !$ACC PARALLEL LOOP COLLAPSE(2) WAIT &
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
    WRITE(*,*) 'fillbatch: ispn1=', ispn, ' nstspin=', nstspin, &
               ' nbatch1=', nbatch1
#endif /* DEBUG */

    ! Zero out b1 batch array
#if defined(_OPENACC)
    !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT( b1 )
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

    ! Fill in b1 batch array
#if defined(_OPENACC) && defined(_PACK_gntuju_)
    !$ACC PARALLEL LOOP COLLAPSE(3) GANG &
    !$ACC   COPYIN( iblock, ikloc, ispn, irowmap_wf1 ) &
    !$ACC   PRIVATE( ic, i1, i2, ibatch, iband, i, ist1, &
    !$ACC            li1, li2, limt, lki, list1, liasw, liass, lig, &
    !$ACC            lispn, libatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nstspin, nmtmax, lmmaxapw, nufrmax, &
    !$ACC            ias2ic, batchidx, spinstidx, idxtranblhloc, bmegqblh, &
    !$ACC            wfsvmt1, sfacgq, b1 )
#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
    !$ACC PARALLEL LOOP COLLAPSE(3) GANG &
    !$ACC   COPYIN( iblock, ikloc, ispn ) &
    !$ACC   PRIVATE( imt, ibatch, iband, i, ist1, &
    !$ACC            li1, li2, limt, lki, list1, liasw, liass, lig, &
    !$ACC            lispn, libatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nstspin, nmtmax, lmmaxapw, nufrmax, &
    !$ACC            ias2ic, batchidx, spinstidx, idxtranblhloc, bmegqblh, &
    !$ACC            wfsvmt1, sfacgq, b1 )   
#elif defined(_OPENMP) && defined(_PACK_gntuju_)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) &
    !$OMP   PRIVATE( i1, i2, ic, ibatch, iband, i, ist1, &
    !$OMP            li1, li2, limt, lki, list1, liasw, liass, lig, &
    !$OMP            lispn, libatch )
#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) &
    !$OMP   PRIVATE( imt, ibatch, iband, i, ist1, &
    !$OMP            li1, li2, limt, lki, list1, liasw, liass, lig, &
    !$OMP            lispn, libatch )
#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot
          DO ki = 1, nstspin

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
                !DO ki = 1, nsize ! Blocked version

                ibatch = batchidx(ias,ig,iblock)

                !iband = k1 + ki - 1     ! Blocked version
                iband = spinstidx( ki ) ! Unblocked version

                i = idxtranblhloc( iband, ikloc )
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
#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

    ! Zero b2 batch array
#ifdef _OPENACC
    !$ACC PARALLEL LOOP COLLAPSE(3) GANG VECTOR
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED)
#endif /* _OPENACC || _OPENMP */
    DO ibatch = 1, SIZE(b2,3)
       DO ist1 = 1, SIZE(b1,2)
          DO imt = 1, SIZE(b1,1)
             b2(imt,ist1,ibatch) = zzero
          END DO ! imt
       END DO ! ist1
    END DO ! ibatch
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP
#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

#if defined(_OPENACC) && defined(_PACK_gntuju_)
    ! Fill in array of device pointers

    !$ACC DATA PRESENT( gntuju_packed, b1, b2 )
    !$ACC HOST_DATA USE_DEVICE( gntuju_packed, b1, b2 )
    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC   COPYIN( iblock ) PRIVATE( ic, ibatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nstspin, nmtmax, batchidx, ias2ic, &
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
    !$ACC END HOST_DATA
    !$ACC END DATA

#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
    ! Fill in array of device pointers

    !$ACC DATA PRESENT( gntuju, b1, b2 )
    !$ACC HOST_DATA USE_DEVICE( gntuju, b1, b2 )
    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC   COPYIN( iblock ) PRIVATE( ic, ibatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nstspin, nmtmax, batchidx, ias2ic, &
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
       n = nstspin
       k = npackdim
       lda = SIZE(gntuju_packed,1)
       ldb = SIZE(b1,1)
       ldc = SIZE(b2,1)
#else
       ! gntuju isn't packed
       m = nmtmax
       n = nstspin
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

       ! Perform batched ZGEMM on device using MAGMA (pointer mode)
       CALL zgemm_batched_gpu_acc_magma_ptr( 'N', 'N', &
                                             m, n, k, &
                                             alpha, dptr_gntuju, lda, &
                                                    dptr_b1,     ldb, &
                                             beta,  dptr_b2,     ldc, &
                                             nbatch )

#if EBUG > 0
       WRITE(*,*) 'batchzgemm: using zgemm_batched'
       WRITE(*,*) 'batchzgemm: nbatch=', nbatch, &
                  ' m=', m, ' n=', n, ' k=', k, &
                  ' lda=', lda, ' ldb=', ldb, ' ldc=', ldc
#endif /* DEBUG */

       ! Note: PARAMETERs don't need to be COPYIN-ed to device

       ! Perform batched ZGEMM on device using MAGMA (pointer mode)
       CALL zgemm_batched_gpu_acc_magma_ptr( 'N', 'N', &
                                             m, n, k, &
                                             alpha, dptr_gntuju, lda, &
                                                    dptr_b1,     ldb, &
                                             beta,  dptr_b2,     ldc, &
                                             nbatch )
       !$ACC WAIT

#ifdef _MAGMA_
       ! Synchronize with device
       CALL magma_queue_sync( queue )
#endif /* _MAGMA_ */

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
       n = nstspin
       k = npackdim
#else
       ! gntuju isn't packed
       m = nmtmax
       n = nstspin
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
       ! b2(1:nmt,1:nstspin) = bgntuju(1:nmt,1:nmt) x b1(1:nmt,1:nstspin)
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
    INTEGER(KIND=dl) :: ki, ist, i1, imt, iblock, ibatch, ias, ic, ig, tid

    ! Debugging variables
    LOGICAL :: li1w, li1b, li2, lki, list1, liasw, lig, libatch

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
!     k1 = (ki-1)*nb + 1
!     IF( iblock == nblock ) THEN
!        k2 = idxhiband
!     ELSE
!        k2 = ki*nb
!     END IF

    iblock = 1 ! Unblocked version

    ! Stub for multi-GPU support
    ! TODO: generalize for AMD GPUs
    !CALL acc_set_device_num( devnum, acc_device_nvidia )

#if defined(_OPENACC) && defined(_PACK_gntuju_)
    ! Zero out and fill in wftmp1mt on device (with unpacking)

    !$ACC PARALLEL LOOP COLLAPSE(2) GANG &
    !$ACC   PRESENT( ngqiq, natmtot, lmmaxapw, nufrmax, nstspin, ias2ic, &
    !$ACC            batchidx, b2, wftmp1mt ) &
    !$ACC   COPYIN( iblock, irowmap_res ) &
    !$ACC   PRIVATE( ibatch, ist, ic, i1, &
    !$ACC            li1w, li1b, lki, list1, liasw, lig, libatch )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ic = ias2ic(ias)
          ibatch = batchidx(ias,ig,iblock)

          !$ACC LOOP COLLAPSE(2) VECTOR
          DO ki = 1, nstspin
             DO imt = 1, lmmaxapw*nufrmax
                i1 = irowmap_res(imt,ic,ig)
                IF( i1 == 0 ) THEN
                   wftmp1mt(imt,ki,ias,ig) = zzero
                ELSE

#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
    ! Fill in wftmp1mt on device

    !$ACC PARALLEL LOOP COLLAPSE(2) GANG &
    !$ACC   PRESENT( ngqiq, natmtot, nmtmax, nstspin, ias2ic, &
    !$ACC            batchidx, b2, wftmp1mt ) &
    !$ACC   COPYIN( iblock ) &
    !$ACC   PRIVATE( ibatch, ist, imt, &
    !$ACC            li1w, li1b, lki, list1, liasw, lig, libatch )
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)

          !$ACC LOOP COLLAPSE(2) VECTOR
          DO ki = 1, nstspin
             DO i1 = 1, nmtmax
                imt = i1

#elif defined(_OPENMP) && defined(_PACK_gntuju_)
    ! Zero out wftmp1mt and copy b2 to wftmp1mt (with unpacking)

    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE( tid, ibatch, ist, ic, i1, &
    !$OMP            li1w, li1b, lki, list1, liasw, lig, libatch )
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
          DO ki = 1, nstspin
             DO imt = 1, lmmaxapw*nufrmax
                i1 = irowmap_res(imt,ic,ig)
                IF( i1 == 0 ) THEN
                   wftmp1mt(imt,ki,ias,ig) = zzero
                ELSE

#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
    ! Copy b2 to wftmp1mt

    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE( tid, ibatch, ist, imt, &
    !$OMP            li1w, li1b, lki, list1, liasw, lig, libatch )
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
          DO ki = 1, nstspin
             DO i1 = 1, nmtmax
                imt = i1

#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */

#if EBUG > 2
                ! Check array bounds
                ! i1, imt
                li1w = ( imt >= LBOUND(wftmp1mt,1) ) .AND. &
                       ( imt <= UBOUND(wftmp1mt,1) )
                li1b = ( i1 >= LBOUND(b2,1) )      .AND. &
                       ( i1 <= UBOUND(b2,1) )
                IF( .NOT. li1w ) THEN
                   WRITE(*,*) 'fillresult: imt ', imt, &
                              ' writing wftmp1mt out of bounds', &
                              LBOUND(wftmp1mt,1), UBOUND(wftmp1mt,1)
                   STOP
                END IF
                IF( .NOT. li1b ) THEN
                   WRITE(*,*) 'fillresult: i1 ', i1, &
                              ' reading b2 out of bounds', &
                              LBOUND(b2,1), UBOUND(b2,1)
                   STOP
                END IF

                ! ki, ist1
                list1 = ( ki >= LBOUND(wftmp1mt,2) ) .AND. &
                        ( ki <= UBOUND(wftmp1mt,2) )
                lki   = ( ki >= LBOUND(b2,2) )       .AND. &
                        ( ki <= UBOUND(b2,2) )
                IF( .NOT. list1 ) THEN
                   WRITE(*,*) 'fillresult: ki ', ki, &
                              ' writing wftmp1mt out of bounds', &
                              LBOUND(wftmp1mt,2), UBOUND(wftmp1mt,2)
                   STOP
                END IF
                IF( .NOT. lki ) THEN
                   WRITE(*,*) 'fillresult: ki ', ki, &
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

                wftmp1mt(imt,ki,ias,ig) = b2(i1,ki,ibatch)

#if defined(_OPENACC) && defined(_PACK_gntuju_)
                END IF ! i1 == 0
             END DO ! imt
             !$ACC END LOOP
#elif defined(_OPENACC) && !defined(_PACK_gntuju_)
             END DO ! i1
             !$ACC END LOOP
#elif defined(_OPENMP) && defined(_PACK_gntuju_)
                END IF ! i1 == 0
             END DO ! imt
             !$OMP END SIMD
#elif defined(_OPENMP) && !defined(_PACK_gntuju_)
             END DO ! i1
             !$OMP END SIMD
#endif /* _OPENACC || _OPENMP && _PACK_gntuju_ */

          END DO ! ki
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
