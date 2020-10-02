MODULE mod_genmegqblh_gpu

  USE ISO_C_BINDING ! for C_PTR
  USE mod_prec
  USE mod_gpu
  USE modmain, ONLY: zzero, zone, natmtot, nstsv, nkptnr, &
                     ngrid, ngrtot, ngkmax, ivg, &
                     lmmaxapw, nufrmax, nspinor, spinpol, igfft, ivgig, cfunir
  USE mod_mpi_grid
  USE mod_timer
  USE mod_papi
  USE mod_prec
  USE mod_addons, ONLY: debug_level, dim_k, &
                        pt_megqblh, pt_megqblh_mt, pt_megqblh_it
  USE mod_addons_q, ONLY: ngq, igqig, sfacgq
  USE mod_nrkp, ONLY: spinor_ud
  USE mod_expigqr, ONLY: expigqr22, gntuju, megqblh, bmegqblh, nmegqblh, &
                         idxkq, nbandblhloc, ltranblhloc, ntranblhloc, &
                         idxtranblhloc
#ifdef _USE_3M_
  USE mod_lapack, ONLY: ZGEMM3M
#else
  USE mod_lapack, ONLY: ZGEMM
#endif /* _USE_3M_ */

  IMPLICIT NONE
  
  ! Parameter for genmegqblh_countspin()
  INTEGER, PARAMETER :: spinup =  1
  INTEGER, PARAMETER :: spindn = -1

  ! Table of spin-up/dn states (replaces l1 check)
  ! Dimension is nstsv, but only the first nstspin elements will be used
  INTEGER, DIMENSION(:), ALLOCATABLE :: spinstidx

  ! Number of 2nd-variational states per spin
  INTEGER :: nstspin

  ! Number of muffin-tin elements
  INTEGER :: nmt

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

    ! Number of muffin-tin elements
    nmt = lmmaxapw * nufrmax

#ifdef _OPENACC

    ! Copy/allocate constants to device
    ! Note: natmtot, nspinor, and nstsv are declared in modmain
    !$ACC ENTER DATA CREATE( nstspin ) &
    !$ACC            COPYIN( natmtot, nstsv, nspinor, lmmaxapw, nufrmax, &
    !$ACC                    nblock1, nbatch1, nband1, nmt, ngqiq )

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
    !CALL cudaMalloc( d_nmt, ... )
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
    !$ACC                    nblock1, nbatch1, nband1, nmt, ngqiq )

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
    USE modmain, ONLY: natmtot
    IMPLICIT NONE

    ! Allocate batching index on CPU
    ALLOCATE( batchidx( natmtot, ngqiq, nblock1 ) )

    ! Allocate temporary array to store results on CPU
    ALLOCATE( wftmp1mt( nmt, nstspin, natmtot, ngqiq ) )

    ! Allocate batch arrays for the temporary matrices on CPU
    ! Note: we don't use bgntuju in the OpenACC implementation
    ALLOCATE( b1( nmt, nstspin, nbatch1 ) )
    ALLOCATE( b2( nmt, nstspin, nbatch1 ) )

#ifdef _OPENACC

    ! Allocate array of device pointers on CPU
    ALLOCATE( dptr_gntuju( nbatch1 ) )
    ALLOCATE( dptr_b1( nbatch1 ) )
    ALLOCATE( dptr_b2( nbatch1 ) )

!--DEBUG
!    ALLOCATE( bgntuju( nmt, nmt, nbatch1 ))
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

    ! Allocate batch array for gntuju
    ALLOCATE( bgntuju( nmt, nmt, nbatch1 ))

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
    !$ACC                   b1, b2, & 
    !$ACC                   dptr_gntuju, dptr_b1, dptr_b2 )

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
#else
    IF( ALLOCATED(bgntuju)     ) DEALLOCATE( bgntuju )
#endif /* _OPENACC */

    RETURN
  END SUBROUTINE genmegqblh_freemodvar_mt

!==============================================================================
! Kernel 0: Counts how many 2nd-variational states are spin up/down,
!           and returns a list of such states as spinstidx
!==============================================================================

  SUBROUTINE genmegqblh_countspin( spinproj, ikloc )

    USE modmain, ONLY: nstsv, spinpol
    USE mod_nrkp, ONLY: spinor_ud
    USE mod_expigqr, ONLY: bmegqblh, nbandblhloc, &
                           idxtranblhloc, ltranblhloc, ntranblhloc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: spinproj, ikloc

    ! Internal variables
    INTEGER :: iband, istsp
    LOGICAL :: lup, ldn, lcond, lpaired

    lup = (spinproj == spinup)
    ldn = (spinproj == spindn)

    !$ACC DATA COPYIN( lup, ldn )

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
    !$ACC KERNELS COPYIN( ikloc ) &
    !$ACC   PRESENT( spinstidx, nstspin, lup, ldn, &
    !$ACC            nstsv, spinor_ud, ltranblhloc ) &
    !$ACC   CREATE( istsp )
#else
    !$OMP MASTER
#endif /* _OPENACC */

    ! Initialize values
    nstspin = 0

    ! Begin search algorithm
    ! TODO: Parallelize
#ifdef _OPENACC
    !$ACC LOOP SEQ PRIVATE( iband, lpaired, lcond )
#endif /* _OPENACC */
    DO iband = 1, nstsv

       ! Test whether the band contributes to transitions
       lpaired = ltranblhloc(iband,ikloc)

       IF( lpaired ) THEN

          IF( spinpol ) THEN
             ! Test the condition (Are we counting spin up or spin down states?)
             lcond = ( lup .AND. ( spinor_ud(1,iband,ikloc) == 1 &
                             .AND. spinor_ud(2,iband,ikloc) == 0 ) ) .OR. &
                     ( ldn .AND. ( spinor_ud(1,iband,ikloc) == 0 &
                             .AND. spinor_ud(2,iband,ikloc) == 1 ) )
          ELSE
             ! When spinpol == .FALSE., no need to test the condition
             lcond = .TRUE.
          END IF ! spinpol

          IF( lcond ) THEN

             ! Increment nstspin, then save band index iband
             ! Note: ATOMICs not needed since LOOP SEQ is in effect
             nstspin = nstspin + 1
             spinstidx(nstspin) = iband

          END IF ! lcond
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
    WRITE(*,*) 'countspin: ', spinproj, ' ikloc=', ikloc, ' nstspin=', nstspin
    WRITE(*,*) spinstidx
#endif /* DEBUG */

    IF( nstspin > nband1 ) THEN
       WRITE(*,*) 'Warning[countspin]: nstspin ', nstspin, ' > nband1 ', nband1
    END IF

    ! lup, ldn
    !$ACC END DATA

    RETURN
  END SUBROUTINE genmegqblh_countspin

!==============================================================================
! Kernel 1: Fill in bgntuju (or dptr_gntuju) and b1 arrays, and zero b2 array
!==============================================================================

  SUBROUTINE genmegqblh_fillbatch( wfsvmt1, ikloc, ispn )
    USE modmain, ONLY: zzero, natmtot, nspinor, nstsv, lmmaxapw, nufrmax
    USE mod_expigqr, ONLY: gntuju, bmegqblh, idxtranblhloc
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
    COMPLEX(KIND=dz), DIMENSION(lmmaxapw,nufrmax,natmtot,nspinor,nstsv), INTENT(IN) :: wfsvmt1

    ! Internal variables
    !INTEGER, PARAMETER :: nb = 64          ! Block size for ZGEMM batching
    INTEGER :: iblock                       ! Block index
    INTEGER :: ibatch                       ! Batch index
    !INTEGER :: iblock                      ! Block index
    INTEGER :: k1, k2, ki, nsize           ! Dummy variables for batching
    INTEGER :: iband, i, ist1, ic, ig, ias, i1, i2, imt  ! Data access and/or loop indices
    INTEGER :: tid                         ! Thread ID

!--DEBUG
    LOGICAL :: li1, li2, limt, lki, list1, liasw, liass, lig, lispn, libatch
    INTEGER :: tmp, intmax
    INTEGER :: ist, iold
    INTEGER :: ltmp, lcond, lup, ldn, lrange, lpaired
    LOGICAL :: lcollapse4
!--DEBUG

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
#if defined(_OPENACC)
    !$ACC PARALLEL LOOP COLLAPSE(2) WAIT &
    !$ACC   COPYIN( iblock ) &
    !$ACC   PRESENT( natmtot, ngqiq, batchidx ) &
    !$ACC   PRIVATE( ig, ias, ibatch )
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
#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

    ! Check for integer overflow
    intmax = HUGE(intmax)
    !$ACC DATA CREATE( lcollapse4 ) COPYIN( intmax )
    !$ACC KERNELS PRESENT( natmtot, ngqiq, nstspin, nmt, intmax )
    lcollapse4 = ( DBLE(ngqiq*natmtot) * DBLE(nstspin*nmt) <= intmax )
    !$ACC END KERNELS
    !$ACC UPDATE HOST( lcollapse4 )

#if EBUG >= 2
     WRITE(*,*) 'fillbatch: ispn1=', ispn, ' nstspin=', nstspin
#endif /* DEBUG */

    ! Fill in b1 batch array
#if defined(_OPENACC)
    !$ACC PARALLEL LOOP COLLAPSE(2) GANG &
    !$ACC   COPYIN( iblock, ikloc, ispn ) &
    !$ACC   PRIVATE( ig, ias, ki, i1, ibatch, iband, i, ist1, &
    !$ACC            li1, li2, limt, lki, list1, liasw, liass, lig, lispn, libatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nstspin, nmt, lmmaxapw, nufrmax, &
    !$ACC            batchidx, spinstidx, idxtranblhloc, bmegqblh, &
    !$ACC            wfsvmt1, sfacgq, b1 )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(5) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ibatch, iband, i, ist1, &
    !$OMP            li1, li2, limt, lki, list1, liasw, liass, lig, lispn, libatch )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot
#if defined(_OPENACC)
          !$ACC LOOP COLLAPSE(3) VECTOR &
          !$ACC   PRIVATE( ki, i1, i2, imt, ibatch, iband, i, ist1, &
          !$ACC            li1, li2, limt, lki, list1, liasw, liass, lig, lispn, libatch )
#endif /* _OPENACC */
          DO ki = 1, nstspin
             DO i2 = 1, nufrmax
                DO i1 = 1, lmmaxapw

                   imt = (i2-1)*lmmaxapw + i1

                   ! Note: putting this here because both OpenMP and OpenACC don't like
                   !       breaking up consecutive DO statements, even with comments
                   !DO ki = 1, nsize ! Blocked version

                   ibatch = batchidx(ias,ig,iblock)
                
                   !iband = k1 + ki - 1     ! Blocked version
                   iband = spinstidx( ki ) ! Unblocked version

                   i = idxtranblhloc( iband, ikloc )
                   ist1 = bmegqblh(1,i,ikloc)

#if EBUG > 2
                   ! Check array bounds
                   ! i1, i2, imt
                   li1  = ( i1 >= LBOUND(wfsvmt1,1) ) .AND. ( i1 <= UBOUND(wfsvmt1,1) )
                   li2  = ( i2 >= LBOUND(wfsvmt1,2) ) .AND. ( i2 <= UBOUND(wfsvmt1,2) )
                   limt = ( imt >= LBOUND(b1,1) )     .AND. ( imt <= UBOUND(b1,1) )
                   IF( .NOT. li1 ) THEN
                      WRITE(*,*) 'fillbatch: i1 ', i1, ' reading wfsvmt1 out of bounds', LBOUND(wfsvmt1,1), UBOUND(wfsvmt1,1)
                   END IF
                   IF( .NOT. li2 ) THEN
                      WRITE(*,*) 'fillbatch: i2 ', i2, ' reading wfsvmt1 out of bounds', LBOUND(wfsvmt1,2), UBOUND(wfsvmt1,2)
                   END IF
                   IF( .NOT. limt ) THEN
                      WRITE(*,*) 'fillbatch: imt ', imt, ' writing b1 out of bounds', LBOUND(b1,1), UBOUND(b1,1)
                   END IF
                   ! ki, ist1
                   list1 = ( ist1 >= LBOUND(wfsvmt1,5) ) .AND. ( ist1 <= UBOUND(wfsvmt1,5) )
                   lki   = ( ki >= LBOUND(b1,2) )        .AND. ( ki <= UBOUND(b1,2) )
                   IF( .NOT. list1 ) THEN
                      WRITE(*,*) 'fillbatch: ist1 ', ist1, ' reading wfsvmt1 out of bounds', LBOUND(wfsvmt1,4), UBOUND(wfsvmt1,4)
                   END IF
                   IF( .NOT. lki ) THEN
                      WRITE(*,*) 'fillbatch: ki ', ki, ' writing b1 out of bounds', LBOUND(b1,2), UBOUND(b1,2)
                   END IF
                   ! ias
                   liasw = ( ias >= LBOUND(wfsvmt1,3) ) .AND. ( ias <= UBOUND(wfsvmt1,3) )
                   liass = ( ias >= LBOUND(sfacgq,2) ) .AND. ( ias <= UBOUND(sfacgq,2) )
                   IF( .NOT. liasw ) THEN
                      WRITE(*,*) 'fillbatch: ias ', ias, ' reading wfsvmt1 out of bounds', LBOUND(wfsvmt1,2), UBOUND(wfsvmt1,2)
                   END IF
                   IF( .NOT. liass ) THEN
                      WRITE(*,*) 'fillbatch: ias ', ias, ' reading sfacgq out of bounds', LBOUND(sfacgq,2), UBOUND(sfacgq,2)
                   END IF
                   ! ig
                   lig = ( ig >= LBOUND(sfacgq,1) ) .AND. ( ig <= UBOUND(sfacgq,1) )
                   IF( .NOT. lig ) THEN
                      WRITE(*,*) 'fillbatch: ig ', ig, ' reading sfacgq out of bounds', LBOUND(sfacgq,1), UBOUND(sfacgq,1)
                   END IF
                   ! ispn
                   lispn = ( ispn >= LBOUND(wfsvmt1,4) ) .AND. ( ispn <= UBOUND(wfsvmt1,4) )
                   IF( .NOT. lispn ) THEN
                      WRITE(*,*) 'fillbatch: ispn ', ispn, ' reading wfsvmt1 out of bounds', LBOUND(wfsvmt1,4), UBOUND(wfsvmt1,4)
                   END IF
                   ! ibatch
                   libatch = ( ibatch >= LBOUND(b1,3) ) .AND. ( ibatch <= UBOUND(b1,3) )
                   IF( .NOT. libatch ) THEN
                      WRITE(*,*) 'fillbatch: ibatch ', ibatch, ' writing b1 out of bounds', LBOUND(b1,3), UBOUND(b1,3)
                   END IF
#endif /* DEBUG */

                   ! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
                   b1( imt, ki, ibatch ) = DCONJG( wfsvmt1(i1,i2,ias,ispn,ist1) * &
                                                  sfacgq(ig,ias) )

                END DO ! i1
             END DO ! i2
          END DO ! ki
#ifdef _OPENACC
    !$ACC END LOOP
#endif /* _OPENACC */
       END DO ! ias
    END DO ! ig
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP
#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

    ! Zero b2 batch array
#ifdef _OPENACC
    !$ACC PARALLEL LOOP COLLAPSE(2) GANG &
    !$ACC   COPYIN( iblock, ikloc, ispn ) &
    !$ACC   PRIVATE( ig, ias, i1, i2, ibatch, &
    !$ACC            limt, li2, libatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nband1, nmt, batchidx, b2 )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(4) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ibatch, &
    !$OMP            limt, li2, libatch )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot
#ifdef _OPENACC
    !$ACC LOOP COLLAPSE(2) VECTOR &
    !$ACC   PRIVATE( i1, i2, ibatch, &
    !$ACC            limt, li2, libatch )
#endif /* _OPENACC */
          DO i2 = 1, nstspin
             DO i1 = 1, nmt

                ibatch = batchidx(ias,ig,iblock)

#if EBUG > 2
                ! Check array bounds
                ! imt
                limt = ( imt >= LBOUND(b2,1) ) .AND. ( imt <= UBOUND(b2,1) )
                IF( .NOT. limt ) THEN
                   WRITE(*,*) 'fillbatch: imt ', imt, ' writing b2 out of bounds', LBOUND(b2,1), UBOUND(b2,1)
                END IF
                ! i2
                li2 = ( i2 >= LBOUND(b2,2) ) .AND. ( i2 <= UBOUND(b2,2) )
                IF( .NOT. li2 ) THEN
                   WRITE(*,*) 'fillbatch: i2 ', i2, ' writing b2 out of bounds', LBOUND(b2,2), UBOUND(b2,2)
                END IF
                ! ibatch
                libatch = ( ibatch >= LBOUND(b1,3) ) .AND. ( ibatch <= UBOUND(b1,3) )
                IF( .NOT. libatch ) THEN
                   WRITE(*,*) 'fillbatch: ibatch ', ibatch, ' writing b2 out of bounds', LBOUND(b2,3), UBOUND(b2,3)
                END IF
#endif /* DEBUG */

                b2(i1,i2,ibatch) = zzero

             END DO ! i1
          END DO ! i2
#ifdef _OPENACC
          !$ACC END LOOP
#endif /* _OPENACC */
       END DO ! ias
    END DO ! ig
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP
#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

    ! lcollapse4, intmax
    !$ACC END DATA

    ! Fill in array of device pointers
#ifdef _OPENACC
    !$ACC HOST_DATA USE_DEVICE( gntuju, b1, b2 )

    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC   COPYIN( iblock ) &
    !$ACC   PRIVATE( ig, ias, ic, ibatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, &
    !$ACC            batchidx, ias2ic, &
    !$ACC            gntuju, b1, b2, dptr_gntuju, dptr_b1, dptr_b2 )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP   PRIVATE( ic, ibatch )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)
          ic = ias2ic(ias)

#ifdef _OPENACC
          dptr_gntuju(ibatch) = C_LOC( gntuju(1,1,ic,ig) )
          dptr_b1(ibatch) = C_LOC( b1(1,1,ibatch) )
          dptr_b2(ibatch) = C_LOC( b2(1,1,ibatch) )
#else
          bgntuju(:,:,ibatch) = gntuju(:,:,ic,ig)
#endif /* _OPENACC */

       END DO ! ias
    END DO ! ig
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP

    ! gntuju, b1, b2
    !$ACC END HOST_DATA

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

#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC | _OPENMP */

!  END DO ! k1

    !$ACC WAIT

#endif /* _CUDA_ */

!--DEBUG
!    WRITE(*,*) 'exiting genmegqblh_fillbatch'
!--DEBUG
    
    RETURN
  END SUBROUTINE genmegqblh_fillbatch

!==============================================================================
! Kernel 2: Perform batched ZGEMM b2(:,:) = bgntuju(:,:) x b1(:,:)
!==============================================================================

  SUBROUTINE genmegqblh_batchzgemm()

    !USE modmain, ONLY: zzero, zone

#ifdef _MAGMA_
    USE mod_magma
#elif defined(_CUBLAS_)
#else
    USE mod_gpu, ONLY: zgemm_strided_batched_omp
#endif /* _MAGMA_ */

    IMPLICIT NONE

    ! Internal variables
    COMPLEX(KIND=dz), PARAMETER :: alpha = (1._dd,0._dd)
    COMPLEX(KIND=dz), PARAMETER :: beta  = (0._dd,0._dd)
    INTEGER :: m, n, k, lda, ldb, ldc, ncolA, ncolB, ncolC, stA, stB, stC

    ! Fill in parameters
    m = nmt
    n = nstspin
    k = nmt
    lda = nmt
    ldb = nmt
    ldc = nmt

#if EBUG > 0
    WRITE(*,*) 'batchzgemm: m =', m, ' n = ', n, 'k = ', k
#endif /* DEBUG */
   
  !-2a-------------------------------------------------------------------------
    IF( usemagma ) THEN
  !----------------------------------------------------------------------------

       !$ACC DATA PRESENT( b1, b2, dptr_gntuju, dptr_b1, dptr_b2 )

       ! Note: PARAMETERs don't need to be COPYIN-ed to device

       ! Perform batched ZGEMM on device using MAGMA (pointer mode)
       CALL zgemm_batched_gpu_acc_magma_ptr( 'N', 'N', m, n, k, &
                                             alpha, dptr_gntuju, lda, &
                                                    dptr_b1,     ldb, &
                                             beta,  dptr_b2,     ldc, &
                                             nbatch1 )
#ifdef _MAGMA_
       ! Synchronize with device
       CALL magma_queue_sync( queue )
#endif /* _MAGMA_ */

       ! b1, b2, dptr_gntuju, dptr_b1, dptr_b2
       !$ACC END DATA
       !$ACC WAIT

  !-2b-------------------------------------------------------------------------
    !ELSE IF( usecublas )
  !----------------------------------------------------------------------------

       ! Perform batched ZGEMM on device using cuBLAS
       !CALL cublasZgemmBatched( blashandle, CUBLAS_OP_N, CUBLAS_OP_N, ... )

  !-2c-------------------------------------------------------------------------
    ELSE ! Fall back to CPU only using OpenMP
  !----------------------------------------------------------------------------

!--DEBUG
       !$ACC UPDATE HOST( b1, bgntuju )
!--DEBUG

       ncolA = k ! No transpose
       ncolB = n ! No transpose
       ncolC = n

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
                                       nbatch1 )

  !----------------------------------------------------------------------------
    END IF ! CPU/GPU method
  !----------------------------------------------------------------------------

    RETURN
  END SUBROUTINE genmegqblh_batchzgemm

!==============================================================================
! Kernel 3: Save results to wftmp1mt and transfer back to CPU (for now)
!==============================================================================

  SUBROUTINE genmegqblh_fillresult( wftmp1mt )
    USE modmain, ONLY: natmtot, nstsv

#ifdef _OPENACC
    USE openacc
#endif /* _OPENACC */

#ifdef _OPENMP
     USE omp_lib
#endif /* _OPENMP */
    
    IMPLICIT NONE

    ! Argument
    COMPLEX(KIND=dz), DIMENSION( nmt, nstspin, &
                                 natmtot, ngqiq ) :: wftmp1mt

    ! Internal variables
    INTEGER :: ki, ist, i1, iblock, ibatch, ias, ig, tid

!--DEBUG
    LOGICAL :: li1w, li1b, li2, lki, list1, liasw, lig, libatch
!--DEBUG

#ifdef _CUDA_

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
    
#ifdef _OPENACC

    ! Stub for multi-GPU support
    ! TODO: generalize for AMD GPUs
    !CALL acc_set_device_num( devnum, acc_device_nvidia )

    ! Fill in wftmp1mt on device
    !$ACC PARALLEL LOOP COLLAPSE(2) GANG &
    !$ACC   PRESENT( b2, ngqiq, natmtot, nmt, nstspin, &
    !$ACC            spinstidx, batchidx, wftmp1mt ) &
    !$ACC   PRIVATE( ibatch, ist, ki, &
    !$ACC            li1w, li1b, lki, list1, liasw, lig, libatch ) &
    !$ACC   COPYIN( iblock )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(4) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ibatch, ist, tid, &
    !$OMP            li1w, li1b, lki, list1, liasw, lig, libatch )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot
#ifdef _OPENACC
          !$ACC LOOP COLLAPSE(2) VECTOR &
          !$ACC   PRIVATE( ibatch, ist, ki, &
          !$ACC            li1w, li1b, lki, list1, liasw, lig, libatch )
#endif /* _OPENACC */
          DO ki = 1, nstspin
             DO i1 = 1, nmt

                ist = spinstidx(ki)
                ibatch = batchidx(ias,ig,iblock)

#if EBUG > 2
#ifdef _OPENACC
                ! Check array bounds
                ! i1
                li1w = ( i1 >= LBOUND(wftmp1mt,1) ) .AND. ( i1 <= UBOUND(wftmp1mt,1) )
                li1b = ( i1 >= LBOUND(b2,1) )       .AND. ( i1 <= UBOUND(b2,1) )
                IF( .NOT. li1w ) THEN
                   WRITE(*,*) 'fillresult: i1 ', i1, ' writing wftmp1mt out of bounds', LBOUND(wftmp1mt,1), UBOUND(wftmp1mt,1)
                END IF
                IF( .NOT. li1b ) THEN
                   WRITE(*,*) 'fillresult: i1 ', i1, ' reading b2 out of bounds', LBOUND(b2,1), UBOUND(b2,1)
                END IF
                ! ki, ist1
                list1 = ( ki >= LBOUND(wftmp1mt,2) ) .AND. ( ki <= UBOUND(wftmp1mt,2) )
                lki   = ( ki >= LBOUND(b2,2) )       .AND. ( ki <= UBOUND(b2,2) )
                IF( .NOT. list1 ) THEN
                   WRITE(*,*) 'fillresult: ki ', ki, ' writing wftmp1mt out of bounds', LBOUND(wftmp1mt,2), UBOUND(wftmp1mt,2)
                END IF
                IF( .NOT. lki ) THEN
                   WRITE(*,*) 'fillresult: ki ', ki, ' reading b2 out of bounds', LBOUND(b2,2), UBOUND(b2,2)
                END IF
                ! ias
                liasw = ( ias >= LBOUND(wftmp1mt,3) ) .AND. ( ias <= UBOUND(wftmp1mt,3) )
                IF( .NOT. liasw ) THEN
                   WRITE(*,*) 'fillresult: ias ', ias, ' writing wftmp1mt out of bounds', LBOUND(wftmp1mt,3), UBOUND(wftmp1mt,3)
                END IF
                ! ig
                lig = ( ig >= LBOUND(wftmp1mt,4) ) .AND. ( ig <= UBOUND(wftmp1mt,4) )
                IF( .NOT. lig ) THEN
                   WRITE(*,*) 'fillresult: ig ', ig, ' writing wftmp1mt out of bounds', LBOUND(wftmp1mt,4), UBOUND(wftmp1mt,4)
                END IF
                ! ibatch
                libatch = ( ibatch >= LBOUND(b2,3) ) .AND. ( ibatch <= UBOUND(b2,3) )
                IF( .NOT. libatch ) THEN
                   WRITE(*,*) 'fillresult: ibatch ', ibatch, ' reading b2 out of bounds', LBOUND(b2,3), UBOUND(b2,3)
                END IF
#else
                tid = omp_get_thread_num()
                WRITE(*,*) 'genmegqblh_fillresult: tid=', tid, &
                     ' ias=', ias, ' ig=', ig, ' ibatch=', ibatch, ' ist=', ist
#endif /* _OPENACC */
#endif /* DEBUG */

                wftmp1mt(i1,ki,ias,ig) = b2(i1,ki,ibatch)

             END DO ! i1
          END DO ! ki
#ifdef _OPENACC
          !$ACC END LOOP
#endif /* _OPENACC  */
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
