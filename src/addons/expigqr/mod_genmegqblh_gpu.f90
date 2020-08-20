MODULE mod_genmegqblh_gpu

  USE ISO_C_BINDING ! for C_PTR
  USE mod_prec
  USE mod_gpu

  IMPLICIT NONE
  
  ! Parameter for genmegqblh_countspin()
  INTEGER, PARAMETER :: spinup =  1
  INTEGER, PARAMETER :: spindn = -1

  ! Table of spin-up/dn states (replaces l1 check)
  ! Dimension is nstsv, but only the first nstspin elements will be used
  INTEGER, DIMENSION(:), ALLOCATABLE :: spinstidx

  ! Number of 2nd-variational states per spin
  INTEGER :: nstspin

  ! Flag for whether the states in spinstidx are contiguous
  LOGICAL :: lcontig

  ! Number of muffin-tin elements
  INTEGER :: nmt

  ! Block size for genmegqblh_batchzgemm()
  !INTEGER, PARAMETER :: nb = 64

  ! Number of blocks
  INTEGER :: nblock
  
  ! Number of batches
  INTEGER :: nbatch

  ! Number of bands associated with the bra state vectors
  INTEGER :: nband1
  
  ! Number of G+q vectors for a particular value of q-vector
  INTEGER :: ngqiq

  ! Translation table for each batch index
  ! Dimensions are natmtot, ngqiq, nblock, respectively
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: batchidx

  ! Matrices for genmegqblh_batchzgemm()
  COMPLEX(KIND=dz), DIMENSION(:,:,:), ALLOCATABLE :: bgntuju, b1, b2

  ! Array of device pointers for genmegqblh_batchzgemm()
  TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: dptr_gntuju, dptr_b1, dptr_b2
  
#ifdef _CUDA_
  ! Device pointers
  !TYPE(C_PTR) :: d_bgntuju, d_b1, d_b2, d_gntuju, d_sfacgq, d_wfsvmt1
#endif /* _CUDA_ */
  
CONTAINS

!==============================================================================
! Counts how many 2nd-variational states are spin up/down,
! and returns a list of such states as spinstidx

  SUBROUTINE genmegqblh_countspin( spinproj, ikloc )

    USE modmain, ONLY: nstsv, spinpol
    USE mod_nrkp, ONLY: spinor_ud
    USE mod_expigqr, ONLY: idxtranblhloc, bmegqblh, &
                           idxlobandblhloc, idxhibandblhloc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: spinproj, ikloc

    ! Internal variables
    INTEGER, DIMENSION(nstsv) :: tmp
    LOGICAL, DIMENSION(nstsv) :: ltmp
    INTEGER :: iband, i, ist, iold, ilo, ihi
    INTEGER :: k1, k2
    LOGICAL :: lcond, lup, ldn, lrange, lpaired

    lup = (spinproj == spinup)
    ldn = (spinproj == spindn)

    ! Zero out arrays
#ifdef _OPENACC
    !$ACC DATA CREATE( tmp, ltmp )
    !$ACC PARALLEL LOOP PRESENT( spinstidx, nstsv )
#else
    !$OMP PARALLEL DO
#endif /* _OPENACC */
    DO iband = 1, nstsv
       spinstidx(iband) = 0
       tmp(iband) = 0
       ltmp(iband) = .FALSE.
    END DO ! nstsv
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP
#else
    !$OMP END PARALLEL DO
#endif /* _OPENACC */

    IF( spinpol ) THEN

       IF( lup ) THEN
          ilo = idxlobandblhloc(1,ikloc)
          ihi = idxhibandblhloc(1,ikloc)
       ELSE IF( ldn ) THEN
          ilo = idxlobandblhloc(2,ikloc)
          ihi = idxhibandblhloc(2,ikloc)
       END IF

       ! Fill in temporary array
#ifdef _OPENACC
       !$ACC PARALLEL LOOP &
       !$ACC   COPYIN( lup, ldn, ikloc, ilo, ihi ) &
       !$ACC   PRESENT( spinor_ud, idxtranblhloc, bmegqblh, nband1 ) &
       !$ACC   PRIVATE( i, ist, lcond, lrange, lpaired )
#else
       !$OMP PARALLEL DO &
       !$OMP   PRIVATE( i, ist, lcond, lrange, lpaired )
#endif /* _OPENACC */
       DO iband = 1, nstsv

          ! Test whether the band is within range
          lrange = ((iband >= ilo) .AND. (iband <= ihi))

          IF( lrange ) THEN

             i = idxtranblhloc(iband,ikloc)
             ist = bmegqblh(1,i,ikloc)

             ! Test whether the band is contributing to transitions
             lpaired = ( ist /= 0 )
             IF( lpaired ) THEN

                ! Test the condition (Are we counting spin up or spin down states?)
                lcond = ( lup .AND. ( spinor_ud(1,ist,ikloc) == 1 &
                                    .AND. spinor_ud(2,ist,ikloc) == 0 ) ) .OR. &
                        ( ldn .AND. ( spinor_ud(1,ist,ikloc) == 0 &
                                    .AND. spinor_ud(2,ist,ikloc) == 1 ) )

                IF( lcond ) THEN
                   ltmp(iband) = .TRUE.
                   tmp(iband) = ist
                END IF ! lcond

             END IF ! lpaired

          END IF ! lrange

       END DO ! nstsv
#ifdef _OPENACC
       !$ACC END PARALLEL LOOP
#else
       !$OMP END PARALLEL DO
#endif /* _OPENACC */

#ifdef _OPENACC
       !$ACC KERNELS COPYIN( nstsv ) CREATE( i ) &
       !$ACC   PRESENT( spinstidx, nstspin, lcontig, ltmp, tmp )
#else
       !$OMP MASTER
#endif

       ! Count number of elements
       nstspin = COUNT( ltmp )

       ! Filter out zeroes
#ifdef _OPENACC
       i = 1
       !$ACC LOOP SEQ PRIVATE( iband, iold )
       DO iband = 1, nstsv
          IF( ltmp(iband) ) THEN
             iold = i
             spinstidx(iold) = tmp(iband)
             i = iold + 1
          END IF
       END DO ! nstsv
       !$ACC END LOOP
#else
       spinstidx(1:nstspin) = PACK( tmp, ltmp )
#endif /* _OPENACC */

       ! Check contiguity of states
       k1 = spinstidx(1)
       k2 = spinstidx(nstspin)
       lcontig = ( (k2-k1+1) == nstspin )

#ifdef _OPENACC
       !$ACC END KERNELS
#else
       !$OMP END MASTER
#endif /* _OPENACC */

    ELSE

#ifdef _OPENACC
       !$ACC KERNELS PRESENT( nband1, spinstidx, nstspin, lcontig )
#else
       !$OMP MASTER
#endif /* _OPENACC */

       ! If spinpol is .FALSE. there is only one spin projection
       nstspin = nband1

       ! Contiguity of states is guaranteed
       lcontig = .TRUE.

#ifdef _OPENACC
       !$ACC END KERNELS
#endif /* _OPENACC */

#ifdef _OPENACC
       !$ACC LOOP SEQ PRIVATE(iband)
#endif /* _OPENACC */
       DO iband = 1, nband1
          i = idxtranblhloc(iband,ikloc)
          ist = bmegqblh(1,i,ikloc)
          spinstidx(iband) = ist
       END DO ! iband
#ifdef _OPENACC
       !$ACC END LOOP
#else
       !$OMP END MASTER
#endif /* _OPENACC */

    END IF ! spinpol

    ! tmp, ltmp
    !$ACC END DATA

    !$ACC UPDATE HOST( spinstidx, nstspin, lcontig )
    !$ACC WAIT

#if EBUG > 1
    WRITE(*,*) 'genmegqblh_countspin: ', spinproj, ' ikloc=', ikloc, ' nstspin=', nstspin
    WRITE(*,*) spinstidx
#endif /* DEBUG */

    IF( nstspin > nband1 ) THEN
       WRITE(*,*) 'Warning[countspin]: nstspin ', nstspin, ' > nband1 ', nband1
    END IF

    RETURN
  END SUBROUTINE genmegqblh_countspin

!==============================================================================
! Kernel 1: Fill in bgntuju and b1, and zero b2
!==============================================================================

  SUBROUTINE genmegqblh_fillbatch( wfsvmt1, ikloc, ispn )
    USE modmain, ONLY: zzero, natmtot, nspinor, nstsv
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
    COMPLEX(KIND=dz), DIMENSION(nmt,natmtot,nspinor,nstsv), INTENT(IN) :: wfsvmt1

    ! Internal variables
    !INTEGER, PARAMETER :: nb = 64          ! Block size for ZGEMM batching
    INTEGER :: iblock                       ! Block index
    INTEGER :: ibatch                       ! Batch index
    !INTEGER :: iblock                      ! Block index
    INTEGER :: k1, k2, ki, nsize           ! Dummy variables for batching
    INTEGER :: iband, i, ist1, ic, ig, ias, i1, i2  ! Data access and/or loop indices
    INTEGER :: tid                         ! Thread ID

!--DEBUG
    LOGICAL :: li1w, li1b, li2, lki, list1, liasw, liass, lig, lispn, libatch
!--DEBUG

#ifdef _CUDA_

    ! Allocate device pointers
    !CALL cudaMalloc( d_wfsvmt1, ... )
    !CALL cudaMalloc( d_sfacgq, ... )
    !CALL cudaMalloc( d_gntuju, ... )
    !CALL cudaMalloc( d_bgntuju, ... )
    !CALL cudaMalloc( d_b1, ... )
    !CALL cudaMalloc( d_b2, ... )

    ! Transfer data from host to device
    !CALL cudaMemcpy( d_wfsvmt1, wfsvmt1, cudaMemcpyHostToDevice )
    !CALL cudaMemcpy( d_sfacgq,  sfacgq,  cudaMemcpyHostToDevice )
    !CALL cudaMemcpy( d_gntuju,  gntuju,  cudaMemcpyHostToDevice )

    ! Call CUDA C++ kernel
    !CALL genmegqblh_fillbatch_cu_( ... )

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

    ! Transfer wfsvmt1 to device
    !$ACC DATA COPYIN( wfsvmt1 )

    ! Stub for multi-GPU support
    ! TODO: generalize for AMD GPUs
    !CALL acc_set_device_num( devnum, acc_device_nvidia )
#endif /* _OPENACC */

    ! Fill in batchidx, the translation table for ibatch <-> {ig,ias,iblock}
#ifdef _OPENACC
    !$ACC PARALLEL LOOP COLLAPSE(2) WAIT &
    !$ACC   COPYIN( iblock ) &
    !$ACC   PRESENT( natmtot, ngqiq, batchidx ) &
    !$ACC   PRIVATE( ig, ias, ibatch )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ig, ias, ibatch )
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
    !lcollapse4 = ( DBLE(ngqiq*natmtot) * DBLE(nstspin*nmt) <= intmax )
!    !$ACC DATA COPYIN( lcollapse4 )

    ! Fill in b1 batch array
#ifdef _OPENACC
    !$ACC PARALLEL LOOP COLLAPSE(2) GANG &
    !$ACC   COPYIN( iblock, ikloc, ispn ) &
    !$ACC   PRIVATE( ig, ias, ki, i1, ibatch, iband, i, ist1, &
    !$ACC            li1w, li1b, lki, list1, liasw, liass, lig, lispn, libatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nstspin, nmt, &
    !$ACC            batchidx, spinstidx, idxtranblhloc, bmegqblh, &
    !$ACC            wfsvmt1, sfacgq, b1 )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(4) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ig, ias, ki, i1, ibatch, iband, i, ist1, &
    !$OMP            li1w, li1b, lki, list1, liasw, liass, lig, lispn, libatch )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot
#ifdef _OPENACC
    !$ACC LOOP COLLAPSE(2) VECTOR &
    !$ACC   PRIVATE( ki, i1, ibatch, iband, i, ist1, &
    !$ACC            li1w, li1b, lki, list1, liasw, liass, lig, lispn, libatch )
#endif /* _OPENACC */
          DO ki = 1, nstspin
             DO i1 = 1, nmt

          ! Note: putting this here because both OpenMP and OpenACC don't like
          !       breaking up consecutive DO statements, even with comments
          !DO ki = 1, nsize ! Blocked version

                ibatch = batchidx(ias,ig,iblock)
                
                !iband = k1 + ki - 1     ! Blocked version
                ist1 = spinstidx( ki ) ! Unblocked version

                !i = idxtranblhloc( iband, ikloc )
                !ist1 = bmegqblh(1,i,ikloc)

#if EBUG > 2
                ! Check array bounds
                ! i1
                li1w = ( i1 >= LBOUND(wfsvmt1,1) ) .AND. ( i1 <= UBOUND(wfsvmt1,1) )
                li1b = ( i1 >= LBOUND(b1,1) )      .AND. ( i1 <= UBOUND(b1,1) )
                IF( .NOT. li1w ) THEN
                   WRITE(*,*) 'fillbatch: i1 ', i1, ' reading wfsvmt1 out of bounds', LBOUND(wfsvmt1,1), UBOUND(wfsvmt1,1)
                END IF
                IF( .NOT. li1b ) THEN
                   WRITE(*,*) 'fillbatch: i1 ', i1, ' writing b1 out of bounds', LBOUND(b1,1), UBOUND(b1,1)
                END IF
                ! ki, ist1
                list1 = ( ist1 >= LBOUND(wfsvmt1,4) ) .AND. ( ist1 <= UBOUND(wfsvmt1,4) )
                lki   = ( ki >= LBOUND(b1,2) )        .AND. ( ki <= UBOUND(b1,2) )
                IF( .NOT. list1 ) THEN
                   WRITE(*,*) 'fillbatch: ist1 ', ist1, ' reading wfsvmt1 out of bounds', LBOUND(wfsvmt1,4), UBOUND(wfsvmt1,4)
                END IF
                IF( .NOT. lki ) THEN
                   WRITE(*,*) 'fillbatch: ki ', ki, ' writing b1 out of bounds', LBOUND(b1,2), UBOUND(b1,2)
                END IF
                ! ias
                liasw = ( ias >= LBOUND(wfsvmt1,2) ) .AND. ( ias <= UBOUND(wfsvmt1,2) )
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
                lispn = ( ispn >= LBOUND(wfsvmt1,3) ) .AND. ( ispn <= UBOUND(wfsvmt1,3) )
                IF( .NOT. lispn ) THEN
                   WRITE(*,*) 'fillbatch: ispn ', ispn, ' reading wfsvmt1 out of bounds', LBOUND(wfsvmt1,3), UBOUND(wfsvmt1,3)
                END IF
                ! ibatch
                libatch = ( ibatch >= LBOUND(b1,3) ) .AND. ( ibatch <= UBOUND(b1,3) )
                IF( .NOT. libatch ) THEN
                   WRITE(*,*) 'fillbatch: ibatch ', ibatch, ' writing b1 out of bounds', LBOUND(b1,3), UBOUND(b1,3)
                END IF
#endif /* DEBUG */

                ! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
                b1( i1, ki, ibatch ) = DCONJG( wfsvmt1(i1,ias,ispn,iband) * &
                                               sfacgq(ig,ias) )
             END DO ! i1
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
    !$ACC            li1b, li2, libatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, nband1, nmt, batchidx, b2 )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(4) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ig, ias, i1, i2, ibatch, &
    !$OMP            li1b, li2, libatch )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot
#ifdef _OPENACC
    !$ACC LOOP COLLAPSE(2) VECTOR &
    !$ACC   PRIVATE( i1, i2, ibatch, &
    !$ACC            li1b, li2, libatch )
#endif /* _OPENACC */
          DO i2 = 1, nstspin
             DO i1 = 1, nmt

                ibatch = batchidx(ias,ig,iblock)

#if EBUG > 2
                ! Check array bounds
                ! i1
                li1b = ( i1 >= LBOUND(b2,1) ) .AND. ( i1 <= UBOUND(b2,1) )
                IF( .NOT. li1b ) THEN
                   WRITE(*,*) 'fillbatch: i1 ', i1, ' writing b2 out of bounds', LBOUND(b2,1), UBOUND(b2,1)
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

    ! lcollapse4
!    !$ACC END DATA

    ! Fill in array of device pointers
#ifdef _OPENACC
    !$ACC HOST_DATA USE_DEVICE( gntuju, b1, b2 )

    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC   COPY( iblock ) &
    !$ACC   PRIVATE( ig, ias, ic, ibatch ) &
    !$ACC   PRESENT( natmtot, ngqiq, &
    !$ACC            batchidx, ias2ic, &
    !$ACC            gntuju, b1, b2, dptr_gntuju, dptr_b1, dptr_b2 )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP   PRIVATE( ig, ias, ic, ibatch )
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

    ! wfsvmt1
    !$ACC END DATA
#elif defined(_OPENMP)
    !$OMP END PARALLEL DO
#endif /* _OPENACC | _OPENMP */

!  END DO ! k1

!    !$ACC WAIT

#endif /* _CUDA_ */

!--DEBUG
!    WRITE(*,*) 'exiting genmegqblh_fillbatch'
!--DEBUG
    
    RETURN
  END SUBROUTINE genmegqblh_fillbatch

!==============================================================================
! Kernel 2: Perform batched ZGEMM b2(:,:) = b1(:,:) x bgntuju(:,:)
!==============================================================================

  SUBROUTINE genmegqblh_batchzgemm()

    !USE modmain, ONLY: zzero, zone

#ifdef _MAGMA_
    USE mod_magma
#endif /* _MAGMA_ */

    IMPLICIT NONE

    ! Internal variables
    COMPLEX(KIND=dz), PARAMETER :: alpha = (1._dd,0._dd)
    COMPLEX(KIND=dz), PARAMETER :: beta  = (0._dd,0._dd)
    INTEGER :: m, n, k, lda, ldb, ldc

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

       !$ACC DATA COPYIN( b1 ) COPY( b2 ) &
       !$ACC   PRESENT( dptr_gntuju, dptr_b1, dptr_b2 )

       ! Note: PARAMETERs don't need to be COPYIN-ed to device

       ! Perform batched ZGEMM on device using MAGMA (pointer mode)
       CALL zgemm_batched_gpu_acc_magma_ptr( 'N', 'N', m, n, k, &
                                             alpha, dptr_gntuju, lda, &
                                                    dptr_b1,     ldb, &
                                             beta,  dptr_b2,     ldc, &
                                             nbatch )
#ifdef _MAGMA_
       ! Synchronize with device
       CALL magma_queue_sync( queue )
#endif /* _MAGMA_ */

       ! b1, b2, dptr_gntuju, dptr_b1, dptr_b2
       !$ACC END DATA

  !-2b-------------------------------------------------------------------------
    !ELSE IF( usecublas )
  !----------------------------------------------------------------------------

       ! Perform batched ZGEMM on device using cuBLAS
       !CALL cublasZgemmBatched( blashandle, CUBLAS_OP_N, CUBLAS_OP_N, ... )

  !-2c-------------------------------------------------------------------------
    ELSE ! Fall back to CPU only using OpenMP
  !----------------------------------------------------------------------------

       ! Perform batched ZGEMM on CPU using OpenMP parallel do
       ! b2(1:nmt,1:nstsvup) = bgntuju(1:nmt,1:nmt) x b1(1:nmt,1:nstsv)
       CALL zgemm_batched_omp( 'N', 'N', m, n, k, &
                               alpha, bgntuju, lda, &
                                      b1,      ldb, &
                               beta,  b2,      ldc, &
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
    USE modmain, ONLY: natmtot, nstsv

#ifdef _OPENACC
    USE openacc
#endif /* _OPENACC */

#ifdef _OPENMP
     USE omp_lib
#endif /* _OPENMP */
    
    IMPLICIT NONE

    ! Argument
    COMPLEX(KIND=dz), DIMENSION( nmt, nstsv, &
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
    !$ACC            spinstidx, batchidx, wftmp1mt, lcontig ) &
    !$ACC   PRIVATE( ibatch, ist, ki, &
    !$ACC            li1w, li1b, lki, list1, liasw, lig, libatch ) &
    !$ACC   COPYIN( iblock )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(4) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ig, ias, ibatch, ist, ki, tid, &
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
!     !$ACC WAIT
#elif defined(_OPENMP)
     !$OMP END PARALLEL DO
#endif /* _OPENACC || _OPENMP */

#endif /* _CUDA_ */

    RETURN
  END SUBROUTINE genmegqblh_fillresult

!==============================================================================

END MODULE mod_genmegqblh_gpu
