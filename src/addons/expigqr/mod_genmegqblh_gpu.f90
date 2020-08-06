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
  
#ifdef _CUDA_
  ! Device pointers
  !TYPE(C_PTR) :: d_bgntuju, d_b1, d_b2, d_gntuju, d_sfacgq, d_wfsvmt1
#endif /* _CUDA_ */
  
CONTAINS

!==============================================================================
! Counts how many 2nd-variational states are spin up/down,
! and returns a list of such states as spinstidx
! For now, the contents of spinstidx should be consecutive
! TODO: Perform on device? (rewrite using OpenACC?)

  SUBROUTINE genmegqblh_countspin( spinproj, ikloc )

    USE modmain, ONLY: nstsv, spinpol
    USE mod_nrkp, ONLY: spinor_ud
    USE mod_expigqr, ONLY: idxtranblhloc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: spinproj, ikloc

    ! Internal variables
    INTEGER :: iband, i
    INTEGER :: k1, k2
    LOGICAL :: cond, lup, ldn

    lup = (spinproj == spinup)
    ldn = (spinproj == spindn)

    spinstidx(:) = 0

    IF( spinpol ) THEN

!       !$OMP ATOMIC WRITE
       nstspin = 0
       DO iband = 1, nband1

          i = idxtranblhloc( iband, ikloc )

          ! Test the condition (Are we counting spin up or spin down states?)
          cond = ( lup .AND. ( spinor_ud(1,i,ikloc) == 1 &
                               .AND. spinor_ud(2,i,ikloc) == 0 ) ) .OR. &
                 ( ldn .AND. ( spinor_ud(1,i,ikloc) == 0 &
                               .AND. spinor_ud(2,i,ikloc) == 1 ) )

          IF( cond ) THEN
!             !$OMP ATOMIC
             nstspin = nstspin + 1
!             !$OMP ATOMIC WRITE
             spinstidx(nstspin) = i
          END IF

       END DO ! iband

       ! Check contiguity of states
       k1 = spinstidx(1)
       k2 = spinstidx(nstspin)
       lcontig = ( (k2-k1+1) == nstspin )

    ELSE

       ! If spinpol is .FALSE. there is only one spin projection
       nstspin = nband1
       DO iband = 1, nstspin
          spinstidx(iband) = iband
       END DO ! iband

       ! Contiguity of states is guaranteed
       lcontig = .TRUE.

    END IF ! spinpol

!!!DEBUG
!    WRITE(*,*) 'genmegqblh_countspin: ', spinproj, ' ikloc=', ikloc, ' nstspin=', nstspin
!!!DEBUG

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
    !COMPLEX(KIND=dz), DIMENSION(nmt,nb) :: myb1      ! Blocked ver.
    COMPLEX(KIND=dz), DIMENSION(nmt,nstspin) :: myb1 ! Unblocked ver.
    INTEGER :: ibatch                      ! Batch index
    !INTEGER :: iblock                      ! Block index
    INTEGER :: k1, k2, ki, nsize           ! Dummy variables for batching
    INTEGER :: iband, i, ist1, ic, ig, ias ! Data access and/or loop indices
    INTEGER :: tid                         ! Thread ID

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

    ! Stub for multi-GPU support
    ! TODO: generalize for AMD GPUs
    !CALL acc_set_device_num( devnum, acc_device_nvidia )

    !$ACC PARALLEL LOOP COLLAPSE(2) WAIT &
    !$ACC   PRESENT( nmt, nbatch, bgntuju, b1, b2, &
    !$ACC            gntuju, sfacgq, ias2ic, &
    !$ACC            bmegqblh, idxtranblhloc, &
    !$ACC            natmtot, ngqiq, batchidx, nstspin, spinstidx) &
    !$ACC   COPYIN( iblock, ikloc, ispn, wfsvmt1 ) &
    !$ACC   PRIVATE( ic, ibatch, i, ist1, iband, ki, myb1 )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ig, ias, ic, ibatch, i, ist1, iband, ki, myb1 )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ic = ias2ic(ias)

          ! ibatch = (iblock-1)*natmtot*ngq(iq) + (ig-1)*natmtot + ias
          ! batchidx(ias,ig,iblock) = ibatch

          ibatch = (ig-1)*natmtot + ias
          batchidx(ias,ig,iblock) = ibatch

          ! Loop for a single batch
          myb1(:,:) = zzero
          !DO ki = 1, nsize   ! Blocked ver.
          DO ki = 1, nstspin  ! Unblocked ver.

             !iband = k1 + ki - 1
             iband = spinstidx( ki )
             i = idxtranblhloc( iband, ikloc )
             ist1 = bmegqblh(1,i,ikloc)

             ! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
             myb1(1:nmt,ki) = DCONJG( wfsvmt1(1:nmt,ias,ispn,ist1) * &
                                      sfacgq(ig,ias) )

          END DO ! ki

#ifndef _OPENACC
          !$OMP CRITICAL
#endif /* _OPENACC */
          bgntuju(:,:,ibatch) = gntuju(:,:,ic,ig)
          b1(:,:,ibatch) = myb1(:,:)
          b2(:,:,ibatch) = zzero
#ifndef _OPENACC
          !$OMP END CRITICAL
#endif /* _OPENACC */

       END DO ! ias
    END DO ! ig
#ifdef _OPENACC
    !$ACC END PARALLEL LOOP
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
    COMPLEX(KIND=dz) :: zzero, zone

    zzero = (0._dd,0._dd)
    zone  = (1._dd,0._dd)
   
  !-2a-------------------------------------------------------------------------
    IF( usemagma ) THEN
  !----------------------------------------------------------------------------

       !$ACC DATA COPYIN( bgntuju, b1 ) COPY( b2 )

!       !$ACC DATA COPYIN( zzero, zone, nmt, nstspin, nbatch )
       !$ACC DATA COPYIN( zzero, zone )
       
       ! Perform batched ZGEMM on device using MAGMA
       CALL zgemm_batched_gpu_acc_magma( 'N', 'N', nmt, nstspin, nmt, &
                                    zone,  bgntuju(1:nmt,1:nmt,1:nbatch), nmt, &
                                           b1(1:nmt,1:nstspin,1:nbatch),  nmt, &
                                    zzero, b2(1:nmt,1:nstspin,1:nbatch),  nmt, &
                                    nbatch )
#ifdef _MAGMA_
       ! Synchronize with device
       CALL magma_queue_sync( queue )
#endif /* _MAGMA_ */

       ! zzero, zone
       !$ACC END DATA
       
       ! bgntuju, b1, b2
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
       CALL zgemm_batched_omp( 'N', 'N', nmt, nstspin, nmt, &
                               zone,  bgntuju(1:nmt,1:nmt,1:nbatch), nmt, &
                                      b1(1:nmt,1:nstspin,1:nbatch),  nmt, &
                               zzero, b2(1:nmt,1:nstspin,1:nbatch),  nmt, &
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
    USE modmain, ONLY: natmtot

#ifdef _OPENACC
    USE openacc
#endif /* _OPENACC */

#ifdef _OPENMP
     USE omp_lib
#endif /* _OPENMP */
    
    IMPLICIT NONE

    ! Argument
    COMPLEX(KIND=dz), DIMENSION( nmt, nband1, &
                                 natmtot, ngqiq ) :: wftmp1mt

    ! Internal variables
    INTEGER :: k1, k2, ki, ist, iblock, ibatch, ias, ig, tid

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
    
    IF( lcontig ) THEN
       k1 = spinstidx(1)
       k2 = spinstidx(nstspin)
    END IF ! lcontig

#ifdef _OPENACC    

    ! Stub for multi-GPU support
    ! TODO: generalize for AMD GPUs
    !CALL acc_set_device_num( devnum, acc_device_nvidia )

    ! Fill in wftmp1mt on device
    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC   PRESENT( b2, ngqiq, natmtot, nmt, nstspin, &
    !$ACC            spinstidx, batchidx, wftmp1mt, lcontig ) &
    !$ACC   PRIVATE( ibatch, ist, ki ) &
    !$ACC   COPYIN( k1, k2, iblock )
#elif defined(_OPENMP)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE( ig, ias, ibatch, ist, ki, tid )
#endif /* _OPENACC || _OPENMP */
    DO ig = 1, ngqiq
       DO ias = 1, natmtot

          ibatch = batchidx(ias,ig,iblock)

          ! If contiguous
          IF( lcontig ) THEN
#if EBUG > 2 && !defined(_OPENACC)
!!!DEBUG
           tid = omp_get_thread_num()
           WRITE(*,*) 'genmegqblh_fillresult: tid=', tid, &
                      ' ias=', ias, ' ig=', ig, ' ibatch=', ibatch, ' k1=', k1, ' k2=', k2
!!!DEBUG
#endif /* DEBUG */
#ifndef _OPENACC
             !$OMP CRITICAL
#endif /* _OPENACC */
             wftmp1mt(1:nmt,k1:k2,ias,ig) = b2(1:nmt,1:nstspin,ibatch)
#ifndef _OPENACC
             !$OMP END CRITICAL
#endif /* _OPENACC */
          ELSE
             DO ist = 1, nstspin
                ki = spinstidx(ist)
#if EBUG > 2 && !defined(_OPENACC)
!!!DEBUG
                tid = omp_get_thread_num()
                WRITE(*,*) 'genmegqblh_fillresult: tid=', tid, &
                     ' ias=', ias, ' ig=', ig, ' ibatch=', ibatch, ' ki=', ki
!!!DEBUG
#endif /* DEBUG */
#ifndef _OPENACC
                !$OMP CRITICAL
#endif /* _OPENACC */
                wftmp1mt(1:nmt,ki,ias,ig) = b2(1:nmt,ist,ibatch)
#ifndef _OPENACC
                !$OMP END CRITICAL
#endif /* _OPENACC */
             END DO ! ist
          END IF ! lcontig
        
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
