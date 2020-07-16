MODULE mod_genmegqblh_gpu

  USE ISO_C_BINDING ! for C_PTR
  USE mod_gpu

  ! Parameter for genmegqblh_countspin()
  INTEGER, PARAMETER :: spinup =  1
  INTEGER, PARAMETER :: spindn = -1

#ifdef _CUDA_
  ! Device pointers
  !TYPE(C_PTR) :: d_bgntuju, d_b1, d_b2, d_gntuju, d_sfacgq, d_wfsvmt1
#endif /* _CUDA_ */

CONTAINS

!==============================================================================
! Counts how many 2nd-variational states are spin up,
! and returns a list of such states as stateidx
! For now, the contents of stateidx should be consecutive
! TODO: Perform on device? (rewrite using OpenACC?)

  FUNCTION genmegqblh_countspin( spinproj, ikloc, nstspin ) RESULT (spinstidx)

    USE modmain     ! for nstsv
    USE mod_nrkp    ! for spinor_ud
    USE mod_expigqr ! for idxhibandblhloc(:), idxtranblhloc(:,:)

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: spinproj, ikloc
    INTEGER, INTENT(OUT) :: nstspin
    INTEGER, DIMENSION(nstsv) :: spinstidx

    ! Internal variables
    INTEGER :: iband, i, ist1
    LOGICAL :: cond

    IF( spinpol ) THEN

       nstspin = 0
       spinstidx(:) = 0
       DO iband = 1, idxhibandblhloc(ikloc)

          i = idxtranblhloc( iband, ikloc )
          ist1 = bmegqblh(1,i,ikloc)

          SELECT CASE( spinproj )
          CASE( spinup )
             cond = spinor_ud(1,i,ikloc) == 1 .AND. spinor_ud(2,i,ikloc) == 0
          CASE( spindn )
             cond = spinor_ud(1,i,ikloc) == 0 .AND. spinor_ud(2,i,ikloc) == 1
          END SELECT ! spinproj

          IF( cond ) THEN
             nstspin = nstspin + 1
             spinstidx(nstspin) = i
          END IF

       END DO ! iband

    ELSE

       ! If spinpol is .FALSE. there is only one spin projection
       nstspin = nstsv
       DO iband = 1, idxhibandblhloc(ikloc)
          spinstidx(iband) = iband
       END DO ! iband

    END IF ! spinpol

!!!DEBUG
    WRITE(*,*) 'genmegqblh_countspin: ', spinproj, ' ikloc=', ikloc, ' nstspin=', nstspin
!!!DEBUG

  END FUNCTION genmegqblh_countspin

!==============================================================================
! Kernel 1: Fill in bgntuju and b1, and zero b2
!==============================================================================

  SUBROUTINE genmegqblh_fillbatch( bgntuju, b1, b2, batchidx, &
                                   wfsvmt1, nmt, nstspin, &
                                   iq, ikloc, ispn1, spinstidx )
    IMPLICIT NONE

    ! (dummy) Arguments
    COMPLEX(KIND=dz), DIMENSION(:,:,:) :: bgntuju, b1, b2
    COMPLEX(KIND=dz), DIMENSION(:,:,:,:) :: wfsvmt1
    INTEGER, DIMENSION(:,:,:) :: batchidx
    INTEGER, DIMENSION(:) :: spinstidx
    INTEGER :: nmt, nstspin, iq, ikloc, ispn1

  !-1a-------------------------------------------------------------------------
    IF( useacc .AND. usemagma ) THEN
  !----------------------------------------------------------------------------

       !$ACC ENTER DATA CREATE( bgntuju, b1, b2, batchidx )
       CALL genmegqblh_fillbatch_acc( bgntuju, b1, b2, batchidx, &
                                      wfsvmt1, nmt, nstspin, &
                                      iq, ikloc, ispn1, spinstidx )

  !-1b-------------------------------------------------------------------------
    !ELSE IF( usecuda .AND. usecublas )
  !----------------------------------------------------------------------------

       !CALL cudaMemcpy( d_wfsvmt1, wfsvmt1, cudaMemcpyHostToDevice )
       !CALL cudaMemcpy( d_sfacgq,  sfacgq,  cudaMemcpyHostToDevice )
       !CALL cudaMemcpy( d_gntuju,  gntuju,  cudaMemcpyHostToDevice )

       !CALL cudaMalloc( d_bgntuju, ... )
       !CALL cudaMalloc( d_b1, ... )
       !CALL cudaMalloc( d_b2, ... )

       !CALL genmegqblh_fillbatch_cuda( d_bgntuju, d_b1, d_b2, batchidx, &
       !                                d_gntuju, d_sfacgq, d_wfsvmt1, &
       !                                nmt, nstspin, &
       !                                iq, ikloc, ispn, spinstidx )

  !-1c-------------------------------------------------------------------------
    ELSE ! Fall back to CPU only using OpenMP
  !----------------------------------------------------------------------------

     ! Fill in bgntuju and b1 on CPU
       CALL genmegqblh_fillbatch_omp( bgntuju, b1, b2, batchidx, &
                                      wfsvmt1, nmt, nstspin, &
                                      iq, ikloc, ispn1, spinstidx )

  !----------------------------------------------------------------------------
    END IF ! CPU/GPU method
  !----------------------------------------------------------------------------

    RETURN
  END SUBROUTINE genmegqblh_fillbatch

!==============================================================================
! Kernel 1a - OpenACC version
! Fill in bgntuju and b1 (and zero b2) using OpenACC parallel loop
! Also outputs "translation table" for (ias,ig,iblock) to ibatch as batchidx

  SUBROUTINE genmegqblh_fillbatch_acc( bgntuju, b1, b2, batchidx, &
                                       wfsvmt1, nmt, nstsvup, &
                                       iq, ikloc, ispn, spinupidx )

    USE modmain      ! for dz, natmtot
    USE mod_expigqr  ! for gntuju, bmegqblh
    USE mod_addons_q ! for sfacgq, ngq(iq)
    USE mod_nrkp     ! for spinor_ud
    USE mod_gpu      ! for devnum

#ifdef _OPENACC
    USE openacc
#endif /* _OPENACC */

    IMPLICIT NONE

    ! Arguments
    COMPLEX(KIND=dz), DIMENSION(:,:,:), INTENT(OUT) :: bgntuju, b1, b2
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: batchidx
    COMPLEX(KIND=dz), DIMENSION(:,:,:,:), INTENT(IN) :: wfsvmt1
    INTEGER, INTENT(IN) :: nmt, nstsvup, iq, ikloc, ispn
    INTEGER, DIMENSION(:), INTENT(IN) :: spinupidx
    ! Note: other ingredients (gntuju,sfacgq) are passed through modules

#ifdef _OPENACC

    ! Internal variables
    !INTEGER, PARAMETER :: nb = 64          ! Block size for ZGEMM batching
    INTEGER :: ibatch, nbatch, nblock      ! Batch index and number of batches
    INTEGER :: k1, k2, ki, nsize           ! Dummy variables for batching
    INTEGER :: iband, i, ist1, ic, ig, ias ! Data access and/or loop indices

    ! Number of batches
    !nblock = CEILING( REAL(idxhiband)/REAL(nb) )
    !nbatch = ngq(iq) * natmtot * nblock
    nbatch = ngq(iq) * natmtot

    ! Batching by block size nb for idxhiband
    ! TODO: Re-enable if needed for larger systems

    !  idxhiband = idxhibandblhloc(ikloc)
    !  ibatch = 0
    !  DO k1 = 1, idxhiband, nb
    !     k2 = MIN( idxhiband, k1+nb-1 )
    !     nsize = k2 - k1 + 1
    !     iblock = ???

    ! For consistency
    ! TODO: generalize for AMD GPUs
    CALL acc_set_device_num( devnum, acc_device_nvidia )

    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC   PRESENT( bgntuju, b1, b2, gntuju, sfacgq, wfsvmt1, &
    !$ACC            bmegqblh(:,:,ikloc), idxtranblhloc(:,ikloc), spinor_ud, &
    !$ACC            ngq(iq), ias2ic, batchidx ) &
    !$ACC   CREATE( ic, ibatch, i, ist1, iband, ki ) &
    !$ACC   COPYIN( natmtot, nstsvup, nmt )
    DO ig = 1, ngq(iq)
       DO ias = 1, natmtot

          ic = ias2ic(ias)

          ! ibatch = ???
          ! batchidx(ias,ig,iblock) = ibatch

          ibatch = (ig-1)*natmtot + ias
          batchidx(ias,ig,1) = ibatch

          bgntuju(:,:,ibatch) = gntuju(:,:,ic,ig)
          b1(:,:,ibatch) = zzero
          b2(:,:,ibatch) = zzero

          ! Loop for a single batch
          !DO ki = 1, nsize
          DO ki = 1, nstsvup

             !iband = k1 + ki - 1
             iband = spinupidx( ki )
             i = idxtranblhloc( iband, ikloc )
             ist1 = bmegqblh(1,i,ikloc)

             ! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
             b1(1:nmt,ki,ibatch) = DCONJG( wfsvmt1(1:nmt,ias,ispn,ist1) * &
                                            sfacgq(ig,ias) )

          END DO ! ki

       END DO ! ias
    END DO ! ig
    !$ACC END PARALLEL LOOP

!  END DO ! k1

#endif /* _OPENACC */

    RETURN
  END SUBROUTINE genmegqblh_fillbatch_acc

!==============================================================================
! Kernel 1b - CUDA version (not implemented)
! Fill in bgntuju and b1 (and zero b2) using CUDA C++ kernel
! Also outputs "translation table" for (ias,ig,iblock) to ibatch as batchidx

  SUBROUTINE genmegqblh_fillbatch_cuda( d_bgntuju, d_b1, d_b2, batchidx, &
                                        d_gntuju, d_sfacgq, d_wfsvmt1, &
                                        nmt, nstsvup, &
                                        iq, ikloc, ispn, spinupidx )
! BIND(C, NAME='genmegblh_fillbatch_cu')

    USE ISO_C_BINDING
    IMPLICIT NONE

    ! Device pointers
    TYPE(C_PTR) :: d_bgntuju, d_b1, d_b2, d_gntuju, d_sfacgq, d_wfsvmt1

    ! Other arguments
    INTEGER(KIND=C_INT), DIMENSION(:,:,:), INTENT(OUT) :: batchidx
    INTEGER(KIND=C_INT), INTENT(IN) :: nmt, nstsvup, iq, ikloc, ispn
    INTEGER(KIND=C_INT), DIMENSION(:), INTENT(IN) :: spinupidx

    ! CALL genmegblh_fillbatch_cu_()

    RETURN
  END SUBROUTINE genmegqblh_fillbatch_cuda

!==============================================================================
! Kernel 1c - Fallback mechanism for no GPUs
! Fill in bgntuju and b1 on CPU using OpenMP parallel do
! Also outputs "translation table" for (ias,ig,iblock) to ibatch as batchidx

  SUBROUTINE genmegqblh_fillbatch_omp( bgntuju, b1, b2, batchidx, &
                                     wfsvmt1, nmt, nstsvup, &
                                     iq, ikloc, ispn, spinupidx )
    USE modmain      ! for dc, natmtot
    USE mod_expigqr  ! for gntuju, bmegqblh
    USE mod_addons_q ! for sfacgq, ngq(iq)
    USE mod_nrkp     ! for spinor_ud

#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    ! Arguments
    COMPLEX(KIND=dz), DIMENSION(:,:,:), INTENT(OUT) :: bgntuju, b1, b2
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: batchidx ! Translation table for (ias,ig,iblock) to ibatch
    COMPLEX(KIND=dz), DIMENSION(:,:,:,:), INTENT(IN) :: wfsvmt1
    INTEGER, INTENT(IN) :: nstsvup, nmt, iq, ikloc, ispn
    INTEGER, DIMENSION(:), INTENT(IN) :: spinupidx
    ! Note: other ingredients (gntuju,sfacgq) are passed through modules

    ! Internal variables
    !INTEGER, PARAMETER :: nb = 64     ! Block size for ZGEMM batching
    !COMPLEX(KIND=dz), DIMENSION(nmt,nb) :: myb1      ! Local thread copy
    COMPLEX(KIND=dz), DIMENSION(nmt,nstsvup) :: myb1 ! Local thread copy
    INTEGER :: ibatch, nbatch, nblock      ! Batch index and number of batches
    INTEGER :: k1, k2, ki, nsize           ! Dummy variables for batching
    INTEGER :: iband, i, ist1, ic, ig, ias ! Data access and/or loop indices
    INTEGER :: tid

    ! Batching by block size nb for idxhiband
    ! TODO: Re-enable if needed for larger systems

    !  idxhiband = idxhibandblhloc(ikloc)
    !  ibatch = 0
    !  DO k1 = 1, idxhiband, nb
    !     k2 = MIN( idxhiband, k1+nb-1 )
    !     nsize = k2 - k1 + 1
    !     iblock = ???

    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
    !$OMP   PRIVATE(ig,ias,ic,ki,iband,i,ispn,ist1,myb1,ibatch)
    DO ig = 1, ngq(iq)
       DO ias = 1, natmtot

          ic = ias2ic(ias)

          ! ibatch = ???
          ! batchidx(ias,ig,iblock) = ibatch

          ibatch = (ig-1)*natmtot + ias
          batchidx(ias,ig,1) = ibatch

#ifdef _OPENMP
!!!DEBUG
          tid = omp_get_thread_num()
          WRITE(*,*) 'genmegqblh_fillbatch_omp: tid=', tid, ' ias=', ias, ' ig=', ig, ' ibatch=', ibatch
!!!DEBUG
#endif /* _OPENMP */

          ! Loop for a single batch
          myb1(:,:) = zzero
          !DO ki = 1, nsize
          DO ki = 1, nstsvup

             !iband = k1 + ki - 1
             iband = spinupidx( ki )
             i = idxtranblhloc( iband, ikloc )
             ist1 = bmegqblh(1,i,ikloc)

             ! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
             myb1(:,ki) = DCONJG( wfsvmt1(:,ias,ispn,ist1) * &
                                   sfacgq(ig,ias) )

          END DO ! ki

          !$OMP CRITICAL
          bgntuju(:,:,ibatch) = gntuju(:,:,ic,ig)
          b1(:,:,ibatch) = myb1(:,:)
          b2(:,:,ibatch) = zzero
          !$OMP END CRITICAL

       END DO ! ias
    END DO ! ig
    !$OMP END PARALLEL DO

!  END DO ! k1

    RETURN
  END SUBROUTINE genmegqblh_fillbatch_omp

!==============================================================================
! Kernel 2: Perform batched ZGEMM b2(:,:) = b1(:,:) x bgntuju(:,:)
!==============================================================================

  SUBROUTINE genmegqblh_batchzgemm( bgntuju, b1, b2, nmt, nstspin, nbatch )

#ifdef _MAGMA_
    USE mod_magma
#endif /* _MAGMA_ */

    IMPLICIT NONE

    ! (dummy) Arguments
    COMPLEX(KIND=dz), DIMENSION(:,:,:) :: bgntuju, b1, b2
    INTEGER :: nmt, nstspin, nbatch

  !-2a-------------------------------------------------------------------------
    IF( useacc .AND. usemagma ) THEN
  !----------------------------------------------------------------------------

       ! Perform batched ZGEMM on device using MAGMA
       CALL zgemm_batched_gpu_acc_magma( 'N', 'N', nmt, nstspin, nmt, &
                                         zone,  bgntuju(:,:,:), nmt, &
                                                b1(:,:,:),      nmt, &
                                         zzero, b2(:,:,:),      nmt, &
                                         nbatch )

       ! Synchronize with device
       CALL magma_queue_sync( queue )

       ! Clean up unneeded device copy
       !$ACC EXIT DATA DELETE( bgntuju, b1 )

  !-2b-------------------------------------------------------------------------
    !ELSE IF( usecuda .AND. usecublas )
  !----------------------------------------------------------------------------

       ! Perform batched ZGEMM on device using cuBLAS
       !CALL cublasZgemmBatched( blashandle, CUBLAS_OP_N, CUBLAS_OP_N, ... )

  !-2c-------------------------------------------------------------------------
    ELSE ! Fall back to CPU only using OpenMP
  !----------------------------------------------------------------------------

       ! Perform batched ZGEMM on CPU using OpenMP parallel do
       ! b2(1:nmt,1:nstsvup) = bgntuju(1:nmt,1:nmt) x b1(1:nmt,1:nstsv
       CALL zgemm_batched_omp( 'N', 'N', nmt, nstspin, nmt, &
                               zone,  bgntuju(:,:,:), nmt, &
                                      b1(:,:,:),      nmt, &
                               zzero, b2(:,:,:),      nmt, &
                               nbatch )

  !----------------------------------------------------------------------------
    END IF ! CPU/GPU method
  !----------------------------------------------------------------------------

    RETURN
  END SUBROUTINE genmegqblh_batchzgemm

!==============================================================================
! Kernel 3: Save results to wftmp1mt and transfer back to CPU (for now)
!==============================================================================

  SUBROUTINE genmegqblh_fillresult( b2, wftmp1mt, &
                                    iq, nmt, nstspin, spinstidx, batchidx )
    IMPLICIT NONE

    ! (dummy) Arguments
    COMPLEX(KIND=dz), DIMENSION(:,:,:) :: b2
    COMPLEX(KIND=dz), DIMENSION(:,:,:,:) :: wftmp1mt
    INTEGER, DIMENSION(:,:,:) :: batchidx    
    INTEGER, DIMENSION(:) :: spinstidx
    INTEGER :: iq, nmt, nstspin

  !-3a-------------------------------------------------------------------------
    IF( useacc .AND. usemagma ) THEN
  !----------------------------------------------------------------------------

       ! Fill in wftmp1mt on device
       CALL genmegqblh_fillresult_acc( b2, wftmp1mt, &
                                       iq, nmt, nstspin, spinstidx, batchidx )

       ! Transfer data to CPU
       !$ACC UPDATE SELF( wftmp1mt )  

       ! Clean up (for now)
       !$ACC EXIT DATA DELETE( b2, nstspin, spinstidx, batchidx )

  !-3b-------------------------------------------------------------------------
     !ELSE IF( usecuda .AND. usecublas )
  !----------------------------------------------------------------------------

       !CALL genmegqblh_fillresult_cuda( d_b2, d_wfsvmt1mt, nmt, nstsvup, spinstidx, batchidx )

       !CALL cudaMemcpy( wftmp1mt, d_wftmp1mt, cudaMemcpyDeviceToHost )

       !CALL cudaFree ...

  !-3c-------------------------------------------------------------------------
    ELSE ! Fall back to CPU only using OpenMP
  !----------------------------------------------------------------------------

       ! Save results to wftmp1mt
       CALL genmegqblh_fillresult_omp( b2, wftmp1mt, &
                                       iq, nmt, nstspin, spinstidx, batchidx )

  !----------------------------------------------------------------------------
    END IF ! CPU/GPU method
  !----------------------------------------------------------------------------

    RETURN
  END SUBROUTINE genmegqblh_fillresult

!==============================================================================
! Kernel 3a - OpenACC version
! Fill in wftmp1mt using OpenACC parallel loop
! TODO: Write directly to wftmp1 after interstitial part is ported

  SUBROUTINE genmegqblh_fillresult_acc( b2, wftmp1mt, &
                                        iq, nmt, nstsvup, spinupidx, batchidx )
    USE modmain      ! for dc, natmtot
    USE mod_addons_q ! for ngq(iq)
    USE mod_gpu      ! for devnum

#ifdef _OPENACC
    USE openacc
#endif /* _OPENACC */

    IMPLICIT NONE

    ! Arguments
    COMPLEX(KIND=dz), DIMENSION(:,:,:), INTENT(IN) :: b2
    COMPLEX(KIND=dz), DIMENSION(:,:,:,:), INTENT(OUT) :: wftmp1mt
    INTEGER, INTENT(IN) :: nmt, nstsvup, iq
    INTEGER, DIMENSION(:), INTENT(IN) :: spinupidx
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: batchidx

#ifdef _OPENACC

    ! Internal variables
    !INTEGER, PARAMETER :: nb = 64
    INTEGER :: nblock
    INTEGER :: k1, k2, ki, ist, ig, ias, ibatch

    ! Blocked version
!  nblock = SIZE(batchidx, 3)
!  DO iblock = 1, nblock
!     k1 = (ki-1)*nb + 1
!     IF( iblock == nblock ) THEN
!        k2 = idxhiband
!     ELSE
!        k2 = ki*nb
!     END IF

    ! These are contiguous, for now
    ! TODO: check behavior of spinor_ud for other systems
    k1 = spinupidx(1)
    k2 = spinupidx(nstsvup)

    ! For consistency
    ! TODO: generalize for AMD GPUs
    CALL acc_set_device_num( devnum, acc_device_nvidia )

    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC   PRESENT( b2, wftmp1mt, ngq(iq), nstsvup, spinupidx, batchidx ) &
    !$ACC   CREATE(ibatch) COPYIN(k1,k2,natmtot)
    DO ig = 1, ngq(iq)
       DO ias = 1, natmtot

          !ibatch = batchidx(ias,ig,iblock)
          ibatch = batchidx(ias,ig,1)

           wftmp1mt(1:nmt,k1:k2,ias,ig) = b2(1:nmt,1:nstsvup,ibatch)

           ! If non-contiguous
           !DO ist = 1, nstsvup
           !   ki = spinupidx(ist)
           !   wftmp1mt(1:nmt,ki,ias,ig) = b2(1:nmt,ist,ibatch)
           !END DO ! ist
        
        END DO ! ias
     END DO ! ig

!  END DO ! iblock

#endif /* _OPENACC */

     RETURN
   END SUBROUTINE genmegqblh_fillresult_acc

!==============================================================================
! Kernel 3b - CUDA version (not implemented)
! Fill in wftmp1mt using CUDA C++ kernel
! TODO: Write directly to wftmp1 after interstitial part is ported

   SUBROUTINE genmegqblh_fillresult_cuda( d_b2, d_wftmp1mt, &
                                          iq, nmt, nstsvup, spinupidx, batchidx )
! BIND(C, NAME='genmegblh_fillresult_cu')

     USE ISO_C_BINDING
     IMPLICIT NONE

     ! Arguments
     TYPE(C_PTR) :: d_b2, d_wftmp1mt ! Device pointers
     INTEGER(KIND=C_INT), DIMENSION(:,:,:), INTENT(OUT) :: batchidx
     INTEGER(KIND=C_INT), INTENT(IN) :: iq, nmt, nstsvup
     INTEGER(KIND=C_INT), DIMENSION(nstsvup), INTENT(IN) :: spinupidx

     ! CALL genmegblh_fillresult_cu_()

     RETURN
   END SUBROUTINE genmegqblh_fillresult_cuda

!==============================================================================
! Kernel 3c - Fallback mechanism for no GPUs
! Fallback mechanism for no GPUs:
! Fill in wftmp1mt on CPU using OpenMP parallel do

   SUBROUTINE genmegqblh_fillresult_omp( b2, wftmp1mt, &
                                        iq, nmt, nstsvup, spinupidx, batchidx )
     USE modmain
     USE mod_addons_q ! for ngq(iq)

#ifdef _OPENMP
     USE omp_lib
#endif /* _OPENMP */

     IMPLICIT NONE

     ! Arguments
     COMPLEX(KIND=dz), DIMENSION(:,:,:), INTENT(IN) :: b2
     COMPLEX(KIND=dz), DIMENSION(:,:,:,:), INTENT(OUT) :: wftmp1mt
     INTEGER, INTENT(IN) :: iq, nmt, nstsvup
     INTEGER, DIMENSION(:), INTENT(IN) :: spinupidx
     INTEGER, DIMENSION(:,:,:), INTENT(IN) :: batchidx

     ! Internal variables
     !INTEGER, PARAMETER :: nb = 64
     INTEGER :: nblock
     INTEGER :: k1, k2, ki, ist, ig, ias, ibatch
     INTEGER :: tid

     ! Blocked version
!  nblock = SIZE(batchidx, 3)
!  DO iblock = 1, nblock
!     k1 = (ki-1)*nb + 1
!     IF( iblock == nblock ) THEN
!        k2 = idxhiband
!     ELSE
!        k2 = ki*nb
!     END IF

     ! These are contiguous, for now
     ! TODO: check behavior of spinor_ud for other systems
     k1 = spinupidx(1)
     k2 = spinupidx(nstsvup)

     !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
     !$OMP   PRIVATE(ibatch,k1,k2,ki)
     DO ig = 1, ngq(iq)
        DO ias = 1, natmtot

           !ibatch = batchidx(ias,ig,iblock)
           ibatch = batchidx(ias,ig,1)

#ifdef _OPENMP
!!!DEBUG
           tid = omp_get_thread_num()
           WRITE(*,*) 'genmegqblh_fillresult_omp: tid=', tid, &
                      ' ias=', ias, ' ig=', ig, ' ibatch=', ibatch, ' k1=', k1, ' k2=', k2
!!!DEBUG
#endif /* _OPENMP */

           !$OMP CRITICAL
           wftmp1mt(1:nmt,k1:k2,ias,ig) = b2(1:nmt,1:nstsvup,ibatch)
           !$OMP END CRITICAL

           ! If non-contiguous
           !DO ist = 1, nstsvup
           !   ki = spinupidx(ist)
           !   wftmp1mt(1:nmt,ki,ias,ig) = b2(1:nmt,ist,ibatch)
           !END DO ! ist
        
        END DO ! ias
     END DO ! ig

!  END DO ! iblock

     RETURN
   END SUBROUTINE genmegqblh_fillresult_omp

!==============================================================================

END MODULE mod_genmegqblh_gpu
