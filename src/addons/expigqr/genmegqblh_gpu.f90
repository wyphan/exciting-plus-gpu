!==============================================================================
! Counts how many 2nd-variational states are spin up,
! and returns a list of such states as stateidx
! For now, the contents of stateidx should be consecutive
! Note: uses Fortran 2003 intrinsic MOVE_ALLOC, following
! https://stackoverflow.com/questions/8264336/how-to-get-priorly-unknown-array-as-the-output-of-a-function-in-fortran/8265903#8265903
! TODO: Perform on device? (rewrite using OpenACC?)

FUNCTION genmegqblh_countspinup( ikloc, nstsvup, spinupidxold ) RESULT( spinupidx )

  USE modmain     ! for nstsv
  USE mod_nrkp    ! for spinor_ud
  USE mod_expigqr ! for idxhibandblhloc(:), idxtranblhloc(:,:)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: ikloc
  INTEGER, INTENT(OUT) :: nstsvup
  INTEGER, DIMENSION(:), INTENT(INOUT) :: spinupidxold
  INTEGER, DIMENSION(:), ALLOCATABLE :: spinupidx

  ! Internal variables
  INTEGER :: iband, i, ist1
  INTEGER, DIMENSION(:), ALLOCATABLE :: temp

  IF( spinpol ) THEN

     nstsvup = 0
     spinupidxold(:) = 0
     DO iband = 1, idxhibandblhloc(ikloc)

        i = idxtranblhloc( iband, ikloc )
        ist1 = bmegqblh(1,i,ikloc)

        IF( spinor_ud(1,i,ikloc) == 1 .AND. spinor_ud(2,i,ikloc) == 0 ) THEN
           nstsvup = nstsvup + 1
           spinupidxold(nstsvup) = i
        END IF

     END DO ! iband

  ELSE

     ! If spinpol is .FALSE. there is only one spin projection
     nstsvup = nstsv
     DO iband = 1, idxhibandblhloc(ikloc)
        spinupidxold(iband) = iband
     END DO ! iband

  END IF ! spinpol

!!!DEBUG
  WRITE(*,*) 'genmegqblh_countspinup: ikloc=', ikloc, ' nstsvup=', nstsvup
!!!DEBUG

  ! Now that we know the extent of the new array, copy array and move_alloc
  ALLOCATE( temp( nstsvup ) )
  temp(1:nstsvup) = spinupidxold(1:nstsvup)
  CALL MOVE_ALLOC( temp, spinupidx )

END FUNCTION genmegqblh_countspinup

!==============================================================================
! Same as above, but for spin down
! Note: uses Fortran 2003 intrinsic MOVE_ALLOC
! TODO: Merge into the above function?

FUNCTION genmegqblh_countspindn( ikloc, nstsvdn, spindnidxold ) RESULT( spindnidx )

  USE modmain     ! for nstsv
  USE mod_nrkp    ! for spinor_ud
  USE mod_expigqr ! for idxhibandblhloc(:), idxtranblhloc(:,:)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: ikloc
  INTEGER, INTENT(OUT) :: nstsvdn
  INTEGER, DIMENSION(:), INTENT(INOUT) :: spindnidxold
  INTEGER, DIMENSION(:), ALLOCATABLE :: spindnidx

  ! Internal variables
  INTEGER :: iband, i, ist1
  INTEGER, DIMENSION(:), ALLOCATABLE :: temp

  IF( spinpol ) THEN

     nstsvdn = 0
     spindnidxold(:) = 0
     DO iband = 1, idxhibandblhloc(ikloc)

        i = idxtranblhloc( iband, ikloc )
        ist1 = bmegqblh(1,i,ikloc)

        ! Other than the variable names, this is the only change from above, TBH
        IF( spinor_ud(1,i,ikloc) == 0 .AND. spinor_ud(2,i,ikloc) == 1 ) THEN
           nstsvdn = nstsvdn + 1
           spindnidxold(nstsvdn) = i
        END IF

     END DO ! iband

  END IF ! spinpol
 ! If spinpol is .FALSE. there is only one spin projection
 ! (Technically, this function shouldn't even be executed!)

!!!DEBUG
! This output should NOT appear when spinpol is .FALSE.
! When spinpol == .TRUE., nstsvdn should be identical to nstsvup
  WRITE(*,*) 'genmegqblh_countspindn: ikloc=', ikloc, ' nstsvdn=', nstsvdn
!!!DEBUG

  ! Now that we know the extent of the new array, copy array and move_alloc
  ALLOCATE( temp( nstsvdn ) )
  temp(1:nstsvdn) = spindnidxold(1:nstsvdn)
  CALL MOVE_ALLOC( temp, spindnidx )

  RETURN
END FUNCTION genmegqblh_countspindn

!==============================================================================
! Fill in bgntuju and b1 (and zero b2) using OpenACC
! Also outputs "translation table" for (ias,ig,iblock) to ibatch as batchidx

SUBROUTINE genmegqblh_fillbatch_acc( bgntuju, b1, b2, batchidx, &
                                     wfsvmt1, nmt, nstsvup, &
                                     iq, ikloc, ispn, spinupidx )

  USE modmain      ! for dc, natmtot
  USE mod_expigqr  ! for gntuju, bmegqblh
  USE mod_addons_q ! for sfacgq, ngq(iq)
  USE mod_nrkp     ! for spinor_ud
  USE mod_gpu      ! for devnum

#ifdef _OPENACC
  USE openacc
#endif /* _OPENACC */

  IMPLICIT NONE

  ! Arguments
  COMPLEX(KIND=dc), DIMENSION(:,:,:), INTENT(OUT) :: bgntuju, b1, b2
  INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: batchidx
  COMPLEX(KIND=dc), DIMENSION(:,:,:,:), INTENT(IN) :: wfsvmt1
  INTEGER, INTENT(IN) :: nmt, nstsvup, iq, ikloc, ispn
  INTEGER, DIMENSION(nstsvup), INTENT(IN) :: spinupidx
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
! Fill in bgntuju and b1 (and zero b2) using CUDA C++ kernel (not implemented)
! Also outputs "translation table" for (ias,ig,iblock) to ibatch as batchidx

SUBROUTINE genmegqblh_fillbatch_cuda( d_bgntuju, d_b1, d_b2, batchidx, &
                                      d_gntuju, d_sfacgq, d_wfsvmt1, &
                                      nmt, nstsvup, iq, ikloc, ispn, spinupidx )
! BIND(C, NAME='genmegblh_fillbatch_cu')

  USE ISO_C_BINDING
  IMPLICIT NONE

  ! Device pointers
  TYPE(C_PTR) :: d_bgntuju, d_b1, d_b2, d_gntuju, d_sfacgq, d_wfsvmt1

  ! Other arguments
  INTEGER(KIND=C_INT), DIMENSION(:,:,:), INTENT(OUT) :: batchidx
  INTEGER(KIND=C_INT), INTENT(IN) :: nmt, nstsvup, iq, ikloc, ispn
  INTEGER(KIND=C_INT), DIMENSION(nstsvup), INTENT(IN) :: spinupidx
  ! Note: other ingredients (gntuju,sfacgq) are passed through modules

  RETURN
END SUBROUTINE genmegqblh_fillbatch_cuda

!==============================================================================
! Fallback mechanism for no GPUs:
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
  COMPLEX(KIND=dc), DIMENSION(:,:,:), INTENT(OUT) :: bgntuju, b1, b2
  INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: batchidx ! Translation table for (ias,ig,iblock) to ibatch
  COMPLEX(KIND=dc), DIMENSION(:,:,:,:), INTENT(IN) :: wfsvmt1
  INTEGER, INTENT(IN) :: nstsvup, nmt, iq, ikloc, ispn
  INTEGER, DIMENSION(nstsvup), INTENT(IN) :: spinupidx
  ! Note: other ingredients (gntuju,sfacgq) are passed through modules

  ! Internal variables
  !INTEGER, PARAMETER :: nb = 64     ! Block size for ZGEMM batching
  !COMPLEX(KIND=dc), DIMENSION(nmt,nb) :: myb1      ! Local thread copy
  COMPLEX(KIND=dc), DIMENSION(nmt,nstsvup) :: myb1 ! Local thread copy
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
! Fill in wftmp1mt using OpenACC
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
  COMPLEX(KIND=dc), DIMENSION(:,:,:), INTENT(IN) :: b2
  COMPLEX(KIND=dc), DIMENSION(:,:,:,:), INTENT(OUT) :: wftmp1mt
  INTEGER, INTENT(IN) :: nmt, nstsvup, iq
  INTEGER, DIMENSION(nstsvup), INTENT(IN) :: spinupidx
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
! Fill in wftmp1mt using CUDA C++ kernel (not implemented)
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

  RETURN
END SUBROUTINE genmegqblh_fillresult_cuda

!==============================================================================
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
  COMPLEX(KIND=dc), DIMENSION(:,:,:), INTENT(IN) :: b2
  COMPLEX(KIND=dc), DIMENSION(:,:,:,:), INTENT(OUT) :: wftmp1mt
  INTEGER, INTENT(IN) :: iq, nmt, nstsvup
  INTEGER, DIMENSION(nstsvup), INTENT(IN) :: spinupidx
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
