subroutine genmegqblh(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
  &wfsvit1,wfsvit2)
use modmain
use mod_addons_q
use mod_nrkp
use mod_expigqr
USE mod_gpu

implicit none
integer, intent(in) :: iq
integer, intent(in) :: ikloc
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvmt1(lmmaxapw*nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(lmmaxapw,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)

integer wfsize
integer ivg1(3)
integer i,j,ik,jk,igkq,n1,ispn1,ispn2,ist1,ist2,ic, j1
integer ig,ig1,ig2,ias,ifg,ir
logical l1
complex(8), allocatable :: wftmp1(:,:)
complex(8), allocatable :: wftmp2(:,:)
complex(8), allocatable :: wfir1(:)

!--begin Convert to true ZGEMM
  INTEGER :: nmt                    ! Number of muffin-tin elements
  !INTEGER, PARAMETER :: nb = 64     ! Block size for ZGEMM batching
  INTEGER :: ibatch, nbatch, nblock ! Batch index and number of batches
  INTEGER :: k1, k2, ki, nsize      ! Dummy variables for batching
  INTEGER :: nstspin                ! Number of 2nd-variational states per spin

! Blocked version
!  COMPLEX(KIND=dc), DIMENSION( lmmaxapw*nufrmax, nb, natmtot, ngq(iq) ) :: wftmp1mt
  ! Allocatable due to 3rd dimension (nbatch) undetermined until runtime
!  COMPLEX(KIND=dc), DIMENSION(:,:,:), ALLOCATABLE :: b1, b2, bgntuju
!  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: batchidx

! Unblocked version
  COMPLEX(KIND=dc), DIMENSION( lmmaxapw*nufrmax, idxhibandblhloc(ikloc), &
                               natmtot, ngq(iq) ) :: wftmp1mt
  COMPLEX(KIND=dc), DIMENSION( lmmaxapw*nufrmax, lmmaxapw*nufrmax, &
                               natmtot*ngq(iq) ) :: bgntuju
  ! b1 and b2 are allocatable due to 2nd dimension (nstspin) undetermined until runtime
  ! TODO: move this to mod_expigqr, maybe? So these can be automatic arrays too?
  COMPLEX(KIND=dc), DIMENSION(:,:,:), ALLOCATABLE :: b1, b2
  INTEGER, DIMENSION(natmtot,ngq(iq),1) :: batchidx

  ! Table of spin-up/dn states (replaces l1 check)
  INTEGER, DIMENSION(:), ALLOCATABLE :: spinstidx 

!--end Convert to true ZGEMM

#if defined(_DEBUG_bmegqblh_) || defined(_DEBUG_megqblh_)
  INTEGER :: dbgcnt1, dbgcnt2, dbgunit1, dbgunit2
#endif /* _DEBUG_bmegqblh_ || _DEBUG_megqblh_ */

INTEGER :: idxhiband, iband, ntran, idxtran
EXTERNAL :: zcopy
INTEGER, EXTERNAL :: genmegqblh_countspinup
INTEGER, EXTERNAL :: genmegqblh_countspindn

wfsize=lmmaxapw*nufrmax*natmtot+ngknr2
allocate(wftmp1(wfsize,ngq(iq))) ! TODO: Change dimensions appropriately
allocate(wftmp2(wfsize,nstsv))   ! TODO: Check contiguity of ZCOPY transfers
allocate(wfir1(ngrtot))
call papi_timer_start(pt_megqblh)

!$ACC ENTER DATA COPYIN( wfsvmt1 ) CREATE( wftmp1, wftmp1mt )

! Note: List of OpenACC variables that are already in device memory 
!       due to inheritance from mod_expigqr::genmegq() :
!         sfacgq, gntuju, bmegqblh, idxhibandblhloc, idxtranblhloc,
!         spinor_ud, ngq(iq), ias2ic

! global k-point
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! jk=k+q-G_q
jk=idxkq(1,ik)
! G_q vector 
igkq=idxkq(2,ik)

#ifdef _DEBUG_bmegqblh_
  dbgunit1 = 1000 + iproc ! Make sure this matches the definition in mod_expigqr::genmegq()
  dbgcnt1 = 1
  WRITE( dbgunit1, '(A,I3,A,I5)' ) 'nmegqblh(ikloc=', ikloc, ') = ', nmegqblh(ikloc)
#endif /* _DEBUG_bmegqblh_ */

#ifdef _DEBUG_megqblh_
  dbgunit2 = 2000 + iproc ! Make sure this matches the definition in mod_expigqr::genmegq()
  !dbgcnt2 = 1 ! Not needed anymore since we have ibatch now
#endif

do ispn1=1,nspinor

  ! expigqr22 is always 1, for now (see mod_expigqr)
  if (expigqr22.eq.1) ispn2=ispn1

  ! Convert to the corresponding ist1 loop in getmeidx() line 56-149
  ! as stored into idxhibandblh(ikloc=1:nkptnr) at getmeidx() line 155,
  ! skipping as necessary (with a warning message... it should NOT happen!)
  idxhiband = idxhibandblhloc(ikloc)
  IF( idxhiband == 0 ) THEN
     ! Unit 151 is either 'CRPA.OUT' or 'RESPONSE.OUT'
     WRITE(151, '( "Warning[genmegqblh]: highest band is zero for iq=", &
                 &I6, " ikloc=", I6, " ispn1=", I1 )' ) iq, ikloc, ispn1
     CYCLE
  END IF

!--begin Convert to true ZGEMM

  call timer_start(3)
  call papi_timer_start(pt_megqblh_mt)

  ! Note that the loop order has been switched
  ! such that iband loop is now the innermost loop

  ! Number of muffin-tin elements
  nmt = lmmaxapw*nufrmax

  ! Number of batches, blocked version
  !nblock = CEILING( REAL(idxhiband)/REAL(nb) )
  !nbatch = ngq(iq) * natmtot * nblock
  
  ! Number of batches, unblocked version
  nbatch = ngq(iq) * natmtot

  ! Blocked version
  !ALLOCATE( bgntuju( nmt, nmt, nbatch ))
  !ALLOCATE( b1( nmt, nb, nbatch ))
  !ALLOCATE( b2( nmt, nb, nbatch ))
  !ALLOCATE( batchidx( natmtot, ngq(iq), nblock ))

  ALLOCATE( spinstidx( nstsv ))
  ! Count spin up states for this particular k-vector (replaces l1 check)
  ! Note: after the function call, spinupidx will be reallocated to (1:nstspin) 
  IF( ispn1 == 1 ) THEN
     ! Spin up
     spinstidx(1:nstspin) = genmegqblh_countspinup( ikloc, nstspin, spinstidx )
  ELSE
     ! Spin down (never executed if spinpol = .FALSE. )
     spinstidx(1:nstspin) = genmegqblh_countspindn( ikloc, nstspin, spinstidx )
  END IF
  !$ACC ENTER DATA COPYIN( nstspin, spinstidx )

  ! Unblocked version
  ALLOCATE( b1( nmt, nstspin, nbatch ))
  ALLOCATE( b2( nmt, nstspin, nbatch ))

!------------------------------------------------------------------------------
  IF( useacc .AND. usemagma ) THEN
!------------------------------------------------------------------------------

     ! Fill in bgntuju and b1 on device
     !$ACC ENTER DATA CREATE( bgntuju, b1, b2, batchidx )
     CALL genmegqblh_fillbatch_acc( bgntuju, b1, b2, batchidx, &
                                    wfsvmt1, nmt, nstspin, &
                                    iq, ikloc, ispn1, spinstidx )

     ! Perform batched ZGEMM on device using MAGMA
     CALL zgemm_batched_gpu_acc_magma( 'N', 'N', nmt, nstspin, nmt, &
                                        zone,  bgntuju(:,:,:), nmt, &
                                               b1(:,:,:),      nmt, &
                                        zzero, b2(:,:,:),      nmt, &
                                        nbatch )

     ! Save results to wftmp1mt and transfer back to CPU (for now)
     CALL genmegqblh_fillresult_acc( b2, wftmp1mt, &
                                     iq, nmt, nstspin, spinstidx, batchidx )
     !$ACC UPDATE SELF( wftmp1mt )  

     ! Clean up (for now)
     !$ACC EXIT DATA DELETE( bgntuju, b1, b2, nstspin, spinstidx, batchidx )

!------------------------------------------------------------------------------
  !ELSE IF( usecuda .AND. usecublas )
!------------------------------------------------------------------------------

     !CALL cudaMemcpy( d_wfsvmt1,   wfsvmt1,   cudaMemcpyHostToDevice )
     !CALL cudaMemcpy( d_sfacgq,    sfacgq,    cudaMemcpyHostToDevice )
     !CALL cudaMemcpy( d_gntuju,    gntuju,    cudaMemcpyHostToDevice )

     !CALL cudaMalloc( d_bgntuju, ... )
     !CALL cudaMalloc( d_b1, ... )
     !CALL cudaMalloc( d_b2, ... )

     !CALL genmegqblh_fillbatch_cuda( d_bgntuju, d_b1, d_b2, batchidx, &
     !                                d_gntuju, d_sfacgq, d_wfsvmt1, &
     !                                nmt, nstspin, iq, ikloc, ispn, spinstidx )

     !CALL cublasZgemmBatched( blashandle, CUBLAS_OP_N, CUBLAS_OP_N, ... )

     !CALL genmegqblh_fillresult_cuda( d_b2, d_wfsvmt1mt, nmt, nstsvup, spinstidx, batchidx )

     !CALL cudaMemcpy( wftmp1mt, d_wftmp1mt, cudaMemcpyDeviceToHost )

     !CALL cudaFree ...

!------------------------------------------------------------------------------
  ELSE ! Fall back to CPU only using OpenMP
!------------------------------------------------------------------------------

     ! Fill in bgntuju and b1 on CPU
     CALL genmegqblh_fillbatch_omp( bgntuju, b1, b2, batchidx, &
                                    wfsvmt1, nmt, nstspin, &
                                    iq, ikloc, ispn1, spinstidx )

     ! Original code (retained for historical purpose)
     !do j=1,ngntuju(ic,ig)
     !  b2(igntuju(2,j,ic,ig))=b2(igntuju(2,j,ic,ig))+&
     !    &b1(igntuju(1,j,ic,ig))*gntuju(j,ic,ig)
     !enddo

     ! Perform batched ZGEMM on CPU using OpenMP parallel do
     ! b2(1:nmt,1:nstsvup) = bgntuju(1:nmt,1:nmt) x b1(1:nmt,1:nstsv
     CALL zgemm_batched_omp( 'N', 'N', nmt, nstspin, nmt, &
                              zone,  bgntuju(:,:,:), nmt, &
                                     b1(:,:,:),      nmt, &
                              zzero, b2(:,:,:),      nmt, &
                              nbatch )

     ! Save results to wftmp1mt
     CALL genmegqblh_fillresult_omp( b2, wftmp1mt, &
                                     iq, nmt, nstspin, spinstidx, batchidx )

  END IF ! CPU/GPU method

  ! Clean up
  !DEALLOCATE( bgntuju )
  DEALLOCATE( b1 )
  DEALLOCATE( b2 )
  !DEALLOCATE( batchidx )
  DEALLOCATE( spinstidx )

  call timer_stop(3)
  call papi_timer_stop(pt_megqblh_mt)

  ! Start the bounded do loop for each band
  ! TODO: Complete removal of l1 check
  DO iband = 1, idxhiband

! left <bra| state 
     wftmp1=zzero

     ! The starting point of the index "i" for accessing bmegqblh(:,i,:)
     ! for each iband and ikloc was stored as idxtranblhloc
     i = idxtranblhloc( iband, ikloc )
     ist1 = bmegqblh(1,i,ikloc)

     ! Same as above
     ! TODO: Complete removal of l1 check
     l1=.true.
     if (spinpol) then
        if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
     endif

     ! Note: wftmp1 combines the muffin-tin and interstitial parts
     !       for each band, to prepare for the second ZGEMM below
     !       Complete removal of wftmp1mt is impossible until
     !       interstitial part also ported to GPU (cuFFT with fallback to FFTW)
     if (l1) then
        DO ig = 1, ngq(iq)
           DO ias = 1, natmtot
              wftmp1( (ias-1)*nmt+1:ias*nmt, ig ) = wftmp1mt( 1:nmt, iband, ias, ig )
           END DO ! ias
        END DO ! ig

!--end Convert to true ZGEMM

! interstitial part
      call papi_timer_start(pt_megqblh_it)
      call timer_start(4)
      wfir1=zzero
      do ig1=1,ngknr1
        ifg=igfft(igkignr1(ig1))
        wfir1(ifg)=wfsvit1(ig1,ispn1,ist1)
      enddo
      call zfftifc(3,ngrid,1,wfir1)
      do ir=1,ngrtot
        wfir1(ir)=wfir1(ir)*cfunir(ir)
      enddo
      call zfftifc(3,ngrid,-1,wfir1)
      do ig=1,ngq(iq)
        do ig2=1,ngknr2
! G1=G2-G-Gkq
          ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,igqig(ig,iq))-ivg(:,igkq)
          ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
          wftmp1(lmmaxapw*nufrmax*natmtot+ig2,ig)=dconjg(wfir1(ifg))
        enddo
      enddo
      call timer_stop(4)      
      call papi_timer_stop(pt_megqblh_it)

    endif !l1

    call timer_start(5)

#ifdef _DEBUG_bmegqblh_
IF( ntran > 0 ) THEN
  WRITE( dbgunit1, '(7(1X,I5))' ) dbgcnt1, ikloc, iq, iband, i, ntran, i+ntran-1
  dbgcnt1 = dbgcnt1 + 1
END IF
#endif /* _DEBUG_bmegqblh_ */

  ! Load number of matching |ist2=n'> ket states for each <ist1=n| bra
  IF( ltranconst ) THEN
     ntran = ntranblhloc(ikloc)
  ELSE
     ! Note: seeing the pattern, this shouldn't happen, but just in case
     IF( iband == idxhiband ) THEN
        ntran = nmegqblh(ikloc) - idxtranblhloc(idxhiband,ikloc) + 1
     ELSE
        ntran = idxtranblhloc(iband+1,ikloc) - idxtranblhloc(iband,ikloc)
     END IF
  END IF ! ltranconst

! collect right |ket> states into matrix wftmp2
    ! Note: ntran can be zero (zero-trip loop)
    DO n1 = 1, ntran

       ist2 = bmegqblh(2,i+n1-1,ikloc) ! Now n1 starts from 1 instead of 0

       ! Following Ed's advice, use ZCOPY() from BLAS instead of memcopy
       ! TODO: check whether it's better to transfer all in one go
       !       or overlap computation & data movement

       ! Muffin tin
       CALL zcopy( lmmaxapw*nufrmax*natmtot, &
                   wfsvmt2(1,1,1,ispn2,ist2), 1, &
                   wftmp2(1,n1), 1 )
       ! Interstitial
       CALL zcopy( ngknr2, &
                   wfsvit2(1,ispn2,ist2), 1, &
                   wftmp2(lmmaxapw*nufrmax*natmtot+1,n1), 1 )

    END DO ! n1; replaced do while loop (i+n1) <= nmegqblh(ikloc)

! update several matrix elements by doing matrix*matrix operation
!  me(ib,ig)=wftmp2(ig2,ib)^{T}*wftmp1(ig2,ig)

    ! This particular ZGEMM() call corresponds with line 9 of Algorithm 2
    ! in the Gordon Bell paper
    CALL zgemm( 'T', 'N', ntran, ngq(iq), wfsize, zone, &
                wftmp2, wfsize, wftmp1, wfsize, zone, &
                megqblh(i,1,ikloc), nstsv*nstsv )

    ! No need to add n1 to i anymore to move on to the next <nk| bra
    ! since it is already stored as ntranblhloc

    CALL timer_stop(5) ! Same as before

 END DO ! iband; replaces do while loop i <= nmegqblh(ikloc)

#ifdef _DEBUG_bmegqblh_
     WRITE( dbgunit1, '(A,I3)') 'highest band = ', idxhiband
#endif /* _DEBUG_bmegqblh_ */

!--end Convert do while into bounded do loop

enddo !ispn

!$ACC EXIT DATA DELETE( wftmp1, wftmp1mt )

deallocate(wftmp1)
deallocate(wftmp2)
deallocate(wfir1)

!$ACC EXIT DATA DELETE( wfsvmt1 )

call papi_timer_stop(pt_megqblh)

return
end subroutine genmegqblh
