subroutine genmegqblh(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
  &wfsvit1,wfsvit2)
use modmain
use mod_addons_q
use mod_nrkp
use mod_expigqr

#ifdef _MAGMA_
USE mod_magma
#endif /* _MAGMA_ */

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
  INTEGER, PARAMETER :: nb = 64     ! Block size for ZGEMM batching
  INTEGER :: ibatch, nbatch, nblock ! Batch index and number of batches
  INTEGER :: k1, k2, ki, nsize      ! Dummy variables for batching
  COMPLEX(KIND((0.D0,1.D0))), DIMENSION(:,:,:), ALLOCATABLE :: b1, b2, bgntuju
  !COMPLEX(KIND((0.D0,1.D0))), DIMENSION( lmmaxapw*nufrmax, nb ) :: b1, b2
  COMPLEX(KIND((0.D0,1.D0))), &
    DIMENSION( lmmaxapw*nufrmax, nb, natmtot, ngq(iq) ) :: wftmp1mt
!--end Convert to true ZGEMM

#if defined(_DEBUG_bmegqblh_) || defined(_DEBUG_megqblh_)
  INTEGER :: dbgcnt1, dbgcnt2, dbgunit1, dbgunit2
#endif /* _DEBUG_bmegqblh_ || _DEBUG_megqblh_ */

INTEGER :: idxhiband, iband, ntran, idxtran
EXTERNAL :: zcopy

wfsize=lmmaxapw*nufrmax*natmtot+ngknr2
allocate(wftmp1(wfsize,ngq(iq))) ! TODO: Change dimensions appropriately
allocate(wftmp2(wfsize,nstsv))   ! TODO: Check contiguity of ZCOPY transfers
allocate(wfir1(ngrtot))
call papi_timer_start(pt_megqblh)

!$ACC ENTER DATA COPYIN( wfsvmt1 )

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
  if (expigqr22.eq.1) ispn2=ispn1

  ! Convert to the corresponding ist1 loop in getmeidx() line 56-149
  ! as stored into idxhibandblh(ikloc=1:nkptnr) at getmeidx() line 155,
  ! skipping as necessary (with a warning message... it should NOT happen!)
  idxhiband = idxhibandblhloc(ikloc)
  IF( idxhiband == 0 ) THEN
     ! Unit 151 is either 'CRPA.OUT' or 'RESPONSE.OUT'i
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

  ! Number of batches
  nblock = CEILING( REAL(idxhiband)/REAL(nb) )
  nbatch = ngq(iq) * natmtot * nblock
  ALLOCATE( bgntuju( nmt, nmt, nbatch ) )
  ALLOCATE( b1( nmt, nb, nbatch ) )
  ALLOCATE( b2( nmt, nb, nbatch ) )
#ifdef _OPENACC
  !$ACC ENTER DATA CREATE( bgntuju, b1, b2 )
#endif /* _OPENACC */

  ! Batching by block size nb and over ig and ias loops
  ibatch = 0
  DO k1 = 1, idxhiband, nb
     k2 = MIN( idxhiband, k1+nb-1 )
     nsize = k2 - k1 + 1

#ifdef _OPENACC
     ! Fill in bgntuju and b1 on device
     !$ACC PARALLEL LOOP COLLAPSE(2) &
     !$ACC   PRESENT( bgntuju, b1, b2, gntuju, sfacgq, wfsvmt1, &
     !$ACC            bmegqblh(:,:,ikloc), idxtranblhloc(:,ikloc), spinor_ud, &
     !$ACC            ngq(iq), ias2ic ) &
     !$ACC   COPYIN( ibatch )

#else
     ! Fill in bgntuju and b1 on CPU
     !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) &
     !$OMP   PRIVATE(ig,ias,ic,ki,iband,i,ist1,ic,l1)
#endif
     do ig=1,ngq(iq)
        do ias=1,natmtot

#ifndef _OPENACC
#ifdef _DEBUG_megqblh_
           ! OpenACC doesn't support WRITE statements in device code
           WRITE( dbgunit2, '(4(A,I4))' ) 'Batch ', ibatch, &
                                          ' ig=', ig, 'ias=', ias, &
                                          ' start=', k1, ' end=', k2, ' nsize=', nsize
#endif /* _DEBUG_megqblh_ */
#endif /* _OPENACC */

           !$OMP ATOMIC
           ibatch = ibatch + 1

           ic = ias2ic(ias)
           bgntuju(:,:,ibatch) = gntuju(:,:,ic,ig)

           b1(:,:,ibatch) = zzero
           b2(:,:,ibatch) = zzero

           ! Loop for a single batch
           DO ki = 1, nsize

              iband = k1 + ki - 1
              i = idxtranblhloc( iband, ikloc )
              ist1 = bmegqblh(1,i,ikloc)

              ! TODO: Dump bmegqblh and inspect, then get rid of the l1 check
              l1=.true.
              if (spinpol) then
                 if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
              endif
              if (l1) then
                 ! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
                 b1(:,ki,ibatch) = DCONJG( wfsvmt1(:,ias,ispn1,ist1) * &
                                           sfacgq(ig,ias) )
              END IF ! l1

           END DO ! ki

        enddo !ias
     enddo !ig
#ifdef _OPENACC
     !$ACC END PARALLEL LOOP
#else
     !$OMP END PARALLEL DO
#endif /* _OPENACC */

  END DO ! k1

#ifndef _OPENACC
  !CALL cudaMemcpy( d_bgntuju, bgntuju, cudaMemcpyHostToDevice )
  !CALL cudaMemcpy( d_b1, b1, cudaMemcpyHostToDevice )
#endif /* _OPENACC */

#ifdef _MAGMA_

  ! Perform batched ZGEMM on device using MAGMA
  CALL magma_zgemm_batched_f( 'N', 'N', nmt, nsize, nmt, &
                               zone,  bgntuju(:,:,:), nmt, &
                                      b1(:,:,:),      nmt, &
                               zzero, b2(:,:,:),      nmt, &
                               nbatch )

  ! Fetch result from device
#ifdef _OPENACC
  !$ACC UPDATE SELF( b2 )
#else
  !CALL cudaMemcpy( b2, d_b2, cudaMemcpyDeviceToHost )
#endif /* _OPENACC */

  ! Collect results into wftmp1mt
  ibatch = 0
  DO ki = 1, nblock

     k1 = (ki-1)*nb + 1
     IF( ki == nblock ) THEN
        k2 = idxhiband
     ELSE
        k2 = ki*nb
     END IF

     DO ig = 1, ngq(iq)
        DO ias = 1, natmtot

           ibatch = ibatch + 1
           wftmp1mt(:,k1:k2,ias,ig) = b2(:,:,ibatch)

        END DO ! ig
     END DO ! ias
  END DO ! ki

#else

  ibatch = 0
  ! Perform batched ZGEMM on CPU using OpenMP
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) &
  !$OMP   PRIVATE(ig,ias,ic,k1,k2,ki)
  DO ki = 1, nblock
     do ig=1,ngq(iq)
        do ias=1,natmtot
           ic = ias2ic(ias)

           ! Original code for historical purpose
           !do j=1,ngntuju(ic,ig)
           !  b2(igntuju(2,j,ic,ig))=b2(igntuju(2,j,ic,ig))+&
           !    &b1(igntuju(1,j,ic,ig))*gntuju(j,ic,ig)
           !enddo

           !$OMP ATOMIC
           ibatch = ibatch + 1

           ! Perform ZGEMM by batch
           ! Reminder: ZGEMM( transA, transB, M, N, K, alpha,
           !                  A, lda, B, ldb, beta, C, ldc )
           !           C := alpha*op(A)*op(B) + beta*C
           ! nmt = lmmaxapw*nufrmax (see line 91)
           ! nsize = number of bands per batch (maximum nb = 64, see line 30)
           CALL zgemm( 'N', 'N', nmt, nsize, nmt, &
                       zone,  gntuju(1,1,ic,ig), nmt, &
                              b1(1,1,ibatch), nmt, &
                       zzero, b2(1,1,ibatch), nmt )

           ! Note: this one is clearer but generates implicit copyin and copyout
           !CALL zgemm( 'N', 'N', nmt, nsize, nmt, &
           !             zone,  gntuju(1:nmt,1:nmt,ic,ig), nmt, &
           !                    b1(1:nmt,1:nsize),  nmt, &
           !             zzero, b2(1:nmt,1:nsize),  nmt )

           k1 = (ki-1)*nb + 1
           IF( ki == nblock ) THEN
              k2 = idxhiband
           ELSE
              k2 = ibatch*nb
           END IF

           ! Collect results into wftmp1mt
           wftmp1mt(:,k1:k2,ias,ig) = b2(:,:,ibatch)

        enddo !ias
     enddo !ig

  END DO ! ibatch
  !$OMP END PARALLEL DO
#endif

  call timer_stop(3)
  call papi_timer_stop(pt_megqblh_mt)

  ! Start the bounded do loop for each band
  DO iband = 1, idxhiband

! left <bra| state 
     wftmp1=zzero

     ! The starting point of the index "i" for accessing bmegqblh(:,i,:)
     ! for each iband and ikloc was stored as idxtranblhloc
     i = idxtranblhloc( iband, ikloc )
     ist1 = bmegqblh(1,i,ikloc)

     ! Same as above
     l1=.true.
     if (spinpol) then
        if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
     endif

     if (l1) then

     ! Fetch muffin-tin part from b2 and store it into wftmp1
     ! TODO: Split wftmp1 into muffin-tin and interstitial parts
     DO ig = 1, ngq(iq)
        DO ias = 1, natmtot
           wftmp1( (ias-1)*nmt+1:ias*nmt, ig ) = wftmp1mt( :, iband, ias, ig )
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
deallocate(wftmp1)
deallocate(wftmp2)
deallocate(wfir1)

!$ACC EXIT DATA DELETE( wfsvmt1 )

call papi_timer_stop(pt_megqblh)

return
end subroutine genmegqblh
