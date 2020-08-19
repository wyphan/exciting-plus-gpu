subroutine genmegqblh(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
                      wfsvit1,wfsvit2)
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
  USE mod_expigqr, ONLY: expigqr22, gntuju, megqblh, bmegqblh, nmegqblh, idxkq, &
                         idxhibandblhloc, ntranblhloc, idxtranblhloc, ltranconst
  USE mod_genmegqblh_gpu

!--DEBUG
  USE mod_lapack, ONLY: ZGEMM
!--DEBUG

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

  ! Temporary array to hold results for muffin-tin calculation
  ! (will be removed after everything is ported to GPU)
  COMPLEX(KIND=dz), DIMENSION(:,:,:,:), ALLOCATABLE :: wftmp1mt
 
#if defined(_DEBUG_bmegqblh_) || defined(_DEBUG_megqblh_)
  INTEGER :: dbgcnt1, dbgcnt2, dbgunit1, dbgunit2
#endif /* _DEBUG_bmegqblh_ || _DEBUG_megqblh_ */

  INTEGER :: idxhiband, iband, ntran, idxtran, ispst
  EXTERNAL :: zcopy

!--DEBUG
  INTEGER :: ibatch, iblock
  COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: mybgntuju, myb1, myb2
!--DEBUG

wfsize=lmmaxapw*nufrmax*natmtot+ngknr2
allocate(wftmp1(wfsize,ngq(iq))) ! TODO: Change dimensions appropriately
allocate(wftmp2(wfsize,nstsv))   ! TODO: Check contiguity of ZCOPY transfers
allocate(wfir1(ngrtot))
call papi_timer_start(pt_megqblh)

  ! Note: List of OpenACC variables that are already in device memory 
  !       due to inheritance from mod_expigqr::genmegq() :
  !         sfacgq, gntuju, bmegqblh, idxhibandblhloc, idxtranblhloc,
  !         spinor_ud, ngq, ias2ic

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
#endif

  ! Number of muffin-tin elements
  nmt = lmmaxapw * nufrmax

  ! Number of G+q vectors for a particular value of q-vector
  ngqiq = ngq(iq)

!--DEBUG
#if EBUG > 0
  WRITE(*,*) 'genmegqblh: iq=', iq, ' ikloc=', ikloc, ' ngq(iq)=', ngq(iq)
#endif
!--DEBUG

  ! Number of blocks and batches, blocked version
  !idxhiband = idxhibandblhloc(ikloc)
  !nblock = CEILING( REAL(idxhiband)/REAL(nb) )
  !nbatch = ngqiq * natmtot * nblock
  
  ! Number of blocks and batches, unblocked version
  nblock = 1
  nbatch = ngqiq * natmtot

  do ispn1=1,nspinor

     ! expigqr22 is always 1, for now (see mod_expigqr)
     if (expigqr22.eq.1) ispn2=ispn1

     ! Convert to the corresponding ist1 loop in getmeidx() line 56-149
     ! as stored into idxhibandblh(ikloc=1:nkptnr) at getmeidx() line 155,
     ! skipping as necessary (with a warning message... it should NOT happen!)
     nband1 = idxhibandblhloc(ispn1,ikloc)
     IF( nband1 == 0 ) THEN
        WRITE(*, '( "Warning[genmegqblh]: highest band is zero for ispn1 "&
                   &" iq=", I6, " ikloc=", I6 )' ) ispn1, iq, ikloc
        RETURN
     END IF

!--begin Convert to true ZGEMM

     call timer_start(3)
     call papi_timer_start(pt_megqblh_mt)

     ! Note that the loop order has been switched
     ! such that iband loop is now the innermost loop

     ! Allocate arrays on CPU memory
     ALLOCATE( wftmp1mt( nmt, nband1, natmtot, ngqiq ))
     ALLOCATE( spinstidx(nstsv) )
     !ALLOCATE( b1( nmt, nb, nbatch ))      ! Blocked version
     !ALLOCATE( b2( nmt, nb, nbatch ))      ! Blocked version
     ALLOCATE( b1( nmt, nband1, nbatch )) ! Unblocked version
     ALLOCATE( b2( nmt, nband1, nbatch )) ! Unblocked version
     ALLOCATE( batchidx( natmtot, ngqiq, nblock ))
#ifdef _OPENACC
     ALLOCATE( dptr_gntuju( nbatch ))
     ALLOCATE( dptr_b1( nbatch ))
     ALLOCATE( dptr_b2( nbatch ))
#else
     ALLOCATE( bgntuju( nmt, nmt, nbatch ))
#endif /* _OPENACC */

     ! Allocate arrays on GPU memory
     !$ACC DATA COPYIN( natmtot, ngqiq, nblock, nmt, nbatch, nstsv, nband1 ) &
     !$ACC      CREATE( spinstidx, nstspin, lcontig, &
     !$ACC              wftmp1mt, b1, b2, batchidx, &
     !$ACC              dptr_gntuju, dptr_b1, dptr_b2 )

     ! Count spin up states for this particular k-vector (replaces l1 check)
     ! Note: spinup and spindn are defined in mod_genmegqblh_gpu
     IF( ispn1 == 1 ) THEN
        ! Spin up
        CALL genmegqblh_countspin( spinup, ikloc )
     ELSE
        ! Spin down (never executed if spinpol = .FALSE. )
        CALL genmegqblh_countspin( spindn, ikloc )
     END IF

!------------------------------------------------------------------------------
! Kernel 1: Fill in bgntuju and b1, and zero b2
!------------------------------------------------------------------------------

!--DEBUG
#if EBUG > 1
     WRITE(*,*) 'genmegqblh: before 1st kernel'
#endif
!--DEBUG

     CALL genmegqblh_fillbatch( wfsvmt1, ikloc, ispn1 )

!--DEBUG
#if EBUG > 1
     WRITE(*,*) 'genmegqblh: after 1st kernel'
#endif
!--DEBUG
     
!------------------------------------------------------------------------------
! Kernel 2: Perform batched ZGEMM b2(:,:) = b1(:,:) x bgntuju(:,:)
!------------------------------------------------------------------------------

     ! Original code (retained for historical purpose)
     !do j=1,ngntuju(ic,ig)
     !  b2(igntuju(2,j,ic,ig))=b2(igntuju(2,j,ic,ig))+&
     !    &b1(igntuju(1,j,ic,ig))*gntuju(j,ic,ig)
     !enddo

!--DEBUG
#if EBUG > 1
     WRITE(*,*) 'genmegqblh: before 2nd kernel'
#endif
!--DEBUG

     CALL genmegqblh_batchzgemm()

!--DEBUG
#if EBUG > 1     
     WRITE(*,*) 'genmegqblh: after 2nd kernel'
#endif
!--DEBUG

!------------------------------------------------------------------------------
! Kernel 3: Save results to wftmp1mt and transfer back to CPU (for now)
!------------------------------------------------------------------------------

!--DEBUG
#if EBUG > 1
     WRITE(*,*) 'genmegqblh: before 3rd kernel'
#endif
!--DEBUG

     CALL genmegqblh_fillresult( wftmp1mt )

!--DEBUG
#if EBUG > 1     
     WRITE(*,*) 'genmegqblh: after 3rd kernel'
#endif
!--DEBUG

!------------------------------------------------------------------------------

     ! Transfer data to CPU (for now)
     !$ACC UPDATE SELF( wftmp1mt )

  call timer_stop(3)
  call papi_timer_stop(pt_megqblh_mt)

  ! Start the bounded do loop for each band
  DO ispst = 1, nband1

     ! left <bra| state
     wftmp1=zzero

     ! The starting point of the index "i" for accessing bmegqblh(:,i,:)
     ! for each iband and ikloc was stored as idxtranblhloc
     ! Note that spinstidx stores the band indices for a single spin projection
     iband = spinstidx( ispst )
     i = idxtranblhloc( iband, ikloc )
     ist1 = bmegqblh(1,i,ikloc)

     ! Note: wftmp1 combines the muffin-tin and interstitial parts for each band,
     !         to prepare for the second ZGEMM below
     !       Complete removal of wftmp1mt is impossible until
     !         interstitial part also ported to GPU (cuFFT with fallback to FFTW)
     DO ig = 1, ngqiq
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
      do ig=1,ngqiq
        do ig2=1,ngknr2
! G1=G2-G-Gkq
          ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,igqig(ig,iq))-ivg(:,igkq)
          ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
          wftmp1(lmmaxapw*nufrmax*natmtot+ig2,ig)=dconjg(wfir1(ifg))
        enddo
      enddo
      call timer_stop(4)      
      call papi_timer_stop(pt_megqblh_it)

    call timer_start(5)

#ifdef _DEBUG_bmegqblh_
IF( ntran > 0 ) THEN
  WRITE( dbgunit1, '(7(1X,I5))' ) dbgcnt1, ikloc, iq, iband, i, ntran, i+ntran-1
  dbgcnt1 = dbgcnt1 + 1
END IF
#endif /* _DEBUG_bmegqblh_ */

  ! Load number of matching |ist2=n'> ket states for each <ist1=n| bra
  IF( ltranconst ) THEN
     ntran = ntranblhloc(ispn1,ikloc)
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
    CALL zgemm( 'T', 'N', ntran, ngqiq, wfsize, zone, &
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

     ! natmtot, ngqiq, nblock, nbatch, nmt, nstsv, nband1,
     ! wf1tmpmt, spinstidx, nstspin, lcontig,
     ! b1, b2, batchidx, dptr_gntuju, dptr_b1, dptr_b2
     !$ACC END DATA     

     ! Clean up
     DEALLOCATE( spinstidx )
     DEALLOCATE( wftmp1mt )
     DEALLOCATE( b1 )
     DEALLOCATE( b2 )
     DEALLOCATE( batchidx )
#ifdef _OPENACC
     DEALLOCATE( dptr_gntuju )
     DEALLOCATE( dptr_b1 )
     DEALLOCATE( dptr_b2 )
#else
     DEALLOCATE( bgntuju )
#endif /* _OPENACC */

  enddo !ispn

deallocate(wftmp1)
deallocate(wftmp2)
deallocate(wfir1)

call papi_timer_stop(pt_megqblh)

return
end subroutine genmegqblh
