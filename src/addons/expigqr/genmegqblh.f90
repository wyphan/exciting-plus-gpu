subroutine genmegqblh(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
  &wfsvit1,wfsvit2)
use modmain
use mod_addons_q
use mod_nrkp
use mod_expigqr
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
integer i,j,ik,jk,igkq,n1,ispn1,ispn2,ist1,ist2,ic
integer ig,ig1,ig2,ias,ifg,ir
logical l1
complex(8), allocatable :: wftmp1(:,:)
complex(8), allocatable :: wftmp2(:,:)
complex(8), allocatable :: wfir1(:)
complex(8) b1(lmmaxapw*nufrmax),b2(lmmaxapw*nufrmax) ! TODO: convert to matrices

#ifdef _DEBUG_bmegqblh_
  INTEGER :: dbgcnt, dbgunit
#endif // _DEBUG_bmegqblh_

!--begin Convert do while into bounded do loop

  INTEGER :: idxhiband, iband, ntran, idxtran
  EXTERNAL :: zcopy

!--end Convert do while into bounded do loop

wfsize=lmmaxapw*nufrmax*natmtot+ngknr2
allocate(wftmp1(wfsize,ngq(iq))) ! TODO: Change dimensions appropriately
allocate(wftmp2(wfsize,nstsv))   ! TODO: Check contiguity of ZCOPY transfers
allocate(wfir1(ngrtot))
call papi_timer_start(pt_megqblh)

! global k-point
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! jk=k+q-G_q
jk=idxkq(1,ik)
! G_q vector 
igkq=idxkq(2,ik)

#ifdef _DEBUG_bmegqblh_
  dbgunit = 1000+iproc ! Make sure this matches the definition in mod_expigqr
  dbgcnt=1

  WRITE( dbgunit, '(A,I3,A,I5)' ) 'nmegqblh(ikloc=', ikloc, ') = ', nmegqblh(ikloc)
#endif // _DEBUG_bmegqblh_

do ispn1=1,nspinor
  if (expigqr22.eq.1) ispn2=ispn1

!--begin Convert do while into bounded do loop

! index of the interband transitions
     !i=1
! go through the interband transitions    
     !do while (i.le.nmegqblh(ikloc))

     ! This do while loop converts easily to a bounded do loop

     ! Convert to the corresponding ist1 loop in getmeidx() line 56-149
     ! as stored into idxhibandblh(ikloc=1:nkptnr) at getmeidx() line 155,
     ! skipping as necessary (with a warning message... it should NOT happen!)
     idxhiband = idxhibandblhloc(ikloc)
     IF( idxhiband == 0 ) THEN
        ! Unit 151 is either 'CRPA.OUT' or 'RESPONSE.OUT'
        WRITE(151, '( "Warning[genmegqblh]: highest band is zero for iq=",&
                    &I6, " ikloc=", I6, " ispn1=", I1 )' ) iq, ikloc, ispn1
        CYCLE
     END IF

     ! Otherwise, start the bounded do loop
     DO iband = 1, idxhiband

        ! The starting point of the index "i" for accessing bmegqblh(:,i,:)
        ! for each iband and ikloc was stored as idxtranblhloc
        i = idxtranblhloc( iband, ikloc )

!--end Convert do while into bounded do loop

! left <bra| state 
    ist1=bmegqblh(1,i,ikloc)
    wftmp1=zzero
    l1=.true.
    if (spinpol) then
      if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
    endif
    if (l1) then
      call timer_start(3)
      call papi_timer_start(pt_megqblh_mt)

      ! TODO: insert OMP directives here
      do ig=1,ngq(iq)
! precompute muffin-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
        do ias=1,natmtot
          b1=dconjg(wfsvmt1(:,ias,ispn1,ist1)*sfacgq(ig,ias))
          ic=ias2ic(ias)
          b2=zzero
          !do j=1,ngntuju(ic,ig)
          !  b2(igntuju(2,j,ic,ig))=b2(igntuju(2,j,ic,ig))+&
          !    &b1(igntuju(1,j,ic,ig))*gntuju(j,ic,ig)
          !enddo

          ! TODO: convert to true ZGEMM
          call zgemm('N','N',lmmaxapw*nufrmax,1,lmmaxapw*nufrmax,&
            &zone,gntuju(1,1,ic,ig),lmmaxapw*nufrmax,b1,lmmaxapw*nufrmax,&
            &zzero,b2,lmmaxapw*nufrmax)
          wftmp1((ias-1)*lmmaxapw*nufrmax+1:ias*lmmaxapw*nufrmax,ig)=b2(:)

        enddo !ias
      enddo !ig  
      call timer_stop(3)
      call papi_timer_stop(pt_megqblh_mt)
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

!--begin Convert do while into bounded do loop

! collect right |ket> states into matrix wftmp2

        !n1=0
        !do while ((i+n1).le.nmegqblh(ikloc))
           !if (bmegqblh(1,i+n1,ikloc).ne.bmegqblh(1,i,ikloc)) exit
           !ist2=bmegqblh(2,i+n1,ikloc)
           !n1=n1+1
           !call memcopy(wfsvmt2(1,1,1,ispn2,ist2),wftmp2(1,n1),16*lmmaxapw*nufrmax*natmtot)
           !call memcopy(wfsvit2(1,ispn2,ist2),wftmp2(lmmaxapw*nufrmax*natmtot+1,n1),16*ngknr2)
        !enddo ! while (i+n1) <= nmegqblh(ikloc)

        ! This do while loop is trickier to convert

        ! For this <ist1=n| bra, load number of matching |ist2=n'> ket states,
        ! and this will be the upper bound for the loop
        ! Note: ntran can be zero (zero-trip loop)
        ! TODO: optimize search algorithm in getmeidx()
        ntran = ntranblhloc(iband,ikloc)

#ifdef _DEBUG_bmegqblh_
IF( ntran > 0 ) THEN
  WRITE( dbgunit, '(7(1X,I5))' ) dbgcnt, ikloc, iq, iband, i, ntran, i+ntran-1
  dbgcnt = dbgcnt + 1
END IF
#endif // _DEBUG_bmegqblh_

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

        !call zgemm('T','N',n1,ngq(iq),wfsize,zone,wftmp2,wfsize,wftmp1,wfsize,&
        !     &zone,megqblh(i,1,ikloc),nstsv*nstsv)
        !i=i+n1
        !call timer_stop(5)
     !enddo ! while i <= nmegqblh(ikloc)

        ! This particular ZGEMM() call corresponds with line 9 of Algorithm 2
        ! in the Gordon Bell paper
        CALL zgemm( 'T', 'N', ntran, ngq(iq), wfsize, zone, &
                    wftmp2, wfsize, wftmp1, wfsize, zone, &
                    megqblh(i,1,ikloc), nstsv*nstsv )

        ! We don't need to add n1 to i anymore to move on to the next <nk| bra
        ! since it is already stored as 
        ! For reference, n1 is now ntran = ntrangqblhloc(iband,ikloc)

        CALL timer_stop(5) ! Same as before

     END DO ! iband; replaces do while loop i <= nmegqblh(ikloc)

#ifdef _DEBUG_bmegqblh_
     WRITE( dbgunit, '(A,I3)') 'highest band = ', idxhiband
#endif // _DEBUG_bmegqblh_

!--end Convert do while into bounded do loop

enddo !ispn
deallocate(wftmp1)
deallocate(wftmp2)
deallocate(wfir1)

call papi_timer_stop(pt_megqblh)

return
end subroutine genmegqblh
