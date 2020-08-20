subroutine getmeidx(allibt)
use modmain
use mod_nrkp
use mod_expigqr
implicit none

logical, intent(in) :: allibt
! local variables
integer i,ik,jk,ist1,ist2,ikloc,n,i1,i2,n1,n2
logical laddibt,ldocc,laddme,lwanibt
logical l11,l12,l21,l22,le1,le2
integer*2, allocatable :: wann_bnd_n(:,:)
integer*2, allocatable :: wann_bnd_k(:,:)
logical, external :: bndint

INTEGER :: idxloband1, idxhiband1, idxloband2, idxhiband2
INTEGER :: nband1, nband2, ntran1, ntran2, ntranchk
LOGICAL :: lwatch, l1stband1, l1stband2

! Reinitialize this value for every iq
! (flag is declared in module mod_expigqr line 56)
ltranconst = .TRUE.

if (wannier_megq) then
  allocate(wann_bnd_n(nstsv,nwantot))
  allocate(wann_bnd_k(nstsv,nkptnr))
  wann_bnd_n=0
  wann_bnd_k=0
! mark all bands that contribute to WF expansion
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do i=1,nstsv
      do n=1,nwantot
        if (abs(wanncnrloc(n,i,ikloc)).gt.1d-10) then
! at least at one k-point the state i contributes to WF n
          wann_bnd_n(i,n)=1
! for k-point ik band i contributes at least to one WF 
          wann_bnd_k(i,ik)=1
        endif
      enddo
    enddo
  enddo
  call mpi_grid_reduce(wann_bnd_n(1,1),nstsv*nwantot,dims=(/dim_k/),&
    &all=.true.,op=op_max)
  call mpi_grid_reduce(wann_bnd_k(1,1),nstsv*nkptnr,dims=(/dim_k/),&
    &all=.true.,op=op_max)
endif !wannier_megq

do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  jk=idxkq(1,ik)
  i=0

  ! Reinitialize these values for every ikloc
  idxloband1 = 1
  idxloband2 = 1
  idxhiband1 = -1
  idxhiband2 = -1
  nband1 = 0
  nband2 = 0
  ntran1 = 0
  ntran2 = 0

  do ist1=1,nstsv

    ! Start watching for first "hit" on laddme
    lwatch = .TRUE.

    ! Reinitialize this value for every ikloc and ist1
    ntranchk = 0
    l1stband1 = .FALSE.
    l1stband2 = .FALSE.

    do ist2=1,nstsv
      lwanibt=.false.
! include transition between bands ist1 and ist2 when:
!  1) both bands ist1 and ist2 fall into energy interval and
!     difference of band occupation numbers in not zero or all
!     interband transitions are required  
      le1=bndint(ist1,evalsvnr(ist1,ik),megq_include_bands(1),&
        megq_include_bands(2))
      le2=bndint(ist2,evalsvnr(ist2,jk),megq_include_bands(1),&
        megq_include_bands(2))
!      ldocc=abs(occsvnr(ist1,ik)-occsvnr(ist2,jk)).gt.1d-6
      ldocc=.true.
      laddibt=(le1.and.le2).or.allibt
!      laddibt=(le1.and.le2.and.ldocc).or.allibt
!  2) this bands are necessary to compute matrix elements in Wannier basis
      if (wannier_megq) then
! if this bands contribute to at least one Wannier function
        if (wann_bnd_k(ist1,ik).ne.0.and.wann_bnd_k(ist2,jk).ne.0) then
! check contribution to each Wannier function
! NOTE: current implementation works with small overhead; in principle we
!  must check Bloch contribution to particular WF at particular k-point, but
!  this requires much bigger array of "Bloch contribution to WF" flags
          do i1=1,megqwantran%nwan
            n1=megqwantran%iwan(i1)
            do i2=1,megqwantran%nwan
              n2=megqwantran%iwan(i2)
              if (megqwantran%wt(n1,n2).eq.1.and.wann_bnd_n(ist1,n1).eq.1.and.&
                  &wann_bnd_n(ist2,n2).eq.1.and.ldocc) then
                laddibt=.true.
                lwanibt=.true.
              endif
            enddo !i2
          enddo !i1
        endif
      endif !wannier_megq
      laddme=.false.

      ! Reinitialize spin pair checks
      l11 = .FALSE.
      l12 = .FALSE.
      l21 = .FALSE.
      l22 = .FALSE.

! final check: don't add matrix element if it is zero
      if (laddibt) then
        if (.not.spinpol) then
          laddme=.true.
          l11 = .TRUE.
        else
          l11=spinor_ud(1,ist1,ik).eq.1.and.spinor_ud(1,ist2,jk).eq.1
          l12=spinor_ud(1,ist1,ik).eq.1.and.spinor_ud(2,ist2,jk).eq.1
          l21=spinor_ud(2,ist1,ik).eq.1.and.spinor_ud(1,ist2,jk).eq.1
          l22=spinor_ud(2,ist1,ik).eq.1.and.spinor_ud(2,ist2,jk).eq.1
          if (expigqr22.eq.1.and.(l11.or.l22)) THEN
             laddme=.true.
          ELSE if (expigqr22.eq.2.and.(l12.or.l21)) THEN
             laddme=.true.
          END IF
        endif
      endif
      if (laddme) then
        i=i+1

        ! When first "hit" on laddme happens for every ist1 and ikloc,
        ! the "watcher" is active ( lwatch .EQV. .TRUE. )
        IF( lwatch ) THEN

           ! Record position of first "hit" for this ist and ikloc
           ! (array is declared in module mod_expigqr line 67
           !       and allocated in init_band_trans() line 31)
           idxtranblhloc(ist1,ikloc) = i
           
           IF( expigqr22 == 1 ) THEN
              IF( l11 ) nband1 = nband1 + 1
              IF( l22 ) nband2 = nband2 + 1
           ELSE IF( expigqr22 == 2 ) THEN
              IF( l12 ) nband1 = nband1 + 1
              IF( l21 ) nband2 = nband2 + 1
           END IF ! expigqr22

           ! Stop watching
           lwatch = .FALSE.

        END IF ! lwatch
        
        IF( expigqr22 == 1 ) THEN
           IF( l11 ) THEN
              ! For the first "hit" at the first band, idxhiband was still -1
              IF( idxhiband1 == -1 ) THEN
                 l1stband1 = .TRUE.
                 idxloband1 = ist1
              END IF
              ! Update the corresponding idxhiband
              ! up+up
              idxhiband1 = ist1
           END IF
           IF( l22 ) THEN
              IF( idxhiband2 == -1 ) THEN
                 l1stband2 = .TRUE.
                 idxloband2 = ist1
              END IF
              ! dn+dn
              idxhiband2 = ist1
           END IF
        ELSE IF( expigqr22 == 2 ) THEN
           IF( l12 ) THEN
              IF( idxhiband1 == -1 ) THEN
                 l1stband1 = .TRUE.
                 idxloband1 = ist1
              END IF
              ! up+dn
              idxhiband1 = ist1
           END IF
           IF( l21 ) THEN
              IF( idxhiband2 == -1 ) THEN
                 l1stband2 = .TRUE.
                 idxloband2 = ist1
              END IF
              ! dn+up
              idxhiband2 = ist1
           END IF
        END IF ! expigqr22

        ! Accumulate number of paired bands
        ! (like i, but only for the first band of each kind)
        IF( l1stband1 .AND. .NOT.l1stband2 ) ntran1 = ntran1 + 1
        IF( .NOT.l1stband1 .AND. l1stband2 ) ntran2 = ntran2 + 1 

        ! Counter to make sure ntran stays the same across bands
        ! (accumulate number of paired bands for each ist1)
        ntranchk = ntranchk + 1

        bmegqblh(1,i,ikloc)=ist1
        bmegqblh(2,i,ikloc)=ist2
        if (lwanibt) then
          nmegqblhwan(ikloc)=nmegqblhwan(ikloc)+1
          imegqblhwan(nmegqblhwan(ikloc),ikloc)=i
        endif

      endif !laddme

    enddo !ist2

    ! Make sure ntran stays the same across different bands
    ltranconst = ( ntranchk == 0 ) .OR. ( ntranchk == ntran1 ) .OR. ( ntranchk == ntran2 )
    IF( .NOT. ltranconst ) THEN
       WRITE(*, '( "Warning[getmeidx]: rank ", I5, ": ntran is different ",&
                   &"for ist1=", I6, " ikloc=", I6 )' ) iproc, ist1, ikloc
    END IF ! ltransconst

  enddo !ist1
  nmegqblh(ikloc)=i

  ! Save number of |n',k+q> Bloch kets paired to each <n=ist1,k| bra
  ntranblhloc(1,ikloc) = ntran1 ! up+up, or when expigqr22=2, up+dn 

  ! Store idxloband, that is, the 'lowest' band index with transitions
  idxlobandblhloc(1,ikloc) = idxloband1

  ! Store idxhiband, that is, the 'highest' band index with transitions
  idxhibandblhloc(1,ikloc) = idxhiband1

  ! Store nband, the number of paired bras
  nbandblhloc(1,ikloc) = nband1

  ! dn+dn or when expigqr22=2, dn+up
  IF( spinpol ) THEN
     ntranblhloc(2,ikloc) = ntran2
     idxlobandblhloc(2,ikloc) = idxloband2
     idxhibandblhloc(2,ikloc) = idxhiband2
     nbandblhloc(2,ikloc) = nband2
  END IF

#if EBUG > 1
  WRITE(*,*) 'getmeidx: ikloc=', ikloc, ' idxlobandblhloc=', idxlobandblhloc(:,ikloc)
  WRITE(*,*) 'getmeidx: ikloc=', ikloc, ' idxhibandblhloc=', idxhibandblhloc(:,ikloc)
  WRITE(*,*) 'getmeidx: ikloc=', ikloc, ' nbandblhloc=', nbandblhloc(:,ikloc)
  WRITE(*,*) 'getmeidx: ikloc=', ikloc, ' ntranblhloc=', ntranblhloc(:,ikloc)
#endif

enddo !ikloc
if (wannier_megq) then
  deallocate(wann_bnd_n)
  deallocate(wann_bnd_k)
endif
return
end
