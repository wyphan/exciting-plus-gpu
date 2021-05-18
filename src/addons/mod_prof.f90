MODULE mod_prof
!===============================================================================
! Implements a basic profiling system using the system_clock() intrinsic
! Last edit: Apr 29, 2021 (ED)
!===============================================================================

  USE mod_timer, ONLY: dclock

  implicit none
  private

  integer, parameter :: maxroutines = 1024
  integer, parameter :: maxlevels = 1024
  integer, parameter :: maxstrlen = 127

  real, dimension(maxroutines) :: dictstart, dicttotal
  integer, dimension(maxroutines) :: dictcount, dictlevel
  integer :: nroutine, nlevels

  character(len=maxstrlen), dimension(maxroutines) :: dictname
  character(len=maxstrlen), dimension(maxlevels) :: lastroutine

  public :: profstart, profend, profstat, profinit

CONTAINS

!===============================================================================

  subroutine profinit()
    IMPLICIT NONE

    nroutine = 0

    dictname(:) = ' '
    dictstart(:) = 0.0
    dictcount(:) = 0.0
    dicttotal(:) = 0.0

    nlevels = 0
    lastroutine(:) = ' '

  end subroutine profinit

!===============================================================================

  subroutine profstart(rname)
    IMPLICIT NONE

    character(len=*),intent(in) :: rname

    character(len=maxstrlen) :: name
    logical :: found,isok
    integer :: i,j,ipos

    name = rname
    nlevels = nlevels + 1

    isok = (1 .le. nlevels).and.(nlevels .le. maxlevels)
    call assert( isok, &
                 '** profstart: invalid nlevels ', nlevels )
        
    lastroutine(nlevels) = name
    found = .false.
    do j=1,nroutine
       i = nroutine - j + 1
       if (dictname(i)(1:1).eq.name(1:1)) then
          found = (dictname(i) .eq. name)
          if (found) then
             ipos = i
             exit
          endif
       endif
    enddo

    if (.not.found) then
       nroutine = nroutine + 1
       isok = (nroutine .le. maxroutines)
       call assert( isok, &
                    '** profstart: nroutine > maxroutines ', nroutine )

       ipos = nroutine
       dictname(ipos) = name
       dictcount(ipos) = 0
       dicttotal(ipos) = 0.0
       dictlevel(ipos) = nlevels
    endif

    dictstart(ipos) = dclock()
    dictcount(ipos) = dictcount(ipos) + 1

    return
  end subroutine profstart

!===============================================================================

  subroutine profend(rname)
    IMPLICIT NONE

    character(len=*),intent(in) :: rname

    character(len=maxstrlen) :: name
    integer :: i,j,ipos
    logical :: found,isok
    real :: tend

    name = rname
    tend = dclock()

    isok = (1.le.nlevels).and.(nlevels.le.maxlevels)
    call assert(isok, &
                '** profend: invalid nlevels ', nlevels )

    isok = (name .eq. lastroutine(nlevels))
    if (.not.isok) then
       print*,'** profend name != lastroutine(',nlevels,') '
       print*,'name: ', name
       print*,'lastroutine(nlevels): ', lastroutine(nlevels)

       stop '** error ** '
    endif

    found = .false.
    do j=1,nroutine
       i = nroutine - j + 1

       if (dictname(i)(1:1) .eq. name(1:1)) then
          found = (dictname(i) .eq. name)
          if (found) then
             ipos = i
             exit
          endif
       endif
    enddo

    if (.not.found) then
       print*,'** profend: routine name not found '
       print*,'name: ',name
       stop '** error ** '
    endif

    dicttotal(ipos) = dicttotal(ipos) + (tend - dictstart(ipos));
    nlevels = nlevels - 1;

    return
  end subroutine profend

!===============================================================================

  subroutine profstat(outdev_in)
    implicit none
    integer, optional, intent(in):: outdev_in 

    character(len=maxstrlen) :: fname,fstr
    character(len=2) :: x
    integer :: i, outdev

    if (present(outdev_in)) then
       outdev = outdev_in
    else
       outdev = 16
    endif

    fname = 'profstat.dat'
    open(outdev, file=fname, form='formatted', &
         access='sequential',status='unknown')
    rewind(outdev)

    do i=1,nroutine
       write(x,'(I2)') dictlevel(i)
       fstr = "("//x//"X,A24,' was called ',i10,' times, total ',f10.2,' secs')"
       write(outdev,fstr) dictname(i), dictcount(i), dicttotal(i)
       write(*,fstr) dictname(i), dictcount(i), dicttotal(i)
    enddo

    close(outdev)
    return
  end subroutine profstat

!===============================================================================

END MODULE mod_prof
