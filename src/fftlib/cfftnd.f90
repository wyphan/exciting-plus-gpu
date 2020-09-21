subroutine cfftnd(nd,n,sgn,c)
!
! DESCRIPTION:
!  In-place fast Fourier transform for complex arrays in $n_d$ dimensions. The
!  forward transform is scaled by one over the size of the array. Uses a
!  modified version of the FFTPACK5 library.
!
! INPUT/OUTPUT PARAMETERS:
!   nd  : number of dimensions (in,integer)
!   n   : mesh size (in,integer(nd))
!   sgn : FFT direction, -1: forward, 1: backward (in,integer)
!   c   : array to transform (inout,complex(n(1)*n(2)*...*n(nd)))
!
!  Copyright (C) 2005 J. K. Dewhurst
!  Distributed under the terms of the GNU General Public License.
!  See the file COPYING for license details.
!
use fftpack5
implicit none
! arguments
integer, intent(in) :: nd
integer, intent(in) :: n(nd)
integer, intent(in) :: sgn
complex(prec), intent(inout) :: c(*)
! local variables
integer i,j,k,l,p,q,iw,iw1,ni
integer lensav,lenwrk
! allocatable arrays
real(prec), allocatable :: wsave(:)
real(prec), allocatable :: work(:)
if (nd.le.0) then
  write(*,*)
  write(*,'("Error(cfftnd): invalid number of dimensions : ",I8)') nd
  write(*,*)
  stop
end if
p=1
lensav=1
do i=1,nd
  if (n(i).le.0) then
    write(*,*)
    write(*,'("Error(cfftnd): invalid n : ",I8)') n(i)
    write(*,'(" for dimension ",I4)') i
    write(*,*)
    stop
  end if
  p=p*n(i)
  lensav=max(lensav,2*n(i)+int(log(real(n(i),prec)))+4)
end do
lenwrk=2*p
allocate(wsave(lensav))
allocate(work(lenwrk))
if (sgn.gt.0) then
  q=1
  do i=1,nd
    ni=n(i)
    if (ni.gt.1) then
      iw=ni+ni+1
      iw1=iw+1
      p=p/ni
      call cfftmi(ni,wsave,lensav)
      j=1
      k=q*ni
      do l=1,p
        call cmfm1b(q,1,ni,q,c(j),work,wsave,wsave(iw),wsave(iw1))
        j=j+k
      end do
      q=k
    end if
  end do
else
  q=1
  do i=1,nd
    ni=n(i)
    if (ni.gt.1) then
      iw=ni+ni+1
      iw1=iw+1
      p=p/ni
      call cfftmi(ni,wsave,lensav)
      j=1
      k=q*ni
      do l=1,p
        call cmfm1f(q,1,ni,q,c(j),work,wsave,wsave(iw),wsave(iw1))
        j=j+k
      end do
      q=k
    end if
  end do
end if
deallocate(wsave,work)
return
end subroutine

