subroutine getmaxgnt(lmaxexp,maxgnt)
use modmain
implicit none
integer, intent(in) :: lmaxexp
integer, intent(out) :: maxgnt
integer l1,l2,l3,m1,m2,m3, gnt2
integer nrf1,nrf2
real(8), external :: gaunt
real(8) t1
! estimate the maximum number of Gaunt-like coefficients 
!maxgnt=0
!do l1=0,lmaxapw
!  nrf1=maxval(nufr(l1,:))
!  do l2=0,lmaxapw
!    nrf2=maxval(nufr(l2,:))
!    do m1=-l1,l1
!      do m2=-l2,l2
!        t1=0.d0
!        do l3=0,lmaxexp
!          do m3=-l3,l3
!            t1=t1+abs(gaunt(l2,l1,l3,m2,m1,m3))
!          enddo
!        enddo
!        if (t1.gt.1d-16) maxgnt=maxgnt+nrf1*nrf2
!      enddo
      !maxgnt=maxgnt+nrf1
    !enddo
!  enddo
!enddo

gnt2 = lmmaxapw + (nufrmax-1)*16

IF( gnt2 > 96 ) THEN
   maxgnt = 128
ELSE
   maxgnt = 96
END IF ! gnt2

!maxgnt=lmmaxapw*nufrmax

return
end
      
