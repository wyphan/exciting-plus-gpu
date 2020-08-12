      program main

      implicit none

      integer lmax, i, l, ld, lm1, lm2
      real(8) ang(3), rot(3,3)
      real(8),  allocatable :: d(:,:)
      character(20) fmt
      
      open(10, file='genlps.in', action='READ', form='FORMATTED')
      read(10,*)lmax
      do i=1,3
        read(10,*)rot(:,i)
      enddo
      close(10)

      ld = 1
      do l=1,lmax
        ld = ld+2*l*l+1
      enddo
      allocate(d(ld,ld))

      call euler(rot, ang)
      call rlmrot(1, ang(1), ang(2), ang(3), lmax, ld, d)

      open(10, file='LPS.OUT', action='WRITE', form='FORMATTED')
      lm1=1
      do l=1,lmax
        write(10,'(I3)')l
        lm1 = lm1+2*(l-1)+1
        lm2 = lm1+2*l
        write(fmt,'("(",I0,"F14.10)")')2*l+1
        do i=lm1,lm2
          write(10,fmt)d(i,lm1:lm2)
        enddo
      enddo

      end program
