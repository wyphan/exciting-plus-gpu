      subroutine readr3d(funit, ndim, nrxyz, coords, rho)

      implicit none

      integer, intent(in) :: funit, ndim
      integer, intent(in) :: nrxyz(3)
      real(4), intent(out) :: coords(3,ndim)
      real(4), intent(out) :: rho(nrxyz(1),nrxyz(2),nrxyz(3))
      integer, allocatable :: ip2x(:),ip2y(:),ip2z(:)
      integer i,ip,ip1,ip2,ip3,np
      real(4) v1(3),r1
      
      ip=0
      do ip3=0,nrxyz(3)-1
        do ip2=0,nrxyz(2)-1
          do ip1=0,nrxyz(1)-1
            ip=ip+1
          end do
        end do
      end do
      np=ip
      
      allocate(ip2x(np),ip2y(np),ip2z(np))
      ip=0
      do ip3=0,nrxyz(3)-1
        do ip2=0,nrxyz(2)-1
          do ip1=0,nrxyz(1)-1
            ip=ip+1
            ip2x(ip)=ip1
            ip2y(ip)=ip2
            ip2z(ip)=ip3
          end do
        end do
      end do
      
      do ip=1,np
        read(80,'(4G18.10)')v1(:),rho(ip2x(ip)+1,ip2y(ip)+1,ip2z(ip)+1)
        if (ip2x(ip).eq.0.and.ip2y(ip).eq.0) coords(3,ip2z(ip)+1)=v1(3)
        if (ip2x(ip).eq.0.and.ip2z(ip).eq.0) coords(2,ip2y(ip)+1)=v1(2)
        if (ip2y(ip).eq.0.and.ip2z(ip).eq.0) coords(1,ip2x(ip)+1)=v1(1)
      end do
      close(80)

      end subroutine
