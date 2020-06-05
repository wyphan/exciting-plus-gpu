      subroutine readbin(fname, nrxyz, wf)

      implicit none

      integer, intent(in) :: nrxyz(3)
      character(256), intent(in) :: fname
      real(4), intent(out) :: wf(nrxyz(1),nrxyz(2),nrxyz(3))
      real(4), allocatable :: func(:)
      integer recl
      integer i,j,k,ir
      integer, parameter :: float_size = 4

      allocate(func(nrxyz(2)*nrxyz(3)))
      recl = float_size*nrxyz(2)*nrxyz(3)

      open(70, file=trim(fname), action='READ', &
           form='UNFORMATTED', access='DIRECT', recl=recl)

!...Read the data in exactly how it is written in wann_plot_3d.f90
      do i=1,nrxyz(1)
        read(70, rec=i)func(:)
        ir = 0
        do j=1,nrxyz(2)
          do k=1,nrxyz(3)
            ir = ir+1
            wf(i,j,k) = func(ir)
          enddo
        enddo
      enddo

      close(70)
      deallocate(func)

      end subroutine
