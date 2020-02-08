      subroutine readdxh(fname, nrxyz, orig, delta)

      implicit none

      character(256), intent(in) :: fname
      integer, intent(out) :: nrxyz(3)
      real(4), intent(out) :: orig(3), delta(3,3)
      integer i

      open(80, file=trim(fname), action='READ', form='FORMATTED')
      read(80, '(35x,3I4)') &
          nrxyz(1),nrxyz(2),nrxyz(3)
      read(80,'(7x,3G18.10)')orig(:)
      do i=1,3
        read(80,'(6x,3G18.10)')delta(:,i)
      enddo
      close(80)

      end subroutine
