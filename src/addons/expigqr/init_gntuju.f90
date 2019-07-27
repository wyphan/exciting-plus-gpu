subroutine init_gntuju(iq,lmaxexp)
use modmain
use mod_wannier
use mod_expigqr
use mod_addons_q
implicit none
integer, intent(in) :: iq
integer, intent(in) :: lmaxexp
call getmaxgnt(lmaxexp,ngntujumax)
if (allocated(ngntuju)) deallocate(ngntuju)
allocate(ngntuju(natmcls,ngq(iq)))
ngntuju=0
if (allocated(igntuju)) deallocate(igntuju)
allocate(igntuju(2,ngntujumax,natmcls,ngq(iq)))
igntuju=0
if (allocated(gntuju)) deallocate(gntuju)
!allocate(gntuju(ngntujumax,ngntujumax,natmcls,ngq(iq)))

!-- TODO: Remove gntujutmp
if( allocated(gntujutmp) ) deallocate(gntujutmp)
allocate(gntujutmp( ngntujumax, ngntujumax, natmcls, ngq(iq) ))
allocate(gntuju( ngntujumax, ngntujumax*ngq(iq), natmcls ))
!--

gntuju=zzero
call gengntuju(iq,lmaxexp)
return
end
