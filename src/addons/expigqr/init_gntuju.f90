subroutine init_gntuju(iq,lmaxexp)
use modmain
use mod_wannier
use mod_expigqr
USE mod_sparse, ONLY: nmt
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
allocate(gntuju(ngntujumax,ngntujumax,natmcls,ngq(iq)))
gntuju=zzero

#ifdef _PACK_gntuju_

if (allocated(nmt)) deallocate(nmt)
allocate( nmt(natmcls,ngq(iq)) )
nmt=0

#endif /* _PACK_gntuju_ */

call gengntuju(iq,lmaxexp)
return
end
