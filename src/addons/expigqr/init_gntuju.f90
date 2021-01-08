subroutine init_gntuju(iq,lmaxexp)

use modmain
use mod_wannier
use mod_expigqr

#ifdef _PACK_gntuju_
  USE mod_genmegqblh_gpu, ONLY: nmt
#endif /*_PACK_gntuju_ */

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

#ifdef _PACK_gntuju_

  if (allocated(nmt)) deallocate(nmt)
  allocate( nmt(natmcls,ngq(iq)) )
  nmt(:,:) = 0

  ALLOCATE( irownz(ngntujumax,natmcls,ngq(iq)) )
  ALLOCATE( icolnz(ngntujumax,natmcls,ngq(iq)) )
  ALLOCATE( irows(2,ngntujumax,natmcls,ngq(iq)) )
  ALLOCATE( lfit(natmcls,ngq(iq)) )

#endif /* _PACK_gntuju_ */

  if (allocated(gntuju)) deallocate(gntuju)
  allocate(gntuju(ngntujumax,ngntujumax,natmcls,ngq(iq)))
  gntuju=zzero

  call gengntuju(iq,lmaxexp)

return
end
