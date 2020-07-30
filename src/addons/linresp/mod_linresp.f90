module mod_linresp
implicit none

! dimension for q-vectors
integer, parameter :: dim_q=2

! type of linear response calculation
!   0 : charge response
!   1 : magnetic response (not yet implemented)
integer lrtype
data lrtype/0/

! temporary array used in genchi0blh
complex(8), allocatable :: megqblh2(:,:)

! interval of bands to include in chi0
real(8) chi0_include_bands(2)
data chi0_include_bands/-100.1d0,100.1d0/

! interval of bands to exclude from chi0
real(8) chi0_exclude_bands(2)
data chi0_exclude_bands/100.1d0,-100.1d0/

! .true. if all k's are used to calculate chi0
logical chi0allk
data chi0allk/.true./

! maximum |k| used in calculation of chi0
real(8) chi0kmax
data chi0kmax/100.0d0/
! k boundary parameters for an oddly shaped region
real(8) chi0kmaxa
data chi0kmaxa/100.0d0/
real(8) chi0kmaxb
data chi0kmaxb/100.0d0/
real(8) chi0kmaxc
data chi0kmaxc/100.0d0/
integer chi0kmaxan
data chi0kmaxan/0/
integer chi0kmaxbn
data chi0kmaxbn/0/

! k offset used if only some k's are used in calculating chi0
real(8) chi0kloff(3)
data chi0kloff/0.0d0,0.0d0,0.0d0/

integer nexwan
integer, allocatable :: iexwan(:)

complex(8), allocatable :: wann_cc(:,:,:)
complex(8), allocatable :: wann_cc2(:,:)

complex(8), allocatable :: chi0loc(:,:,:)
complex(8), allocatable :: chi0wanloc(:,:,:)
complex(8), allocatable :: jdosloc(:)

! number of energy-mesh points
integer lr_nw
data lr_nw/201/
! first energy point (Ha)
real(8) lr_w0
data lr_w0/0.d0/
! last energy point (Ha)
real(8) lr_w1
data lr_w1/1.d0/
! energy step
real(8) lr_dw
! energy mesh
complex(8), allocatable :: lr_w(:)
! broadening parameter (Ha)
real(8) lr_eta
data lr_eta/0.01d0/
! QP correction: 0 for G0W0, 1 for GW0 and 2 for GW
integer gw_mode
data gw_mode/0/
! energy step used in G0W0 calculation
real(8) del_e
data del_e/0.01d0/
! inverse temperature for the matsubara frequency in eV^-1
real(8) lr_beta
data lr_beta/30.d0/
! .true. if imaginary frequency mesh is required
logical timgw
data timgw/.false./
! first imaginary frequency
real(8) lr_iw0
data lr_iw0/0.d0/
! last imaginary frequency
real(8) lr_iw1
data lr_iw1/80.d0/

real(8) fxca0
data fxca0/0.d0/
real(8) fxca1
data fxca1/0.d0/
integer nfxca
data nfxca/1/
integer fxctype
data fxctype/0/

! high-level switch: compute chi0 and chi in Wannier functions basis
logical wannier_chi0_chi 
data wannier_chi0_chi/.false./
! high-level switch: .true. if chi0 should be multiplied by 2
logical wannier_chi0_afm
data wannier_chi0_afm/.false./

! indices of response functions in global array f_response(:,:,:)
integer, parameter :: f_chi0                 = 1
integer, parameter :: f_chi                  = 2
integer, parameter :: f_chi_scalar           = 3
integer, parameter :: f_chi_pseudo_scalar    = 4
integer, parameter :: f_epsilon_matrix_GqGq  = 5
integer, parameter :: f_epsilon_scalar_GqGq  = 6
integer, parameter :: f_inv_epsilon_inv_GqGq = 7
integer, parameter :: f_epsilon_eff          = 8
integer, parameter :: f_epsilon_eff_scalar   = 9
integer, parameter :: f_sigma                = 10
integer, parameter :: f_sigma_scalar         = 11
integer, parameter :: f_loss                 = 12
integer, parameter :: f_loss_scalar          = 13
integer, parameter :: f_chi0_wann            = 14
integer, parameter :: f_chi_wann             = 15
integer, parameter :: f_epsilon_eff_wann     = 16
integer, parameter :: f_sigma_wann           = 17
integer, parameter :: f_loss_wann            = 18
integer, parameter :: f_epsilon_inv_GqGq     = 19
integer, parameter :: f_jdos                 = 20

integer, parameter :: nf_response            = 20
complex(8), allocatable :: f_response(:,:,:)

complex(8), allocatable :: u4(:,:,:,:)
logical screenu4
data screenu4/.true./

logical wann_trans
data wann_trans/.false./

complex(8), allocatable :: gw_self_energy(:,:,:)
complex(8), allocatable :: self_energy_x(:,:)
contains

subroutine genchi0blh(ikloc,ngqidx,w,chi0w,jdosw,omegap)
use modmain
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
! arguments
integer, intent(in) :: ikloc
integer, intent(in) :: ngqidx
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0w(ngqidx,ngqidx)
complex(8), optional, intent(out) :: jdosw
complex(8), optional, intent(inout) :: omegap(3,3)
! local variables
logical l1
integer i,ist1,ist2,ik,jk,ig,i1,i2
real(8) t1,t2,e1,e2,x
complex(8) pvec(3)
complex(8), allocatable :: wt(:)
! external functions
real(8), external :: sdelta
logical, external :: bndint
! 
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
allocate(wt(nmegqblh(ikloc)))
wt(:)=zzero
do i=1,nmegqblh(ikloc)
  ist1=bmegqblh(1,i,ikloc)
  ist2=bmegqblh(2,i,ikloc)
! default : include all interband transitions         
  l1=.true.
! cRPA case : don't include bands in energy window [crpa_e1,crpa_e2]
  if (bndint(ist1,evalsvnr(ist1,ik),chi0_exclude_bands(1),&
      chi0_exclude_bands(2)).and.bndint(ist2,evalsvnr(ist2,jk),&
      chi0_exclude_bands(1),chi0_exclude_bands(2))) l1=.false.
  if (l1) then
    t1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
    e1=evalsvnr(ist1,ik)
    e2=evalsvnr(ist2,jk)
    if (abs(t1).gt.1d-6.and.(ist1.eq.ist2.and.intraband.or.ist1.ne.ist2)) then
      t2=sign(scissor,t1)
      wt(i)=t1/(e1-e2-t2+w)
      if (present(jdosw)) jdosw=jdosw+wt(i)
    else if (ist1.eq.ist2.and.intraband) then
      x=(efermi-e1)/swidth
      t1=occmax*sdelta(stype,x)/swidth
      pvec(:)=real(pmatnrloc(:,ist1,ist1,ikloc))
      do i1=1,3
        do i2=1,3
          omegap(i1,i2)=omegap(i1,i2)+t1*pvec(i1)*pvec(i2)
        enddo
      enddo
      wt(i)=t1*(dot_product(vq0c(:),pvec(:))/w)**2
!      if (mpi_grid_root().and.t1.gt.1d-10) then
!        write(*,*)'ist1,ist2,t1,pmat,w', &
!                 & ist1,ist2,t1,pmatnrloc(1,ist1,ist1,ikloc),w
!      endif
    endif
  endif
enddo !i
call papi_timer_start(pt_megqblh2)
do ig=1,ngqidx
  megqblh2(1:nmegqblh(ikloc),ig)=dconjg(megqblh(1:nmegqblh(ikloc),ig,ikloc))*wt(1:nmegqblh(ikloc))
enddo
!if (mpi_grid_root().and.ikloc.eq.1) then
!  open(367,file='MEGQBLH.OUT',form='formatted',status='replace')
!  do i=1,nmegqblh(ikloc)
!    do ig=1,ngq_
!      write(367,'(2I6,2G14.6)')i,ig,real(megqblh(i,ig,ikloc)),imag(megqblh(i,ig,ikloc))
!    enddo
!    write(367,*)
!  enddo
!  close(367)
!endif
call papi_timer_stop(pt_megqblh2)
call papi_timer_start(pt_chi0_zgemm)
call zgemm('T','N',ngqidx,ngqidx,nmegqblh(ikloc),zone,megqblh(1,1,ikloc),nstsv*nstsv,&
  &megqblh2(1,1),nstsv*nstsv,zone,chi0w(1,1),ngqidx)
call papi_timer_stop(pt_chi0_zgemm)
deallocate(wt)
return
end subroutine

end module
