subroutine gengntuju(iq,lmaxexp)
use modmain
use mod_addons_q
use mod_expigqr
use mod_util
implicit none
! arguments
integer, intent(in) :: iq
integer, intent(in) :: lmaxexp

integer ig,is,ir,n,ias,io1,io2,l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,ic
! for parallel
integer igloc,ngqloc
integer lmmaxexp
complex(8) zt1
real(8) fr(nrmtmax)
real(8), allocatable :: jl(:,:)
real(8), external :: rfinteg
complex(8), allocatable :: zm(:,:,:)
integer, parameter :: ngvb=2
integer i,j,nuju,nujuloc,i1,imat
real(8), allocatable :: uju(:,:,:,:,:,:)
logical, allocatable :: ujuflg(:,:,:,:,:)

#ifdef _DUMPGNTUJU_
character(LEN=32) :: refile, imfile
character(LEN=11) :: fmt
character(LEN=128) :: cmd
integer :: irow, icol
#endif

lmmaxexp=(lmaxexp+1)**2
allocate(jl(nrmtmax,0:lmaxexp))
allocate(zm(0:lmaxexp,lmmaxapw,lmmaxapw))
igntuju=0
ngntuju=0
gntuju=zzero

!-- TODO: Remove gntujutmp
gntujutmp=zzero
!--

! total number of MT groups of radial integrals
nuju=ngqsh(iq)*natmcls
nujuloc=mpi_grid_map(nuju,dim_k)
allocate(uju(0:lmaxexp,0:lmaxapw,0:lmaxapw,nufrmax,nufrmax,nuju))
allocate(ujuflg(0:lmaxexp,0:lmaxapw,0:lmaxapw,nufrmax,nufrmax))
uju=0.d0
do i=1,nujuloc
  i1=mpi_grid_map(nuju,dim_k,loc=i)
  ic=int((i1-1)/ngqsh(iq))+1
  j=mod(i1-1,ngqsh(iq))+1
  ias=ic2ias(ic)
  is=ias2is(ias)
! generate Bessel functions j_l(|G+q|x)
  do ir=1,nrmt(is)
    call sbessel(lmaxexp,gqshlen(j,iq)*spr(ir,is),jl(ir,:))
  enddo
  ujuflg=.false.
! compute radial integrals <u_{l1,io1} | j_{l3}(|G+q|x) | u_{l2,io2}>
  do l3=0,lmaxexp
    do l1=0,lmaxapw
      do l2=0,lmaxapw
        do io1=1,nufr(l1,is)
          do io2=1,nufr(l2,is)
            if (.not.ujuflg(l3,l1,l2,io1,io2)) then
              do ir=1,nrmt(is)
                fr(ir)=ufr(ir,l1,io1,ic)*ufr(ir,l2,io2,ic)*jl(ir,l3)
              enddo !ir
              uju(l3,l1,l2,io1,io2,i1)=rintegrate(nrmt(is),spr(1,is),fr)
              uju(l3,l2,l1,io2,io1,i1)=uju(l3,l1,l2,io1,io2,i1)
              ujuflg(l3,l1,l2,io1,io2)=.true.
              ujuflg(l3,l2,l1,io2,io1)=.true.
            endif
          enddo !io2
        enddo !io1
      enddo !l2
    enddo !l1
  enddo !l3
enddo 
deallocate(ujuflg)
call mpi_grid_reduce(uju(0,0,0,1,1,1),(lmaxexp+1)*(lmaxapw+1)*(lmaxapw+1)*nufrmax*nufrmax*nuju,dims=(/dim_k/),all=.true.)

ngqloc=mpi_grid_map(ngq(iq),dim_k)
do igloc=1,ngqloc
  ig=mpi_grid_map(ngq(iq),dim_k,loc=igloc)
! precompute atom-independent array
  do l1=0,lmaxapw
    do m1=-l1,l1 
      lm1=idxlm(l1,m1)
      do l2=0,lmaxapw
        do m2=-l2,l2
          lm2=idxlm(l2,m2)
          do l3=0,lmaxexp
            zt1=zzero
            do lm3=l3**2+1,(l3+1)**2
              zt1=zt1+gntyyy(lm3,lm2,lm1)*ylmgq(lm3,ig)
            enddo !m3
            zm(l3,lm2,lm1)=zt1*fourpi*dconjg(zi**l3)
          enddo !l3
        enddo
      enddo
    enddo
  enddo
! loop over atom classes
  do ic=1,natmcls
    ias=ic2ias(ic)
    is=ias2is(ias)
! compute muffin-tin integrals
!  1) sfacgq and ylmgq are generated for exp^{+i(G+q)x}
!     expansion of a plane-wave: 
!       exp^{+igx}=4\pi \sum_{l_3 m_3} i^{l_3} j_{l_3}(gr)Y_{l_3 m_3}^{*}(\hat g)Y_{l_3 m_3}(\hat r)
!     but we need exp^{-i(G+q)x}, so expansion terms will be conjugated
!  2) angular part of integral:
!     <Y_{l_1 m_1} | e^{-i{G+x}x} | Y_{l_2 m_2}> =
!       = \int d \Omega Y_{l_1 m_1}^{*}Y_{l_3 m_3}^{*} Y_{l_2 m_2} = gaunt coeff, which is real
!     so we can conjugate the integral:
!     \int d \Omega Y_{l_1 m_1} Y_{l_3 m_3} Y_{l_2 m_2}^{*} = gaunt(lm2,lm3,lm1)
!  2*) precomputed gaunt array has different order of indices: gntyyy(lm3,lm2,lm1)=gaunt(lm2,lm3,lm1)
!  3) we can sum over lm3 index of a plane-wave expansion
!  3*) we can sum over m3 index and get atom-idependent array
!  4) structure factor sfacgq of a (G+q) plane wave is taken into 
!     account in genmegqblh subroutine; this allows to keep radial integrals
!     for atom classes only (not for all atoms)
    do l1=0,lmaxapw
      do m1=-l1,l1 
        lm1=idxlm(l1,m1)
        do l2=0,lmaxapw
          do m2=-l2,l2
            lm2=idxlm(l2,m2)
            do io1=1,nufr(l1,is)
              do io2=1,nufr(l2,is)
                zt1=zzero
                do l3=0,lmaxexp
                  zt1=zt1+zm(l3,lm2,lm1)*uju(l3,l1,l2,io1,io2,(ic-1)*ngqsh(iq)+gqshidx(ig,iq)) 
                enddo !l3
                !if (abs(zt1).gt.1d-12) then
                  !ngntuju(ic,ig)=ngntuju(ic,ig)+1
                  !n=ngntuju(ic,ig)

                  !gntuju(lm2+(io2-1)*lmmaxapw,lm1+(io1-1)*lmmaxapw,ic,ig)=zt1

!-- TODO: Remove gntujutmp
                gntujutmp(lm2+(io2-1)*lmmaxapw,lm1+(io1-1)*lmmaxapw,ic,ig)=zt1
!--
                
                !igntuju(1,n,ic,ig)=lm1+(io1-1)*lmmaxapw
                  !igntuju(2,n,ic,ig)=lm2+(io2-1)*lmmaxapw
                !endif
              enddo !io2
            enddo !io1
          enddo !m2
        enddo !l2
      enddo !m1
    enddo !l1
  enddo !ic
enddo !ig
! syncronize all values along auxiliary k-direction
!call mpi_grid_reduce(gntuju(1,1,1),ngntujumax*natmcls*ngvecme,dims=(/dim_k/),all=.true.)
!call mpi_grid_barrier(dims=(/dim_k/))
!call mpi_grid_reduce(igntuju(1,1,1,1),2*ngntujumax*natmcls*ngvecme,dims=(/dim_k/),all=.true.)
!call mpi_grid_barrier(dims=(/dim_k/))
!call mpi_grid_reduce(ngntuju(1,1),natmcls*ngvecme,dims=(/dim_k/),all=.true.)    
!call mpi_grid_barrier(dims=(/dim_k/))
!do ig=1,ngvecme
!  call mpi_grid_reduce(gntuju(1,1,ig),ngntujumax*natmcls,dims=(/dim_k/),&
!    all=.true.)
!  call mpi_grid_reduce(igntuju(1,1,1,ig),2*ngntujumax*natmcls,dims=(/dim_k/),&
!    all=.true.)
!  call mpi_grid_reduce(ngntuju(1,ig),natmcls,dims=(/dim_k/),all=.true.)    
!  call mpi_grid_barrier(dims=(/dim_k/))
!enddo

! synchronize blocks of G-vectors arrays
i=ngq(iq)/ngvb
do ig=1,i

  !call mpi_grid_reduce(gntuju(1,1,1,(ig-1)*ngvb+1),ngvb*ngntujumax*ngntujumax*natmcls,&
  !  &dims=(/dim_k/),all=.true.)

!-- TODO: Remove gntujutmp
   call mpi_grid_reduce( gntujutmp( 1, 1, 1, (ig-1)*ngvb + 1 ), &
                         ngvb*ngntujumax*ngntujumax*natmcls, &
                         dims=(/dim_k/), all=.true. )
!--

  call mpi_grid_reduce(igntuju(1,1,1,(ig-1)*ngvb+1),ngvb*2*ngntujumax*natmcls,&
    &dims=(/dim_k/),all=.true.)
  call mpi_grid_reduce(ngntuju(1,(ig-1)*ngvb+1),ngvb*natmcls,dims=(/dim_k/),&
    &all=.true.)    
  call mpi_grid_barrier(dims=(/dim_k/))
enddo
do ig=i*ngvb+1,ngq(iq)

   !call mpi_grid_reduce(gntuju(1,1,1,ig),ngntujumax*ngntujumax*natmcls,dims=(/dim_k/),&
   !&all=.true.)

!-- TODO: Remove gntujutmp
   call mpi_grid_reduce( gntujutmp(1, 1, 1, ig), ngntujumax*ngntujumax*natmcls,&
                         dims=(/dim_k/), all=.true. )
!--

  call mpi_grid_reduce(igntuju(1,1,1,ig),2*ngntujumax*natmcls,dims=(/dim_k/),&
    &all=.true.)
  call mpi_grid_reduce(ngntuju(1,ig),natmcls,dims=(/dim_k/),all=.true.)    
  call mpi_grid_barrier(dims=(/dim_k/))
enddo

#ifdef _DUMPGNTUJU_
! Dump gntujutmp
IF( wproc ) THEN
   IF( iq == 1 ) THEN
      WRITE(*,*)
      WRITE(*,'("gntujutmp( lm2+(io2-1)*lmmaxapw, lm1+(io1-1)*lmmaxapw, &
                           &ic, ig )")')
      WRITE(*,'("lmaxapw = ",I2,", lmmaxapw = (lmaxapw+1)**2 = ",I4)') &
                 lmaxapw, lmmaxapw
      WRITE(*,'("nufrmax = ",I4)') nufrmax
      WRITE(*,'("ngntujumax = lmmaxapw*nufrmax = ",I4)') ngntujumax
      WRITE(*,*)
   END IF
   WRITE(fmt, '("(",I4.4,"F10.6)")') ngntujumax
   DO ig = 1, ngq(iq)
      DO ic = 1, natmcls
         WRITE(refile, '("gntuju.iq",I2.2,".ic",I2.2,".ig",I4.4,"_re.csv")') &
              iq, ic, ig
         WRITE(imfile, '("gntuju.iq",I2.2,".ic",I2.2,".ig",I4.4,"_im.csv")') &
              iq, ic, ig
         OPEN(UNIT = 666, FILE = TRIM(refile))
         DO irow = 1, ngntujumax
            WRITE(666,fmt)( DREAL(gntujutmp(irow,icol,ic,ig)),icol=1,ngntujumax)
         END DO
         CLOSE(666)
         OPEN(UNIT = 777, FILE = TRIM(imfile))
         DO irow = 1, ngntujumax
            WRITE(777,fmt)( DIMAG(gntujutmp(irow,icol,ic,ig)),icol=1,ngntujumax)
         END DO
         CLOSE(777)
         !-- TODO: There should be a better way than this
         WRITE(cmd, '("zip -9 -m gntuju.zip", 2(1X,A))') &
                     TRIM(refile), TRIM(imfile)
         CALL EXECUTE_COMMAND_LINE(TRIM(cmd))
         !--
      END DO ! ic
   END DO ! ig
END IF
#endif

!-- TODO: Get rid of this hack
! Rearrange gntujutmp( ngntujumax, ngntujumax, natmcls, ngq(iq) ) into
!           gntuju(    ngntujumax, ngntujumax*ngq(iq), natmcls  )
do ic = 1, natmcls
   do ig = 1, ngq(iq)
      do imat = 1, ngntujumax
         gntuju( :, (ig-1)*ngntujumax +imat, ic ) = gntujutmp( :, imat, ic, ig )
      end do ! imat
   end do ! ig
end do ! ic
!--

deallocate(jl)
deallocate(uju)
deallocate(zm)
return
end
