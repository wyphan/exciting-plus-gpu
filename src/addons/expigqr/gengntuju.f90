subroutine gengntuju(iq,lmaxexp)
use modmain
use mod_addons_q
use mod_expigqr
use mod_util

#ifdef _HDF5_
USE mod_hdf5
#endif /* HDF5 */

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
integer i,j,nuju,nujuloc,i1
real(8), allocatable :: uju(:,:,:,:,:,:)
logical, allocatable :: ujuflg(:,:,:,:,:)

#ifdef _HDF5_

CHARACTER(LEN=100) :: fname, pathq, pathqc, pathqcg
CHARACTER(LEN=2) :: c2
CHARACTER(LEN=3) :: c3
CHARACTER(LEN=4) :: c4
LOGICAL :: exist
REAL(KIND(1.D0)), PARAMETER :: toKiB = 2.D0**(-10)
REAL(KIND(1.D0)), PARAMETER :: toMiB = 2.D0**(-20)
INTEGER :: bytes
INTEGER, DIMENSION(3) :: gntyyydim, gntyyychunk
INTEGER, DIMENSION(6) :: ujudim, ujuchunk
INTEGER, DIMENSION(2) :: gntujudim, gntujuchunk

#else

#if defined(_DUMP_gntyyy_) || defined(_DUMP_uju_) || defined(_DUMP_gntuju_)

CHARACTER(LEN=32) :: refile
CHARACTER(LEN=128) :: cmd

#ifdef _DUMP_gntyyy_
CHARACTER(LEN=20) :: fmt1
#endif /* DUMP_gntyyy */

#ifdef _DUMP_uju_
CHARACTER(LEN=20) :: fmt2
#endif /* DUMP_uju */

#ifdef _DUMP_gntuju_
CHARACTER(LEN=32) :: imfile
CHARACTER(LEN=17) :: fmt3
INTEGER :: irow, icol
#endif /* DUMP_gntuju */

#endif /* DUMP_{gntyyy,uju,gntuju} */

#endif /* HDF5 */

lmmaxexp=(lmaxexp+1)**2
allocate(jl(nrmtmax,0:lmaxexp))
allocate(zm(0:lmaxexp,lmmaxapw,lmmaxapw))
igntuju=0
ngntuju=0
gntuju=zzero

#ifdef _DUMP_gntyyy_
IF( mpi_grid_root() ) THEN

#ifdef _HDF5_
   fname = "gntuju.hdf5"
   INQUIRE( FILE=TRIM(fname), EXIST=exist )

   ! Create file and populate data structure
   IF ( (iq == 1) .AND. (.NOT. exist) ) CALL writegntujuheader( fname )

   ! Dump gntyyy with HDF5
   IF ( iq == 1 ) THEN
      gntyyydim = (/ lmmaxvr, lmmaxapw, lmmaxapw /)
      gntyyychunk(:) = gntyyydim(:)
      bytes = hdf5_calc_chunksz( 'd', 3, gntyyychunk )
      WRITE(*,*) 'Dumping gntyyy (', INT(REAL(bytes)*toKiB), ' KiB)'
      CALL hdf5_gzwrite_array_d( gntyyy(1,1,1), 3, gntyyydim, gntyyychunk, 9, &
           fname, "/", "gntyyy" )
   END IF

#else

   ! Dump gntyyy
   IF( iq == 1 ) THEN
      WRITE(*,*)
      WRITE(*,'("gntyyy( lm2, lm1, lm3 )")')
      WRITE(*,'("lmmaxvr = ",I2)') lmmaxvr
      WRITE(*,'("lmmaxapw = ",I2)') lmmaxapw
      WRITE(*,*)

      OPEN(UNIT = 666, FILE = "gntyyy.csv")
      fmt1 = '(I2,5(",",I2),F10.6)'
      WRITE(666,'("# L1, m1, L2, m2, L3, m3, gntyyy")')
      DO L1 = 0, lmaxapw
         DO m1 = -L1, L1
            lm1 = idxlm(l1,m1)
            DO L2 = 0, lmaxapw
               DO m2 = -L2, L2
                  lm2 = idxlm(l2,m2)
                  DO L3 = 0, lmaxexp
                     DO m3 = -L3, L3
                        lm3 = idxlm(l3, m3)
                        WRITE(666,fmt) L1, m1, L2, m2, L3, m3, &
                                       gntyyy(lm3, lm2, lm1)
                     END DO ! m3
                     WRITE(666,*)
                  END DO ! L3
               END DO ! m2
            END DO ! L2
         END DO ! m1
      END DO ! L1
      CLOSE(666)
      
      WRITE(cmd, '("zip -9 -m -q gntyyy.zip gntyyy.csv)')
      CALL EXECUTE_COMMAND_LINE(TRIM(cmd))
      
   END IF ! iq = 1

#endif /* HDF5 */
END IF ! wproc
#endif /* DUMP_gntyyy */

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

#ifdef _DUMP_uju_
IF( mpi_grid_root() ) THEN
   
#ifdef _HDF5_

   ! For now, only rank 0 writes using serial hdf5
   ! TODO: Use MPI-IO (parallel hdf5)

   fname = "gntuju.hdf5"
   INQUIRE( FILE=TRIM(fname), EXIST=exist )

   ! Create file and populate data structure
   IF ( (iq == 1) .AND. (.NOT. exist) ) CALL writegntujuheader( fname )

   WRITE( c3, '(I3.3)' ) iq
   pathq = "/qpoints/" // c3

   ! Dump uju( l3, l1, l2, io1, io2, (ic-1)*ngqsh(iq)+gqshidx(ig,iq) )
      ujudim = (/ lmaxexp+1, lmaxapw+1, lmaxapw+1, nufrmax, nufrmax, nuju /)
      ujuchunk(:) = ujudim(:)
      bytes = hdf5_calc_chunksz( 'd', 6, ujuchunk )
      WRITE(*,*) 'Dumping uju (', INT(REAL(bytes)*toKiB), ' KiB)'
      CALL hdf5_gzwrite_array_d( uju(0,0,0,1,1,1), 6, ujudim, ujuchunk, 9, &
           fname, pathq, "uju" )

#else

   IF( iq == 1 ) THEN
      WRITE(*,*)
      WRITE(*,'("uju( l3, l1, l2, io1, io2, &
                &(ic-1)*ngqsh(iq)+gqshidx(ig,iq) )")')
      WRITE(*,'("lmaxapw = ",I2)') lmaxapw
      WRITE(*,'("nufrmax = ",I4)') nufrmax
      WRITE(*,'("nuju = ngqsh(iq)*natmcls = ",I4)') nuju
      WRITE(*,*)
   END IF ! iq = 1

   ! Dump uju
   fmt2 = '(I2,4(",",I2),F10.6)'
   DO ig = 1, ngq(iq)
      DO ic = 1, natmcls
         ias = ic2ias(ic) ! atom index
         is = ias2is(ias) ! species index 
         WRITE(refile, '("uju.iq",I2.2,".ig",I4.4,".ic",I2.2,".csv")') &
                        iq, ig, ic
         OPEN(UNIT = 666, FILE = TRIM(refile))
         WRITE(666,'("# L1, L2, io1, io2, L3, uju")')
         DO L1 = 0, lmaxapw
            DO L2 = 0, lmaxapw
               DO io1 = 1, nufr(L1,is)
                  DO io2 = 1, nufr(L2,is)
                     DO L3 = 0, lmaxexp
                        WRITE(666,fmt2) L1, L2, io1, io2, L3, &
                    uju( L3, L1, L2, io1, io2, (ic-1)*ngqsh(iq)+gqshidx(ig,iq) )
                     END DO ! L3
                     WRITE(666,*)
                  END DO ! io2
               END DO ! io1
            END DO ! L2
         END DO ! L1
         CLOSE(666)
         !-- TODO: There should be a better way than this
         WRITE(cmd, '("zip -9 -m -q uju.zip ", A)') TRIM(refile)
         CALL EXECUTE_COMMAND_LINE(TRIM(cmd))
         !--
      END DO ! ic
   END DO ! ig

#endif /* HDF5 */

END IF ! wproc
#endif /* DUMP_uju */

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
                  gntuju(lm2+(io2-1)*lmmaxapw,lm1+(io1-1)*lmmaxapw,ic,ig)=zt1
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
  call mpi_grid_reduce(gntuju(1,1,1,(ig-1)*ngvb+1),ngvb*ngntujumax*ngntujumax*natmcls,&
    &dims=(/dim_k/),all=.true.)
  ! Right now these are unused
!  call mpi_grid_reduce(igntuju(1,1,1,(ig-1)*ngvb+1),ngvb*2*ngntujumax*natmcls,&
!    &dims=(/dim_k/),all=.true.)
!  call mpi_grid_reduce(ngntuju(1,(ig-1)*ngvb+1),ngvb*natmcls,dims=(/dim_k/),&
!    &all=.true.)    
  call mpi_grid_barrier(dims=(/dim_k/))
enddo
do ig=i*ngvb+1,ngq(iq)
  call mpi_grid_reduce(gntuju(1,1,1,ig),ngntujumax*ngntujumax*natmcls,dims=(/dim_k/),&
       &all=.true.)
  ! Right now these are unused
!  call mpi_grid_reduce(igntuju(1,1,1,ig),2*ngntujumax*natmcls,dims=(/dim_k/),&
!    &all=.true.)
!  call mpi_grid_reduce(ngntuju(1,ig),natmcls,dims=(/dim_k/),all=.true.)    
  call mpi_grid_barrier(dims=(/dim_k/))
enddo

#ifdef _DUMP_gntuju_
! Dump gntuju
IF( wproc ) THEN

#ifdef _HDF5_

   ! For now, only rank 0 writes using serial hdf5
   ! TODO: Use MPI-IO (parallel hdf5)

   fname = "gntuju.hdf5"
   INQUIRE( FILE=TRIM(fname), EXIST=exist )

   ! Create file and populate data structure
   IF ( (iq == 1) .AND. (.NOT. exist) ) CALL writegntujuheader( fname )

   WRITE( c3, '(I3.3)' ) iq
   pathq = "/qpoints/" // c3

   ! Dump gntuju( lm2+(io2-1)*lmmaxapw, lm1+(io1-1)*lmmaxapw, ic, ig )
   DO ic = 1, natmcls
      WRITE( c2, '(I2.2)' ) ic

      DO ig = 1, ngq(iq)
         WRITE( c4, '(I4.4)' ) ig
         pathqcg = TRIM(pathq) // "/class/" // c2 // "/gvectors/" // c4

         gntujudim = (/ ngntujumax, ngntujumax /)
         gntujuchunk = gntujudim
         bytes = hdf5_calc_chunksz( 'z', 2, gntujuchunk )
         WRITE(*,*) 'Dumping gntuju(:,:,ic=', ic, ',ig=', ig, ') (', &
                    INT(REAL(bytes)*tokiB), ' kiB)'
         CALL hdf5_gzwrite_array_z( gntuju(1,1,ic,ig), 2, &
                                    gntujudim, gntujuchunk, 9, &
                                    fname, pathqcg, "gntuju" )

      END DO ! ig
   END DO ! ic

#else
   
   IF( iq == 1 ) THEN
      WRITE(*,*)
      WRITE(*,'("gntuju( lm2+(io2-1)*lmmaxapw, lm1+(io1-1)*lmmaxapw, &
                           &ic, ig )")')
      WRITE(*,'("lmaxapw = ",I2,", lmmaxapw = (lmaxapw+1)**2 = ",I4)') &
                 lmaxapw, lmmaxapw
      WRITE(*,'("nufrmax = ",I4)') nufrmax
      WRITE(*,'("ngntujumax = lmmaxapw*nufrmax = ",I4)') ngntujumax
      WRITE(*,*)
   END IF ! iq = 1
   WRITE(*,'("ngq(iq=",I2,")=",I4)') iq, ngq(iq)
   WRITE(fmt3, '("(",I4.4,"(F10.6,'',''))")') ngntujumax
   DO ig = 1, ngq(iq)
      DO ic = 1, natmcls
         WRITE(refile, '("gntuju.iq",I2.2,".ic",I2.2,".ig",I4.4,"_re.csv")') &
              iq, ic, ig
         WRITE(imfile, '("gntuju.iq",I2.2,".ic",I2.2,".ig",I4.4,"_im.csv")') &
              iq, ic, ig
         OPEN(UNIT = 666, FILE = TRIM(refile))
         DO irow = 1, ngntujumax
            WRITE(666,fmt3)( DREAL(gntuju(irow,icol,ic,ig)),icol=1,ngntujumax )
         END DO
         CLOSE(666)
         OPEN(UNIT = 777, FILE = TRIM(imfile))
         DO irow = 1, ngntujumax
            WRITE(777,fmt3)( DIMAG(gntuju(irow,icol,ic,ig)),icol=1,ngntujumax )
         END DO
         CLOSE(777)
         !-- TODO: There should be a better way than this
         WRITE(cmd, '("zip -9 -m -q gntuju.zip ", 2(1X,A))') &
                     TRIM(refile), TRIM(imfile)
         CALL EXECUTE_COMMAND_LINE(TRIM(cmd))
         !--
      END DO ! ic
   END DO ! ig

#endif /* HDF5 */

END IF ! wproc
#endif /* DUMP_gntuju */

DO ig = 1, ngq(iq)
   DO ic = 1, natmcls

      ! Find nonzero rows of gntuju
      igntujunz(:,ic,ig) = 0
      CALL zsy2sp_findnnz( ngntujumax, gntuju(1,1,ic,ig), &
                           nmt(ic,ig), igntujunz(1,ic,ig) )

   END DO ! ic
END DO ! ig

nmtmax = MAXVAL( nmt )
#if EBUG >= 1
  WRITE(*,*) 'gengntuju: iproc=', iproc, ' nmtmax=', nmtmax
#endif /* DEBUG */

ALLOCATE( gntuju_packed( nmtmax, nmtmax, natmcls, ngq(iq) ))
ALLOCATE( nareanz(       nufrmax, natmcls, ngq(iq) ))
ALLOCATE( tblgntujunz(   nufrmax, natmcls, ngq(iq) ))

DO ig = 1, ngq(iq)
   DO ic = 1, natmcls

      ! Pack gntuju(ngntujumax,ngntujumax,ic,ig) into gntuju_packed(nmt,nmt,ic,ig)
      CALL zsy2sp_pack( ngntujumax, gntuju(1,1,ic,ig), &
                        nmtmax, igntujunz(1,ic,ig), &
                        nareanz(1,ic,ig), tblgntujunz(1,ic,ig), &
                        gntuju_packed(1,1,ic,ig) )

   END DO ! ic
END DO ! ig

deallocate(jl)
deallocate(uju)
deallocate(zm)

DEALLOCATE( nareanz )
DEALLOCATE( tblgntujunz )

return
end

!-------------------------------------------------------------------------------
#ifdef _HDF5_

SUBROUTINE writegntujuheader( fname )

  USE modmain
  USE mod_addons_q
  USE mod_expigqr
  USE mod_util
  USE mod_hdf5

  IMPLICIT NONE

  ! Input argument
  CHARACTER(*) :: fname

  ! Internal variables
  CHARACTER(LEN=100) :: pathp, pathq, pathqc
  CHARACTER(LEN=2) :: c2
  CHARACTER(LEN=3) :: c3
  CHARACTER(LEN=4) :: c4
  INTEGER :: iq, ig, ic, lmaxexp, nuju

  CALL hdf5_create_file( fname )

  CALL hdf5_create_group( fname, "/", "parameters" )
  pathp = "/parameters"
  CALL hdf5_write( fname, pathp, "ngridk", ngridk(1), (/3/) )
  CALL hdf5_write( fname, pathp, "lmaxapw", lmaxapw )
  lmaxexp = lmaxvr
  CALL hdf5_write( fname, pathp, "lmaxexp", lmaxexp )
  CALL hdf5_write( fname, pathp, "natmcls", natmcls )
  CALL hdf5_write( fname, pathp, "nspecies", nspecies )
  CALL hdf5_write( fname, pathp, "nvq", nvq )
  CALL hdf5_write( fname, pathp, "gqsh", gqsh )
  CALL hdf5_write( fname, pathp, "gqmax", gqmax )
  CALL hdf5_write( fname, pathp, "nufr", nufr(0,1), (/ lmaxapw+1, nspecies /) )

#ifdef _DUMP_gntyyy_
  CALL hdf5_write( fname, pathp, "lmmaxvr", lmmaxvr )
#endif /* _DUMP_gntyyy_ */

#ifdef _DUMP_gntuju_
  CALL hdf5_write( fname, pathp, "nufrmax", nufrmax )
  CALL hdf5_write( fname, pathp, "lmmaxapw", lmmaxapw )
  CALL hdf5_write( fname, pathp, "ngntujumax", ngntujumax )
#endif /* _DUMP_gntuju_ */

  CALL hdf5_create_group( fname, "/", "qpoints" )
  CALL hdf5_write( fname, "/qpoints", "vqm", vqm(1,1), (/ 3, nvq /) )

  DO iq = 1, nvq

     WRITE( c3, '(I3.3)' ) iq
     CALL hdf5_create_group( fname, "/qpoints", c3 )
     pathq = "/qpoints/" // c3

     CALL hdf5_write( fname, pathq, "ngq", ngq(iq) )
     CALL hdf5_write( fname, pathq, "ngqsh", ngqsh(iq) )
     CALL hdf5_write( fname, pathq, "gqshidx", gqshidx(1,iq), (/ ngq(iq) /) )
     nuju=ngqsh(iq)*natmcls
     CALL hdf5_write( fname, pathq, "nuju", nuju )

#ifdef _DUMP_gntuju_

     CALL hdf5_create_group( fname, pathq, "class" )
     DO ic = 1, natmcls
        WRITE( c2, '(I2.2)' ) ic
        CALL hdf5_create_group( fname, TRIM(pathq)//"/class", c2 )
        pathqc = TRIM(pathq) // "/class/" // c2
        CALL hdf5_create_group( fname, TRIM(pathqc), "gvectors" )
        DO ig = 1, ngq(iq)
           WRITE( c4, '(I4.4)' ) ig
           CALL hdf5_create_group( fname, TRIM(pathqc)//"/gvectors", c4 )
        END DO ! ig
     END DO ! ic

#endif /* DUMP_gntuju */

  END DO ! iq

  RETURN

END SUBROUTINE writegntujuheader

#endif /* HDF5 */

!-------------------------------------------------------------------------------
! Finds the non-zero rows of a symmetric double complex matrix mat(nrows,nrows)
! as nrownz, and fills in the translation table icolnz(nrownz)
! TODO: OpenACC port

SUBROUTINE zsy2sp_findnnz( nrows, mat, nrownz, icolnz )

  USE IEEE_ARITHMETIC, ONLY :: IEEE_CLASS, IEEE_CLASS_TYPE, IEEE_IS_FINITE
  USE mod_prec, ONLY: dd, dz

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: nrows
  COMPLEX(KIND=dz), DIMENSION(nrows,nrows), INTENT(IN) :: mat
  INTEGER, INTENT(OUT) :: nrownz
  INTEGER, DIMENSION(nrownz), INTENT(OUT) :: icolnz

  ! Internal variables
  INTEGER, DIMENSION(nrows,nrows) :: tblnz
  INTEGER, DIMENSION(nrows) :: col
  REAL(KIND=dd) :: val
  TYPE(IEEE_CLASS_TYPE) :: ieeeclass
  INTEGER :: i, j, logval

  ! Find nonzero rows on the lower triangular matrix
  ! TODO: Parallelize
  tblnz(:,:) = 0
  DO j = 1, nrows
     DO i = j, nrows

        ! Take the logarithm
        val = LOG10( ABS( mat(i,j) ))

        ! Skip infinities from log10( 0.0 ) = -Inf
        ieeeclass = IEEE_CLASS( val )
        IF( IEEE_IS_FINITE( ieeeclass ) ) tblnz(i,j) = FLOOR( val )

     END DO ! i
  END DO ! j

  ! Sum over columns in tblnz
  ! TODO: Parallelize
  col(:) = 0
  DO j = 1, nrows

     logval = 0
     DO i = j, nrows
        logval = logval + tblnz(i,j)
     END DO ! i

     col(j) = logval

  END DO ! j

  ! Finally, count nrownz and
  ! fill in colnz, the translation table from symmetric to sparse
  ! TODO: Parallelize
  nrownz = 0
  icolnz(:) = 0
  DO j = 1, nrows
     IF( col(j) /= 0 ) THEN
        nrownz = nrownz + 1
        icolnz(nrownz) = j
     END IF
  END DO

  RETURN
END SUBROUTINE zsy2sp_countnnz

!-------------------------------------------------------------------------------
! Permutes a sparse symmetric double complex matrix mat(nrows,nrows) such that
! only the first nrownz rows and columns are filled in matnz(nrownz,nrownz).
! The number of distinct non-zero areas is output to nareanz,
! and the start indices of each area is output to tblcolnz

SUBROUTINE zsy2sp_pack( nrows, mat, nrownz, icolnz, nareanz, tblcolnz, matnz )

  USE mod_prec, ONLY: dz

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: nrows, nrownz
  COMPLEX(KIND=dz), DIMENSION(nrows,nrows), INTENT(IN) :: mat
  INTEGER, DIMENSION(nrownz), INTENT(IN) :: icolnz
  INTEGER, INTENT(OUT) :: nareanz
  INTEGER, DIMENSION(0:nareanz), INTENT(OUT) :: tblcolnz
  COMPLEX(KIND=dz), DIMENSION(nrownz,nrownz), INTENT(OUT) :: matnz

  ! Internal variables
  INTEGER :: i, j, irow, iarea
  LOGICAL :: toggle

  ! Count number of distinct non-zero areas
  nareanz = 1
  irow = 0
  tblcolnz(:) = 0
  toggle = .FALSE.
  tblcolnz(0) = 1
  DO i = 1, nrownz

     IF( icolnz(i) /= irow ) THEN

        toggle = .TRUE.
        nareanz = nareanz + 1
        tblcolnz(nareanz) = icolnz(i)
        irow = icolnz(i)

     ELSE

        toggle = .FALSE.
        icol = icol + 1

     END IF ! icol
  END DO ! i

#if EBUG >= 3
  IF( toggle ) WRITE(*,*) 'zsy2sp_pack: iproc=', iproc, ' narea=', narea, ' irow=', irow
#endif /* DEBUG */

  ! Permute the data
  IF( nareanz > 1 ) THEN
     DO iarea = 1, nareanz

        DO icol = tblcolnz(iarea-1), tblcolnz(iarea)

           j = icolnz(icol)
           DO irow = tblcolnz(iarea-1), tblcolnz(iarea)

              i = icolnz(irow)
              matnz(irow,icol) = mat(i,j)

           END DO ! irow
        END DO ! icol
     END DO ! iarea
  END IF ! nareanz

  RETURN
END SUBROUTINE zsy2sp_pack
