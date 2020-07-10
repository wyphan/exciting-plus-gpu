subroutine wann_plot_3d
use modmain
use mod_nrkp
implicit none

real(8) orig(3)
complex(4), allocatable :: wf(:,:,:)
real(4), allocatable :: wf2(:)
integer i,j,nrtot
integer i1,i2,i3,ir,n,ispn
integer nrxloc,i1loc,ikloc
real(8) vrc(3)
real(4), allocatable :: veff(:)
real(4), allocatable :: func(:,:)
real(8) t1,dens(1)
complex(8) zval(2)
complex(8), allocatable :: zfft(:)
complex(8), allocatable :: funcz(:,:)
character*40 fname
logical, parameter :: wfprod=.false.
integer recl
integer, parameter :: float_size = 4

call init0
call init1
if (.not.mpi_grid_in()) return
wproc=mpi_grid_root()

! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp

call genwfnr(-1,.false.)
call elk_m_init

if (task.eq.863) then
  nrtot=nrxyz(1)*nrxyz(2)*nrxyz(3)
  orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2)+bound3d(:,3))/2.d0
endif

allocate(func(nrxyz(1),nrxyz(2)*nrxyz(3))) 
allocate(funcz(nrxyz(1),nrxyz(2)*nrxyz(3))) 

! Fourier transform potential to G-space
allocate(zfft(ngrtot))
zfft(:)=veffir(:)
call zfftifc(3,ngrid,-1,zfft)

recl=float_size*nrxyz(2)*nrxyz(3)

nrxloc=mpi_grid_map(nrxyz(1),dim2)

!if (mpi_grid_root()) write(*,*)'OK1'

do j=1,nwfplot
  n=j+firstwf-1 
!  call elk_load_wann_unk(n)
  write(fname,'("wf_",I3.3,".dx")')n
  if (wproc) then
    open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
    write(70,'("object 1 class gridpositions counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
    write(70,'("origin ",3G18.10)')orig(:)
    do i=1,3
      write(70,'("delta ",3G18.10)')bound3d(:,i)/(nrxyz(i)-1)
    enddo
    write(70,'("object 2 class gridconnections counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
    write(70,'("object 3 class array type float rank 1 shape 1 items ",&
      &I8," lsb ieee data file wf",I3.3,".bin,0")')nrxyz(1)*nrxyz(2)*nrxyz(3),n
    write(70,'("object ""regular positions regular connections"" class field")')
    write(70,'("component ""positions"" value 1")')
    write(70,'("component ""connections"" value 2")')
    write(70,'("component ""data"" value 3")')
    write(70,'("end")')
    close(70)
  endif
  write(fname,'("wf",I3.3,".bin")')n
  if (wproc) then
    open(70,file=trim(fname),form="unformatted",status="replace",access="direct",recl=recl)
  endif
  funcz=zzero
  do ikloc=1,nkptnrloc
    if (mpi_grid_root()) write(*,'("kpt ",I4," out of ",I4)')ikloc,nkptnrloc
    do i1loc=1,nrxloc
      i1=mpi_grid_map(nrxyz(1),dim2,loc=i1loc)-1
!    if (mpi_grid_root()) write(*,'("slab ",I4," out of ",I4)')i1+1,nrxyz(1)
      ir=0
!    func=0.0
      do i2=0,nrxyz(2)-1
        do i3=0,nrxyz(3)-1
          vrc(:)=orig(:)+i1*bound3d(:,1)/(nrxyz(1)-1)+&
                         i2*bound3d(:,2)/(nrxyz(2)-1)+&
                         i3*bound3d(:,3)/(nrxyz(3)-1)
          call get_wanval(vrc,zval,n,ikloc)
!        dens=0.d0
          ir=ir+1
          do ispn=1,nspinor
            funcz(i1+1,ir)=funcz(i1+1,ir)+zval(ispn)
!          dens(1)=dens(1)+abs(zval(ispn))**2
!          dens(1)=dens(1)+sign(abs(zval(ispn)),real(zval(ispn)))
          enddo
!        ir=ir+1
!        func(i1+1,ir)=func(i1+1,ir)+sngl(dens(1))
        enddo
      enddo
    !call mpi_grid_reduce(wf(1,1,1),nspinor*nwfplot*nrxyz(2)*nrxyz(3),dims=(/dim_k/))
    enddo
  !call mpi_grid_barrier()
!  call mpi_grid_reduce(func(1,1),nrxyz(1)*nrxyz(2)*nrxyz(3),dims=(/dim2/))
  enddo
  call mpi_grid_reduce(funcz(1,1),nrxyz(1)*nrxyz(2)*nrxyz(3),dims=(/dim_k,dim2/))
  do i1=1,nrxyz(1)
    ir=0
    do i2=1,nrxyz(2)
      do i3=1,nrxyz(3)
        ir=ir+1
        func(i1,ir)=sign(abs(funcz(i1,ir)),real(funcz(i1,ir)))
      enddo
    enddo
  enddo
!  call mpi_grid_reduce(func(1,1),nrxyz(1)*nrxyz(2)*nrxyz(3),dims=(/dim_k,dim2/))
  if (wproc) then
    do i1=1,nrxyz(1)
      write(70,rec=i1)func(i1,:)
    enddo
  endif
  close(70)
enddo !j
return
end
