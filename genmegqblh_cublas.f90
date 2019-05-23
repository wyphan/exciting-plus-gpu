subroutine genmegqblh(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
  &wfsvit1,wfsvit2)
use modmain
use mod_addons_q
use mod_nrkp
use mod_expigqr
use ISO_C_BINDING
use cublas_f
implicit none
integer, intent(in) :: iq
integer, intent(in) :: ikloc
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvmt1(lmmaxapw*nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(lmmaxapw,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)

integer wfsize
integer ivg1(3)
integer i,j,ik,jk,igkq,n1,ispn1,ispn2,ist1,ist2,ic
integer ig,ig1,ig2,ias,ifg,ir
logical l1, wprocrank
complex(8), allocatable :: wftmp1(:,:)
complex(8), allocatable :: wftmp2(:,:)
complex(8), allocatable :: wfir1(:)
complex(8) b1(lmmaxapw*nufrmax),b2(lmmaxapw*nufrmax)

complex(8),dimension(:),pointer :: b1Batch, b2Batch
complex(8),dimension(:),pointer :: gntujuBatch
type(C_PTR) :: d_b1, d_b2, d_gntuju
type(C_PTR), dimension(:), pointer :: h_d_b1, h_d_b2, h_d_gntuju
type(C_PTR) :: handle,stream
integer :: sizeof_complex,sizeof_ptr,idx,idx2,stat,batch_count
integer(8) :: b1Size,b2Size,gntujuSize,bytes,address
parameter (sizeof_complex=16,sizeof_ptr=8)
wfsize=lmmaxapw*nufrmax*natmtot+ngknr2

allocate(wftmp1(wfsize,ngq(iq)))
allocate(wftmp2(wfsize,nstsv))
allocate(wfir1(ngrtot))
call papi_timer_start(pt_megqblh)

wprocrank=.false.
if (mpi_grid_root((/dim_k/))) then
  wprocrank=.true.
endif

!call mpi_world_barrier

!if(wprocrank) then
!write(*,*) 'wprocrank is ',wprocrank
!flush(6)
!endif
!write(*,*) 'Entered into genmegqblh_cublas.f90 and have allocated wftmp1, wftmp2, and wfir1'
!flush(6)
!endif
! global k-point
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! jk=k+q-G_q
jk=idxkq(1,ik)
! G_q vector 
igkq=idxkq(2,ik)

batch_count = ngq(iq)*natmtot

!allocate(gntujuBatch(lmmaxapw*nufrmax, lmmaxapw*nufrmax, batch_count))
allocate(b1Batch(lmmaxapw*nufrmax*batch_count))
allocate(b2Batch(lmmaxapw*nufrmax*batch_count))
allocate(gntujuBatch(lmmaxapw*nufrmax*lmmaxapw*nufrmax*batch_count))

!if(wprocrank) then
!write(*,*) '  Allocated the b1Batch, b2Batch, and gntujuBatch arrays.'
!flush(6)
!endif
!call mpi_world_barrier
allocate(h_d_b1(batch_count))
allocate(h_d_b2(batch_count))
allocate(h_d_gntuju(batch_count))

!if(wprocrank) then
!write(*,*) '  Allocated the h_d_b1, h_d_b2, and h_d_gntuju arrays.'
!flush(6)
!endif
!call mpi_world_barrier
b1Size = lmmaxapw*nufrmax*batch_count*sizeof_complex
b2Size = b1Size
gntujuSize = lmmaxapw*nufrmax*b1Size

!if(wprocrank) then
!  write(*,*) '  b1Size is ', b1Size
!  write(*,*) '  b2Size is ', b2Size
!  write(*,*) '  gntujuSize is ', gntujuSize
!  flush(6)
!endif
!call mpi_world_barrier
!if (wprocrank) then
!write(*,*) 'Allocated the local batch arrays and the device pointer arrays'
!flush(6)
!endif
!Allocate space on the device to hold the batch arrays as a contiguous memory block
stat = cudaMalloc(h_d_b1(1), b1Size)
stat = cudaMalloc(h_d_b2(1), b2Size)
stat = cudaMalloc(h_d_gntuju(1), gntujuSize)

!if(wprocrank) then
!  write(*,*) '  cudaMalloced h_d_b1(1) with address  '
!  i = printValue( C_LOC(h_d_b1(1)) )
!  !address = C_LOC( h_d_b1(1) )
!  write(*,*) '  cudaMalloced h_d_b2(1) with address  '
!  i = printValue( C_LOC(h_d_b2(1)) )
!  !printValue( C_LOC(h_d_b2(1)) )
!  write(*,*) '  cudaMalloced h_d_gntuju(1) with address  '
!  i = printValue( C_LOC(h_d_gntuju(1)) )
  !printValue( C_LOC(h_d_gntuju(1)) )
!  flush(6)
!endif
!call mpi_world_barrier
!if(wprocrank) then
!  write(*,*) '  Before the loop to add the offset:'
!  do i=1,batch_count
!    write(*,*) ' ITERATION ',i
!    write(*,*) '    h_d_b1(i) = '
!    j = printValue( C_LOC(h_d_b1(i)) )
!    write(*,*) '    h_d_b2(i) = '
!    j = printValue( C_LOC(h_d_b2(i)) )
!    write(*,*) '    h_d_gntuju(i) = '
!    j = printValue( C_LOC(h_d_gntuju(i)) )
!  enddo
!  flush(6)
!endif
!call mpi_world_barrier
do i=2,batch_count
  !h_d_b1(i) = h_d_b1(i) + b1Size
  !h_d_b2(i) = h_d_b2(i) + b2Size
  !h_d_gntuju(i) = h_d_gntuju(i) + gntujuSize

   h_d_b1(i) = h_d_b1(i-1)
   h_d_b2(i) = h_d_b2(i-1)
   h_d_gntuju(i) = h_d_gntuju(i-1)

   stat = addOffsetToPtr( C_LOC(h_d_b1(i)), b1Size )
   stat = addOffsetToPtr( C_LOC(h_d_b2(i)), b2Size )
   stat = addOffsetToPtr( C_LOC(h_d_gntuju(i)), gntujuSize )
enddo
!call mpi_world_barrier

!if(wprocrank) then
!  write(*,*) '  After the loop to add the offset:'
!  do i=1,batch_count
!    write(*,*) 'ITERATION ',i
!    write(*,*) '    h_d_b1(i) = '
!    j = printValue( C_LOC(h_d_b1(i)) )
!    write(*,*) '    h_d_b2(i) = '
!    j = printValue( C_LOC(h_d_b2(i)) )
!    write(*,*) '    h_d_gntuju(i) = '
!    j = printValue( C_LOC(h_d_gntuju(i)) )
!  enddo
!  flush(6)
!endif
!call mpi_world_barrier
!if (wprocrank) then
!write(*,*) 'cudaMalloced the device ptr arrays'
!flush(6)
!endif

!Now we copy the host gntuju array into a contiguous block.
gntujuSize = lmmaxapw*lmmaxapw*nufrmax*nufrmax*sizeof_complex
do ig=1,ngq(iq)
  do ias=1,natmtot
    idx= ( (ig-1)*natmtot + ias - 1 )*(lmmaxapw*lmmaxapw*nufrmax*nufrmax) + 1
    ic = ias2ic(ias)

    call memcopy(gntuju(1,1,ic,ig), gntujuBatch(idx),gntujuSize)
  enddo
enddo
gntujuSize = lmmaxapw*nufrmax*b1Size
!call mpi_world_barrier
!if(wprocrank) then
!  write(*,*) ' Copied the gntuju batches into gntujuBatch array.'
!  flush(6)
!endif
!call mpi_world_barrier
bytes = batch_count*sizeof_ptr

stat = cudaMalloc(d_b1, bytes)
stat = cudaMalloc(d_b2, bytes)
stat = cudaMalloc(d_gntuju, bytes)

!if(wprocrank) then
!  write(*,*) ' cudaMalloced the device pointer arrays.'
!  flush(6)
!endif
!call mpi_world_barrier
stat = cudaMemcpy(d_b1, C_LOC(h_d_b1(1)),bytes,cudaMemcpyHostToDevice); 
stat = cudaMemcpy(d_b2, C_LOC(h_d_b2(1)),bytes,cudaMemcpyHostToDevice); 
stat = cudaMemcpy(d_gntuju, C_LOC(h_d_gntuju(1)),bytes,cudaMemcpyHostToDevice); 

!if(wprocrank) then
!  write(*,*) ' cudaMemcpyed the device pointer arrays over th GPU'
!  flush(6)
!endif
!call mpi_world_barrier
!if (wprocrank) then
!write(*,*) 'cudaMalloced the arrays and copied over the batch array ptrs'
!flush(6)
!endif
stat = cublasCreate(handle)
stat = cudaStreamCreate(stream)

! Transfer gntujuBatch array once.
stat = cudaMemcpy(h_d_gntuju(1), C_LOC(gntujuBatch(1)), gntujuSize, cudaMemcpyHostToDevice)

!if(wprocrank) then
!  write(*,*) ' Copied the gntujuBatch array over to GPU.'
!  write(*,*) ' sent over gntujuSize bytes of data: ', gntujuSize
!  flush(6)
!endif
!call mpi_world_barrier
!if (wprocrank) then
!write(*,*) 'ngkmax=',ngkmax
!write(*,*) 'nstsv=',nstsv
!write(*,*) 'iq=',iq
!write(*,*) 'ngq(iq)=',ngq(iq)
!write(*,*) 'ngrtot=',ngrtot
!write(*,*) 'lmmaxapw=',lmmaxapw
!write(*,*) 'nufrmax=',nufrmax
!write(*,*) 'natmtot=',natmtot
!write(*,*) 'wfsize=',wfsize
!write(*,*) 'batch_count=',batch_count
!write(*,*) 'b1Size=',b1Size
!write(*,*) 'gntujuSize=',gntujuSize
!write(*,*) 'bytes=',bytes
!write(*,*) 'nspinor=',nspinor
!write(*,*) 'ig=1..',ngq(iq)
!write(*,*) 'ias=1,..',natmtot
!endif
do ispn1=1,nspinor
  if (expigqr22.eq.1) ispn2=ispn1
! index of the interband transitions
  i=1
! go through the interband transitions    
  do while (i.le.nmegqblh(ikloc))
! left <bra| state 
    ist1=bmegqblh(1,i,ikloc)
    wftmp1=zzero
    l1=.true.
    if (spinpol) then
      if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
    endif
    if (l1) then
      call timer_start(3)
      call papi_timer_start(pt_megqblh_mt)
      do ig=1,ngq(iq)
! precompute muffint-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
        do ias=1,natmtot
          idx = ( (ig-1)*natmtot + ias - 1 )*(lmmaxapw*nufrmax) + 1
          do j=1,lmmaxapw*nufrmax
            b1Batch(idx+j)=dconjg(wfsvmt1(j,ias,ispn1,ist1)*sfacgq(ig,ias))
          enddo
         !ic=ias2ic(ias)
          !b2Batch(:,1,idx)=zzero
          !b2=zzero
          !do j=1,ngntuju(ic,ig)
          !  b2(igntuju(2,j,ic,ig))=b2(igntuju(2,j,ic,ig))+&
          !    &b1(igntuju(1,j,ic,ig))*gntuju(j,ic,ig)
          !enddo

          !call memcopy(b1(1), b1Batch(1,1,idx), b1Size) 
          !call memcopy(gntuju(1,1,ic,ig), gntujuBatch(1,1,idx),gntujuSize) 
          !stat  = cublasSetMatrixAsync(lmmaxapw*nufrmax,1,sizeof_complex,&
          ! &C_LOC(b1Batch(1,1,idx)),lmmaxapw*nufrmax, h_d_b1(idx), &
          ! &lmmaxapw*nufrmax, stream)

          !stat = cublasSetMatrixAsync(lmmaxapw*nufrmax,lmmaxapw*nufrmax,&
          ! &sizeof_complex,C_LOC(gntuju(1,1,ic,ig)),lmmaxapw*nufrmax,&
          ! &h_d_gntuju(idx),lmmaxapw*nufrmax, stream)

          !stat  = cublasSetMatrixAsync(lmmaxapw*nufrmax,1,sizeof_complex,&
          ! &C_LOC(b2Batch(1,1,idx)),lmmaxapw*nufrmax, h_d_b2(idx), &
          ! &lmmaxapw*nufrmax, stream)

          !stat = cudaMemset( h_d_b2(idx), 0, lmmaxapw*nufrmax*sizeof_complex)

          !!call zgemm('N','N',lmmaxapw*nufrmax,1,lmmaxapw*nufrmax,&
          !!  &zone,gntuju(1,1,ic,ig),lmmaxapw*nufrmax,b1,lmmaxapw*nufrmax,&
          !!  &zzero,b2,lmmaxapw*nufrmax)
          !!wftmp1((ias-1)*lmmaxapw*nufrmax+1:ias*lmmaxapw*nufrmax,ig)=b2(:)
        enddo !ias
      enddo !ig  

      ! Set all b2Batch to zero.
      stat = cudaMemset( h_d_b2(1), 0, b2Size)

      ! Transfer all of b1Batch at once.
      stat = cudaMemcpy(h_d_b1(1), C_LOC(b1Batch(1)), b1Size, cudaMemcpyHostToDevice)

      !!write(*,*) 'Finished calling SetMatrixAsync'
      !!flush(6)
      !do ig=1, batch_count
      !  cublasSetMatrixAsync(lmmaxapw*nufrmax,1,sizeof_complex,&
      !   &C_LOC(b1Batch(1,1,ig)),lmmaxapw*nufrmax, h_d_b1, &
      !   &lmmaxapw*nufrmax, stream)
      !  cublasSetMatrixAsync(lmmaxapw*nufrmax,lmmaxapw*nufrmax,sizeof_complex,&
      !   &C_LOC(gntujuBatch(1,1,ig)),lmmaxapw*nufrmax, h_d_gntuju, &
      !   &lmmaxapw*nufrmax, stream)
      !enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call timer_stop(3)
      call papi_timer_stop(pt_megqblh_mt)
! interstitial part
      call papi_timer_start(pt_megqblh_it)
      call timer_start(4)
      wfir1=zzero
      do ig1=1,ngknr1
        ifg=igfft(igkignr1(ig1))
        wfir1(ifg)=wfsvit1(ig1,ispn1,ist1)
      enddo
      call zfftifc(3,ngrid,1,wfir1)
      do ir=1,ngrtot
        wfir1(ir)=wfir1(ir)*cfunir(ir)
      enddo
      call zfftifc(3,ngrid,-1,wfir1)
      do ig=1,ngq(iq)
        do ig2=1,ngknr2
! G1=G2-G-Gkq
          ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,igqig(ig,iq))-ivg(:,igkq)
          ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
          wftmp1(lmmaxapw*nufrmax*natmtot+ig2,ig)=dconjg(wfir1(ifg))
        enddo
      enddo
      call timer_stop(4)      
      call papi_timer_stop(pt_megqblh_it)
    endif !l1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   call cudaStreamSynchronize(stream)

    stat = cublasZgemmBatched(handle,CUBLAS_OP_N,CUBLAS_OP_N,lmmaxapw*nufrmax,&
       &1,lmmaxapw*nufrmax,zone,d_gntuju,lmmaxapw*nufrmax,d_b1,&
       &lmmaxapw*nufrmax,zzero,d_b2,lmmaxapw*nufrmax,batch_count)
    
    ! Get the result back from the  GPU.
    stat = cudaMemcpy( C_LOC(b2Batch), h_d_b2(1), b2Size, cudaMemcpyDeviceToHost)

    call timer_start(5)
    n1=0
! collect right |ket> states into matrix wftmp2
    do while ((i+n1).le.nmegqblh(ikloc))
      if (bmegqblh(1,i+n1,ikloc).ne.bmegqblh(1,i,ikloc)) exit
      ist2=bmegqblh(2,i+n1,ikloc)
      n1=n1+1
      call memcopy(wfsvmt2(1,1,1,ispn2,ist2),wftmp2(1,n1),16*lmmaxapw*nufrmax*natmtot)
      call memcopy(wfsvit2(1,ispn2,ist2),wftmp2(lmmaxapw*nufrmax*natmtot+1,n1),16*ngknr2)
    enddo !while

    call cudaStreamSynchronize(stream)
      do ig=1,ngq(iq)
        do ias=1,natmtot
           idx = ( (ig-1)*natmtot + ias - 1 )*(lmmaxapw*nufrmax) + 1
           idx2 = ( (ig-1)*natmtot + ias )*(lmmaxapw*nufrmax)
           wftmp1((ias-1)*lmmaxapw*nufrmax+1:ias*lmmaxapw*nufrmax,ig)=b2Batch(idx:idx2)
        enddo
      enddo

! update several matrix elements by doing matrix*matrix operation
!  me(ib,ig)=wftmp2(ig2,ib)^{T}*wftmp1(ig2,ig)
    call zgemm('T','N',n1,ngq(iq),wfsize,zone,wftmp2,wfsize,wftmp1,wfsize,&
      &zone,megqblh(i,1,ikloc),nstsv*nstsv)
    i=i+n1
    call timer_stop(5)
!    if (wprocrank) then
!    write(*,*) '  Finished iteration of do while loop'
!    flush(6)
!    endif
  enddo !while
enddo !ispn
deallocate(wftmp1)
deallocate(wftmp2)
deallocate(wfir1)
!if (wprocrank) then
!write(*,*) 'Finished do while loop and deallocated wftmp1, wftmp2, wfir1'
!flush(6)
!endif

call cudaFree(h_d_b1(1))
call cudaFree(h_d_b2(1))
call cudaFree(h_d_gntuju(1))

call cudaFree(d_b1)
call cudaFree(d_b2)
call cudaFree(d_gntuju)

!if (wprocrank) then
!write(*,*) 'Finished calling cudaFree over the batch count of host to device ptr arrays'
!flush(6)
!endif
deallocate(b1Batch)
deallocate(b2Batch)
deallocate(gntujuBatch)

deallocate(h_d_b1)
deallocate(h_d_b2)
deallocate(h_d_gntuju)
!if (wprocrank) then
!write(*,*) 'Deallocated b1Batch and b2Batch'
!flush(6)
!endif
!deallocate(gntujuBatch)
call cublasDestroy(handle)
!if (wprocrank) then
!write(*,*) 'Finished call cublasDestroy on handle'
!flush(6)
!endif
stat =  cudaStreamDestroy(stream)
!if (wprocrank) then
!write(*,*) 'Finished cudaStreamDestroy(stream)'
!flush(6)
!endif
call papi_timer_stop(pt_megqblh)
return
end
