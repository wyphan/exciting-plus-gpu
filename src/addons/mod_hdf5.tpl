module mod_hdf5

interface hdf5_write
  module procedure hdf5_write_i4,hdf5_write_d,hdf5_write_z
end interface

interface hdf5_read
  module procedure hdf5_read_i4,hdf5_read_d,hdf5_read_z
end interface

public hdf5_initialize
public hdf5_finalize
public hdf5_create_file
public hdf5_create_group

!--gzip
public hdf5_check_gzip_support
public hdf5_calc_chunksz
public hdf5_gzwrite_array_d
public hdf5_gzwrite_array_z
public hdf5_gzread_array_d
public hdf5_gzread_array_z
!--gzip

contains

subroutine hdf5_initialize
#ifdef _HDF5_
use hdf5
implicit none
integer ierr
call h5open_f(ierr)
#endif
end subroutine

subroutine hdf5_finalize
#ifdef _HDF5_
use hdf5
implicit none
integer ierr
call h5close_f(ierr)
#endif
end subroutine

subroutine hdf5_create_file(fname)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
#ifdef _HDF5_
integer ierr
integer(hid_t) h5_root_id
call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_file): h5fcreate_f returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_file): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,'("  fname: ",A)')trim(fname)
stop
#endif
end subroutine

subroutine hdf5_create_group(fname,path,gname)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: gname
#ifdef _HDF5_
integer(hid_t) h5_root_id,h5_group_id,h5_new_group_id
integer ierr

call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5fopen_f returned ",I6)')ierr
  goto 10
endif
call h5gopen_f(h5_root_id,trim(path),h5_group_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5gopen_f returned ",I6)')ierr
  goto 10
endif
call h5gcreate_f(h5_group_id,trim(gname),h5_new_group_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5gcreate_f returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(h5_new_group_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5gclose_f for the new group returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(h5_group_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5gclose_f for the existing path returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path  : ",A)')trim(path)
write(*,'("  gname : ",A)')trim(gname)  
stop
#endif
end subroutine

!----------------------!
!     hdf5_write_i4    !
!----------------------!
subroutine hdf5_write_i4(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
integer, intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call hdf5_write_array_i4(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!---------------------!
!     hdf5_write_d    !
!---------------------!
subroutine hdf5_write_d(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
real(8), intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call hdf5_write_array_d(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!---------------------!
!     hdf5_write_z    !
!---------------------!
subroutine hdf5_write_z(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
complex(8), intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)+1
  allocate(dims_(ndims))
  dims_(1)=2
  dims_(2:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=2
endif
call hdf5_write_array_d(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!---------------------!
!     hdf5_read_i4    !
!---------------------!
subroutine hdf5_read_i4(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
integer, intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call hdf5_read_array_i4(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!--------------------!
!     hdf5_read_d    !
!--------------------!
subroutine hdf5_read_d(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
real(8), intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call hdf5_read_array_d(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!--------------------!
!     hdf5_read_z    !
!--------------------!
subroutine hdf5_read_z(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
complex(8), intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)+1
  allocate(dims_(ndims))
  dims_(1)=2
  dims_(2:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=2
endif
call hdf5_read_array_d(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!--gzip------------------------------------------------------------------------
SUBROUTINE hdf5_check_gzip_support( decode, encode )
#ifdef _HDF5_
  USE hdf5
#endif
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: decode, encode
#ifdef _HDF5_
  ! Internal variables
  LOGICAL :: gzavail
  INTEGER :: gzsupport, ierr

  ! First, check if gzip/deflate filter is available
  CALL h5zfilter_avail_f( H5Z_FILTER_DEFLATE_F, gzavail, ierr )
  IF( gzavail ) THEN
     ! Then, check for decode and encode support
     CALL h5zget_filter_info_f( H5Z_FILTER_DEFLATE_F, gzsupport, ierr )
     decode = ( 0 /= IAND( H5Z_FILTER_DECODE_ENABLED_F, gzsupport ) )
     encode = ( 0 /= IAND( H5Z_FILTER_ENCODE_ENABLED_F, gzsupport ) )
  ELSE
     WRITE(*,*) "(hdf5_check_gzip_support): GZIP support unavailable"
     decode = .FALSE.
     encode = .FALSE.
  END IF ! gzavail
#else
  decode = .FALSE.
  encode = .FALSE.
#endif
  RETURN
END SUBROUTINE hdf5_check_gzip_support

! The default chunk cache size is 1 MB.
! This function checks if that limit is exceeded.
INTEGER FUNCTION hdf5_calc_chunksz( type, ndims, chunk )
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN) :: type
  INTEGER, INTENT(IN) :: ndims
  INTEGER, DIMENSION(ndims), INTENT(IN) :: chunk
  ! Internal variables
  INTEGER, PARAMETER :: MiB = 1024**2
  INTEGER :: sz, bytes

  ! Parse the data type
  SELECT CASE(type)
  CASE( 'r', 'R' )
     sz = 4
  CASE( 'd', 'D' )
     sz = 8
  CASE( 'c', 'C' )
     sz = 8
  CASE( 'z', 'Z' )
     sz = 16
  CASE DEFAULT
     WRITE(*,*) 'Error(hdf5_calc_chunksz): Unknown data type ' // type
     hdf5_calc_chunksz = -1
     RETURN
  END SELECT

  ! Calculate the chunk size in bytes
  bytes = sz * PRODUCT(chunk)

  ! Warn if bigger than 1MB
  IF ( bytes > MiB ) THEN
     WRITE(*,*) 'Warning(hdf5_calc_chunksz): Chunk size is larger than 1 MiB'
     WRITE(*,*) 'Please set a smaller chunk, or increase the chunk cache size'
  END IF

  hdf5_calc_chunksz = bytes
  RETURN
END FUNCTION hdf5_calc_chunksz

! It's a mess to pass assumed-shape arrays without Fortran 2008 support
! Please just pass the first element of a, for now
SUBROUTINE hdf5_gzwrite_array_z( a, ndims, dims, chunk, level, fname, path, nm )
  USE ISO_C_BINDING
#ifdef _HDF5_
  USE hdf5
#endif
  IMPLICIT NONE
  COMPLEX(KIND((0.D0,1.D0))), TARGET, INTENT(IN) :: a
  INTEGER, INTENT(IN) :: ndims, level
  INTEGER, DIMENSION(ndims), INTENT(IN) :: dims, chunk
  CHARACTER(LEN=*), INTENT(IN) :: fname, path, nm
#ifdef _HDF5_
  ! Internal variables
  INTEGER(HID_T) :: h5root, h5group, h5space, h5set, h5prop
  INTEGER(HSIZE_T), DIMENSION(ndims+1) :: h5dims, h5chunk
  REAL(KIND(1.D0)), DIMENSION(:), POINTER :: buf
  INTEGER :: ierr, sz_d
  LOGICAL :: gzavail, dummy
  CHARACTER(LEN=100) :: errmsg

  ! Convert dims and chunk to HSIZE_T
  h5dims(1) = 2
  h5dims(2:ndims+1) = dims(:)
  h5chunk(1) = 2
  h5chunk(2:ndims+1) = chunk(:)

  ! Cast complex->void*->double
  sz_d = PRODUCT(h5dims)
  ALLOCATE( buf( sz_d ) )
  CALL C_F_POINTER( C_LOC(a), buf, (/sz_d/) )

  ! Check gzip encoding support
  CALL hdf5_check_gzip_support( dummy, gzavail )
  IF( .NOT. gzavail ) THEN
     WRITE(*,*) "Error(hdf5_gzwrite_array_z): GZIP compression unavailable"
     ! Fall back to writing uncompressed data
     CALL hdf5_write_array_d( buf, ndims+1, INT(h5dims), fname, path, nm )
     ! Clean up
     NULLIFY( buf )
     RETURN
  END IF

  ! Open file
  CALL h5fopen_f( TRIM(fname), H5F_ACC_RDWR_F, h5root, ierr )
  ! Access group
  CALL h5gopen_f( h5root, TRIM(path), h5group, ierr )
  ! Allocate dataspace
  CALL h5screate_simple_f( ndims+1, h5dims, h5space, ierr )
  ! Create dataset creation property list (DCPL)
  CALL h5pcreate_f( H5P_DATASET_CREATE_F, h5prop, ierr )
  ! Enable data chunking
  CALL h5pset_chunk_f( h5prop, ndims+1, h5chunk, ierr )
  ! Enable gzip/deflate compression
  CALL h5pset_deflate_f( h5prop, level, ierr )
  ! Create dataset
  CALL h5dcreate_f( h5group, TRIM(nm), H5T_NATIVE_DOUBLE, h5space, h5set, ierr,&
                    dcpl_id=h5prop )
  ! Write COMPLEX data through void* buffer
  CALL h5dwrite_f( h5set, H5T_NATIVE_DOUBLE, buf, h5dims, ierr )
  ! Deallocate hdf5 resources
  CALL h5sclose_f( h5space, ierr )
  CALL h5pclose_f( h5prop, ierr )
  CALL h5dclose_f( h5set, ierr )
  CALL h5gclose_f( h5group, ierr )
  CALL h5fclose_f( h5root, ierr )
  NULLIFY( buf )
#endif
  RETURN
END SUBROUTINE hdf5_gzwrite_array_z

! It's a mess to pass assumed-shape arrays without Fortran 2008 support
! Please just pass the first element of a, for now
SUBROUTINE hdf5_gzread_array_z( a, ndims, dims, fname, path, nm )
  USE ISO_C_BINDING
#ifdef _HDF5_
  USE hdf5
#endif
  IMPLICIT NONE
  COMPLEX(KIND((0.D0,1.D0))), INTENT(INOUT), TARGET :: a
  INTEGER, INTENT(IN) :: ndims
  INTEGER, DIMENSION(ndims), INTENT(IN) :: dims
  CHARACTER(LEN=*), INTENT(IN) :: fname, path, nm
#ifdef _HDF5_
  ! Internal variables
  INTEGER(HID_T) :: h5root, h5group, h5set, h5prop
  INTEGER(HSIZE_T), DIMENSION(ndims+1) :: h5dims
  REAL(KIND(1.D0)), DIMENSION(:), POINTER :: buf
  INTEGER :: ierr, sz_d
  LOGICAL :: gzavail, dummy
  CHARACTER(LEN=100) :: errmsg

  ! Convert dims to HSIZE_T
  h5dims(1) = 2
  h5dims(2:ndims+1) = dims(:)

  ! Cast complex->void*->double*
  sz_d = PRODUCT(h5dims)
  ALLOCATE( buf( sz_d ) )
  CALL C_F_POINTER( C_LOC(a), buf, (/sz_d/) )

  ! Check gzip decoding support
  CALL hdf5_check_gzip_support( gzavail, dummy )
  IF( .NOT. gzavail ) THEN
     WRITE(*,*) "Error(hdf5_gzread_array_z): GZIP decompression unavailable"
     STOP
  END IF

  ! Open file
  CALL h5fopen_f( TRIM(fname), H5F_ACC_RDONLY_F, h5root, ierr )
  ! Access group
  CALL h5gopen_f( h5root, TRIM(path), h5group, ierr )
  ! Access dataset
  CALL h5dopen_f( h5group, TRIM(nm), h5set, ierr )
  ! Access dataset creation property list (DCPL), which should be present
  ! if the data is compressed
  CALL h5dget_create_plist_f( h5set, h5prop, ierr )
  IF (ierr /= 0) THEN
     WRITE(*,*) "Error(hdf5_gzread_array_z): Cannot read DCPL. &
                &Are you sure "//TRIM(nm)//" is compressed?"
     ! Clean up and halt
     NULLIFY( buf )
     STOP
  END IF
  ! Read COMPLEX data through buffer
  CALL h5dread_f( h5set, H5T_NATIVE_DOUBLE, buf, h5dims, ierr )
  ! Deallocate hdf5 resources
  CALL h5pclose_f( h5prop, ierr )
  CALL h5dclose_f( h5set, ierr )
  CALL h5gclose_f( h5group, ierr )
  CALL h5fclose_f( h5root, ierr )
  NULLIFY( buf )
#endif
  RETURN
END SUBROUTINE hdf5_gzread_array_z

! It's a mess to pass assumed-shape arrays without Fortran 2008 support
! Please just pass the first element of a, for now
SUBROUTINE hdf5_gzwrite_array_d( a, ndims, dims, chunk, level, fname, path, nm )
  USE ISO_C_BINDING
#ifdef _HDF5_
  USE hdf5
#endif
  IMPLICIT NONE
  REAL(KIND(1.D0)), TARGET, INTENT(IN) :: a
  INTEGER, INTENT(IN) :: ndims, level
  INTEGER, DIMENSION(ndims), INTENT(IN) :: dims, chunk
  CHARACTER(LEN=*), INTENT(IN) :: fname, path, nm
#ifdef _HDF5_
  ! Internal variables
  INTEGER(HID_T) :: h5root, h5group, h5space, h5set, h5prop
  INTEGER(HSIZE_T), DIMENSION(ndims) :: h5dims, h5chunk
  REAL(KIND(1.D0)), DIMENSION(:), POINTER :: buf
  INTEGER :: ierr, sz_d
  LOGICAL :: gzavail, dummy
  CHARACTER(LEN=100) :: errmsg

  ! Convert dims and chunk to HSIZE_T
  h5dims(:) = dims(:)
  h5chunk(:) = chunk(:)

  ! Cast double->void*->double*
  sz_d = PRODUCT(dims)
  ALLOCATE( buf( sz_d ) )
  CALL C_F_POINTER( C_LOC(a), buf, (/sz_d/) )

  ! Check gzip encoding support
  CALL hdf5_check_gzip_support( dummy, gzavail )
  IF( .NOT. gzavail ) THEN
     WRITE(*,*) "Error(hdf5_gzwrite_array_d): GZIP compression unavailable"
     ! Fall back to writing uncompressed data
     CALL hdf5_write_array_d( buf, ndims, dims, fname, path, nm )
     ! Clean up
     NULLIFY( buf )
     RETURN
  END IF

  ! Open file
  CALL h5fopen_f( TRIM(fname), H5F_ACC_RDWR_F, h5root, ierr )
  ! Access group
  CALL h5gopen_f( h5root, TRIM(path), h5group, ierr )
  ! Allocate dataspace
  CALL h5screate_simple_f( ndims, h5dims, h5space, ierr )
  ! Create dataset creation property list (DCPL)
  CALL h5pcreate_f( H5P_DATASET_CREATE_F, h5prop, ierr )
  ! Enable data chunking
  CALL h5pset_chunk_f( h5prop, ndims, h5chunk, ierr )
  ! Enable gzip/deflate compression
  CALL h5pset_deflate_f( h5prop, level, ierr )
  ! Create dataset
  CALL h5dcreate_f( h5group, TRIM(nm), H5T_NATIVE_DOUBLE, h5space, h5set, ierr,&
                    dcpl_id=h5prop )
  ! Write REAL(dp) data through void* buffer
  CALL h5dwrite_f( h5set, H5T_NATIVE_DOUBLE, buf, h5dims, ierr )
  ! Deallocate hdf5 resources
  CALL h5sclose_f( h5space, ierr )
  CALL h5pclose_f( h5prop, ierr )
  CALL h5dclose_f( h5set, ierr )
  CALL h5gclose_f( h5group, ierr )
  CALL h5fclose_f( h5root, ierr )
  NULLIFY( buf )
#endif
  RETURN
END SUBROUTINE hdf5_gzwrite_array_d

! It's a mess to pass assumed-shape arrays without Fortran 2008 support
! Please just pass the first element of a, for now
SUBROUTINE hdf5_gzread_array_d( a, ndims, dims, fname, path, nm )
  USE ISO_C_BINDING
#ifdef _HDF5_
  USE hdf5
#endif
  IMPLICIT NONE
  REAL(KIND(1.D0)), INTENT(INOUT), TARGET :: a
  INTEGER, INTENT(IN) :: ndims
  INTEGER, DIMENSION(ndims), INTENT(IN) :: dims
  CHARACTER(LEN=*), INTENT(IN) :: fname, path, nm
#ifdef _HDF5_
  ! Internal variables
  INTEGER(HID_T) :: h5root, h5group, h5set, h5prop
  INTEGER(HSIZE_T), DIMENSION(ndims) :: h5dims
  REAL(KIND(1.D0)), DIMENSION(:), POINTER :: buf
  INTEGER :: ierr, sz_d
  LOGICAL :: gzavail, dummy
  CHARACTER(LEN=100) :: errmsg

  ! Convert dims to HSIZE_T
  h5dims(:) = dims(:)

  ! Cast double->void*->double*
  sz_d = PRODUCT(h5dims)
  ALLOCATE( buf( sz_d ) )
  CALL C_F_POINTER( C_LOC(a), buf, (/sz_d/) )

  ! Check gzip decoding support
  CALL hdf5_check_gzip_support( gzavail, dummy )
  IF( .NOT. gzavail ) THEN
     WRITE(*,*) "Error(hdf5_gzread_array_z): GZIP decompression unavailable"
     ! Clean up and halt
     NULLIFY( buf )
     STOP
  END IF

  ! Open file
  CALL h5fopen_f( TRIM(fname), H5F_ACC_RDONLY_F, h5root, ierr )
  ! Access group
  CALL h5gopen_f( h5root, TRIM(path), h5group, ierr )
  ! Access dataset
  CALL h5dopen_f( h5group, TRIM(nm), h5set, ierr )
  ! Access dataset creation property list (DCPL), which should be present
  ! if the data is compressed
  CALL h5dget_create_plist_f( h5set, h5prop, ierr )
  IF (ierr /= 0) THEN
     WRITE(*,*) "Error(hdf5_gzread_array_z): Cannot read DCPL. &
                &Are you sure "//TRIM(nm)//" is compressed?"
     ! Clean up and halt
     NULLIFY( buf )
     STOP
  END IF
  ! Read REAL(dp) data through buffer
  CALL h5dread_f( h5set, H5T_NATIVE_DOUBLE, buf, h5dims, ierr )
  ! Deallocate hdf5 resources
  CALL h5pclose_f( h5prop, ierr )
  CALL h5dclose_f( h5set, ierr )
  CALL h5gclose_f( h5group, ierr )
  CALL h5fclose_f( h5root, ierr )
  NULLIFY( buf )
#endif
  RETURN
END SUBROUTINE hdf5_gzread_array_d
!--gzip------------------------------------------------------------------------

end module

#ifdef _HDF5_

@python ftypes=["integer(4)","real(8)"]
@python fsuffixes=["_i4","_d"]
@python fhdf5types=["H5T_NATIVE_INTEGER","H5T_NATIVE_DOUBLE"]
@python ntypes=2

@template begin
@template variable fsuffix
@template variable ftype
@template variable fhdf5type
@python for i in range(ntypes): fsuffix=fsuffixes[i]; ftype=ftypes[i]; fhdf5type=fhdf5types[i];
subroutine hdf5_write_array#fsuffix(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
#ftype, intent(in) :: a(*)
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(hsize_t), dimension(ndims) :: h_dims
character*100 errmsg

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5fopen_f returned ",I6)')ierr
  goto 10
endif
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5screate_simple_f returned ",I6)')ierr
  goto 10
endif
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5gopen_f returned ",I6)')ierr
  goto 10
endif
call h5dcreate_f(group_id,trim(nm),#fhdf5type,dataspace_id,dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5dcreate_f returned ",I6)')ierr
  goto 10
endif 
call h5dwrite_f(dataset_id,#fhdf5type,a,h_dims,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5dwrite_f returned ",I6)')ierr
  goto 10
endif 
call h5dclose_f(dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5dclose_f returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5gclose_f returned ",I6)')ierr
  goto 10
endif
call h5sclose_f(dataspace_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5sclose_f returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array#fsuffix): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,'(A)')trim(errmsg)
write(*,'("  ndims : ",I4)')ndims
write(*,'("  dims  : ",10I4)')dims
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path  : ",A)')trim(path)
write(*,'("  nm    : ",A)')trim(nm)
stop
end subroutine

@template end

@template begin
@template variable fsuffix
@template variable ftype
@template variable fhdf5type
@python for i in range(ntypes): fsuffix=fsuffixes[i]; ftype=ftypes[i]; fhdf5type=fhdf5types[i];
subroutine hdf5_read_array#fsuffix(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
#ftype, intent(out) :: a(*)
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims
character*100 errmsg

do i=1,ndims
  h_dims(i)=dims(i)
enddo

call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array#fsuffix): h5fopen_f returned ",I6)')ierr
  goto 10
endif
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array#fsuffix): h5gopen_f returned ",I6)')ierr
  goto 10
endif
call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array#fsuffix): h5dopen_f returned ",I6)')ierr
  goto 10
endif
call h5dread_f(dataset_id,#fhdf5type,a,h_dims,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array#fsuffix): h5dread_f returned ",I6)')ierr
  goto 10
endif
call h5dclose_f(dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array#fsuffix): h5dclose_f returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array#fsuffix): h5gclose_f returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array#fsuffix): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,*)
write(*,'(A)')trim(errmsg)
write(*,'("  ndims : ",I4)')ndims
write(*,'("  dims : ",10I4)')dims
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path : ",A)')trim(path)
write(*,'("  nm : ",A)')trim(nm)
stop
end subroutine

@template end

#endif

