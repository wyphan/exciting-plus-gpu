MODULE mod_magma

  USE ISO_C_BINDING
  IMPLICIT NONE

  ! MAGMA queue
  TYPE(C_PTR) :: queue

  ! Number of GPUs
  INTEGER(C_INT) :: ngpus

!==============================================================================
! These are copied from the following MAGMA 2.5.3 source files:
! - fortran/magma2.F90
! - fortran/magma2_common.F90
! - fortran/magma2_zfortran.F90
!==============================================================================

  !! Parameter constants from magma_types.h
  integer(c_int), parameter ::   &
    MagmaFalse         = 0,    &
    MagmaTrue          = 1,    &
    
    MagmaRowMajor      = 101,  &
    MagmaColMajor      = 102,  &
    
    MagmaNoTrans       = 111,  &
    MagmaTrans         = 112,  &
    MagmaConjTrans     = 113,  &
    
    MagmaUpper         = 121,  &
    MagmaLower         = 122,  &
    MagmaGeneral       = 123,  &
    MagmaFull          = 123,  &  !! deprecated, use MagmaGeneral
    
    MagmaNonUnit       = 131,  &
    MagmaUnit          = 132,  &
    
    MagmaLeft          = 141,  &
    MagmaRight         = 142,  &
    MagmaBothSides     = 143

!! =====================================================================
!! Parameter constants
  real(c_float),             parameter :: sdummy = 0
  real(c_double),            parameter :: ddummy = 0
  complex(c_float_complex),  parameter :: cdummy = 0
  complex(c_double_complex), parameter :: zdummy = 0
  integer(c_int),            parameter :: idummy = 0
  type(c_ptr),               parameter :: ptr_dummy = c_null_ptr

!! Intel ifort chokes on c_sizeof here, so use extension sizeof
!! see https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/495001
  integer(c_size_t), parameter :: &
       sizeof_real      = sizeof(sdummy), &
       sizeof_double    = sizeof(ddummy), &
       sizeof_complex   = sizeof(cdummy), &
       sizeof_complex16 = sizeof(zdummy), &
       sizeof_int       = sizeof(idummy), &
       sizeof_ptr       = sizeof(ptr_dummy)

#ifdef _MAGMA_
!! =============================================================================
!! Fortran interfaces to C functions
  interface

!! -------------------------------------------------------------------------
!! initialize
     subroutine magma_init() &
          bind(C, name="magma_init")
       use iso_c_binding
     end subroutine magma_init

     subroutine magma_finalize() &
          bind(C, name="magma_finalize")
       use iso_c_binding
     end subroutine magma_finalize

!! -------------------------------------------------------------------------
!! version
     subroutine magma_version( major, minor, micro ) &
          bind(C, name="magma_version")
       use iso_c_binding
       integer(c_int), target :: major, minor, micro
     end subroutine magma_version

    subroutine magma_print_environment() &
         bind(C, name="magma_print_environment")
      use iso_c_binding
    end subroutine magma_print_environment

!! -------------------------------------------------------------------------
!! timing
    real(c_double) function magma_wtime() &
         bind(C, name="magma_wtime")
      use iso_c_binding
    end function magma_wtime

    real(c_double) function magma_sync_wtime( queue ) &
         bind(C, name="magma_sync_wtime")
      use iso_c_binding
      type(c_ptr), value :: queue
    end function magma_sync_wtime

!! -------------------------------------------------------------------------
!! device support
    integer(c_int) function magma_num_gpus() &
         bind(C, name="magma_num_gpus")
      use iso_c_binding
    end function magma_num_gpus

    integer(c_int) function magma_get_device_arch() &
         bind(C, name="magma_getdevice_arch")
      use iso_c_binding
    end function magma_get_device_arch

    subroutine magma_get_device( dev ) &
         bind(C, name="magma_getdevice")
      use iso_c_binding
      integer(c_int), target :: dev
    end subroutine magma_get_device

    subroutine magma_set_device( dev ) &
         bind(C, name="magma_setdevice")
      use iso_c_binding
      integer(c_int), value :: dev
    end subroutine magma_set_device

    integer(c_size_t) function magma_mem_size( queue ) &
         bind(C, name="magma_mem_size")
      use iso_c_binding
      type(c_ptr), value :: queue
    end function magma_mem_size

!! -------------------------------------------------------------------------
!! queue support
    subroutine magma_queue_create_internal( dev, queue_ptr, func, file, line ) &
         bind(C, name="magma_queue_create_internal")
      use iso_c_binding
      integer(c_int), value :: dev
      type(c_ptr), target :: queue_ptr  !! queue_t*
      character(c_char) :: func, file
      integer(c_int), value :: line
    end subroutine magma_queue_create_internal

    subroutine magma_queue_destroy_internal( queue, func, file, line ) &
         bind(C, name="magma_queue_destroy_internal")
      use iso_c_binding
      type(c_ptr), value :: queue  !! queue_t
      character(c_char) :: func, file
      integer(c_int), value :: line
    end subroutine magma_queue_destroy_internal

    subroutine magma_queue_sync_internal( queue, func, file, line ) &
         bind(C, name="magma_queue_sync_internal")
      use iso_c_binding
      type(c_ptr), value :: queue  !! queue_t
      character(c_char) :: func, file
      integer(c_int), value :: line
    end subroutine magma_queue_sync_internal

    integer(c_int) function magma_queue_get_device( queue ) &
         bind(C, name="magma_queue_get_device")
      use iso_c_binding
      type(c_ptr), value :: queue  !! queue_t
    end function magma_queue_get_device

!! -------------------------------------------------------------------------
!! offsets pointers -- 1D vectors with inc
!! see offset.c
    type(c_ptr) function magma_soffset_1d( ptr, inc, i ) &
         bind(C, name="magma_soffset_1d")
      use iso_c_binding
      type(c_ptr),    value :: ptr
      integer(c_int), value :: inc, i
    end function magma_soffset_1d

    type(c_ptr) function magma_doffset_1d( ptr, inc, i ) &
         bind(C, name="magma_doffset_1d")
      use iso_c_binding
      type(c_ptr),    value :: ptr
      integer(c_int), value :: inc, i
    end function magma_doffset_1d

    type(c_ptr) function magma_coffset_1d( ptr, inc, i ) &
         bind(C, name="magma_coffset_1d")
      use iso_c_binding
      type(c_ptr),    value :: ptr
      integer(c_int), value :: inc, i
    end function magma_coffset_1d

    type(c_ptr) function magma_zoffset_1d( ptr, inc, i ) &
         bind(C, name="magma_zoffset_1d")
      use iso_c_binding
      type(c_ptr),    value :: ptr
      integer(c_int), value :: inc, i
    end function magma_zoffset_1d

    type(c_ptr) function magma_ioffset_1d( ptr, inc, i ) &
         bind(C, name="magma_ioffset_1d")
      use iso_c_binding
      type(c_ptr),    value :: ptr
      integer(c_int), value :: inc, i
    end function magma_ioffset_1d

!! -------------------------------------------------------------------------
!! offsets pointers -- 2D matrices with lda
!! see offset.c
    type(c_ptr) function magma_soffset_2d( ptr, lda, i, j ) &
         bind(C, name="magma_soffset_2d")
      use iso_c_binding
      type(c_ptr),    value:: ptr
      integer(c_int), value :: lda, i, j
    end function magma_soffset_2d

    type(c_ptr) function magma_doffset_2d( ptr, lda, i, j ) &
         bind(C, name="magma_doffset_2d")
      use iso_c_binding
      type(c_ptr),    value:: ptr
      integer(c_int), value :: lda, i, j
    end function magma_doffset_2d

    type(c_ptr) function magma_coffset_2d( ptr, lda, i, j ) &
         bind(C, name="magma_coffset_2d")
      use iso_c_binding
      type(c_ptr),    value:: ptr
      integer(c_int), value :: lda, i, j
    end function magma_coffset_2d

    type(c_ptr) function magma_zoffset_2d( ptr, lda, i, j ) &
         bind(C, name="magma_zoffset_2d")
      use iso_c_binding
      type(c_ptr),    value:: ptr
      integer(c_int), value :: lda, i, j
    end function magma_zoffset_2d

    type(c_ptr) function magma_ioffset_2d( ptr, lda, i, j ) &
         bind(C, name="magma_ioffset_2d")
      use iso_c_binding
      type(c_ptr),    value:: ptr
      integer(c_int), value :: lda, i, j
    end function magma_ioffset_2d

!! -------------------------------------------------------------------------
!! magma_malloc (GPU memory)
    integer(c_int) function magma_malloc( ptr, bytes ) &
         bind(C, name="magma_malloc")
      use iso_c_binding
      type(c_ptr), target :: ptr  !! void**
      integer(c_size_t), value :: bytes
    end function magma_malloc

    integer(c_int) function magma_free_internal( ptr, func, file, line ) &
         bind(C, name="magma_free_internal")
      use iso_c_binding
      type(c_ptr), value :: ptr  !! void*
      character(c_char) :: func, file
      integer(c_int), value :: line
    end function magma_free_internal

!! -------------------------------------------------------------------------
!! magma_malloc_cpu (CPU main memory)
!! these are aligned to 32-byte boundary
    integer(c_int) function magma_malloc_cpu( ptr, bytes ) &
         bind(C, name="magma_malloc_cpu")
      use iso_c_binding
      type(c_ptr), target :: ptr  !! void**
      integer(c_size_t), value :: bytes
    end function magma_malloc_cpu

    integer(c_int) function magma_free_cpu( ptr ) &
         bind(C, name="magma_free_cpu")
      use iso_c_binding
      type(c_ptr), value :: ptr  !! void*
    end function magma_free_cpu

!! -------------------------------------------------------------------------
!! magma_malloc_pinned (pinned CPU main memory)
    integer(c_int) function magma_malloc_pinned( ptr, bytes ) &
         bind(C, name="magma_malloc_pinned")
      use iso_c_binding
      type(c_ptr), target :: ptr  !! void**
      integer(c_size_t), value :: bytes
    end function magma_malloc_pinned

    integer(c_int) function magma_free_pinned_internal( ptr, func, file, line ) &
         bind(C, name="magma_free_pinned_internal")
      use iso_c_binding
      type(c_ptr), value :: ptr  !! void*
      character(c_char), value :: func, file
      integer(c_int), value :: line
    end function magma_free_pinned_internal

!! -------------------------------------------------------------------------
!! set/get
    subroutine magma_setmatrix_internal( &
         m, n, elemsize, hA_src, lda, dB_dst, ldb, queue, func, file, line ) &
         bind(C, name="magma_setmatrix_internal")
      use iso_c_binding
      integer(c_int),    value  :: m, n, elemsize, lda, ldb
      type(c_ptr),       value  :: hA_src
      type(c_ptr),       value  :: dB_dst
      type(c_ptr),       value  :: queue
      character(c_char), value  :: func, file
      integer(c_int),    value  :: line
    end subroutine magma_setmatrix_internal

    subroutine magma_getmatrix_internal( &
         m, n, elemsize, dA_src, lda, hB_dst, ldb, queue, func, file, line ) &
    bind(C, name="magma_getmatrix_internal")
      use iso_c_binding
      integer(c_int),    value  :: m, n, elemsize, lda, ldb
      type(c_ptr),       value  :: dA_src
      type(c_ptr),       value  :: hB_dst
      type(c_ptr),       value  :: queue
      character(c_char), value  :: func, file
      integer(c_int),    value  :: line
    end subroutine magma_getmatrix_internal
    
    subroutine magma_setvector_internal( &
         n, elemsize, hx_src, incx, dy_dst, incy, queue, func, file, line ) &
         bind(C, name="magma_setvector_internal")
      use iso_c_binding
      integer(c_int),    value  :: n, elemsize, incx, incy
      type(c_ptr),       value  :: hx_src
      type(c_ptr),       value  :: dy_dst
      type(c_ptr),       value  :: queue
      character(c_char), value  :: func, file
      integer(c_int),    value  :: line
    end subroutine magma_setvector_internal

    subroutine magma_getvector_internal( &
         n, elemsize, dx_src, incx, hy_dst, incy, queue, func, file, line ) &
         bind(C, name="magma_getvector_internal")
      use iso_c_binding
      integer(c_int),    value  :: n, elemsize, incx, incy
      type(c_ptr),       value  :: dx_src
      type(c_ptr),       value  :: hy_dst
      type(c_ptr),       value  :: queue
      character(c_char), value  :: func, file
      integer(c_int),    value  :: line
    end subroutine magma_getvector_internal

!! -------------------------------------------------------------------------
!! CPU interfaces (matrix in CPU memory)
    subroutine magma_zgetrf( m, n, A, lda, ipiv, info ) &
         bind(C, name="magma_zgetrf")
      use iso_c_binding
      integer(c_int),            value  :: m, n, lda
      complex(c_double_complex), target :: A(lda,*)
      integer(c_int),            target :: ipiv(*)
      integer(c_int),            target :: info  !! int*
    end subroutine magma_zgetrf

    subroutine magma_zpotrf( uplo, n, A, lda, info ) &
         bind(C, name="magma_zpotrf")
      use iso_c_binding
      integer(c_int),            value  :: uplo
      integer(c_int),            value  :: n, lda
      complex(c_double_complex), target :: A(lda,*)
      integer(c_int),            target :: info  !! int*
    end subroutine magma_zpotrf

!! -------------------------------------------------------------------------
!! GPU interfaces (matrix in GPU memory)
    subroutine magma_zgetrf_gpu( m, n, dA, lda, ipiv, info ) &
         bind(C, name="magma_zgetrf_gpu")
      use iso_c_binding
      integer(c_int), value  :: m, n, lda
      type(c_ptr),    value  :: dA
      integer(c_int), target :: ipiv(*)
      integer(c_int), target :: info  !! int*
    end subroutine magma_zgetrf_gpu

    subroutine magma_zpotrf_gpu( uplo, n, dA, lda, info ) &
         bind(C, name="magma_zpotrf_gpu")
      use iso_c_binding
      integer(c_int), value  :: uplo, n, lda
      type(c_ptr),    value  :: dA
      integer(c_int), target :: info  !! int*
    end subroutine magma_zpotrf_gpu

!! -------------------------------------------------------------------------
!! batched GPU interfaces (all arrays in GPU memory)
    subroutine magma_zgetrf_batched( &
         m, n, dA_array, lda, ipiv_array, info_array, batchcount, queue ) &
         bind(C, name="magma_zgetrf_batched")
      use iso_c_binding
      integer(c_int), value  :: m, n, lda, batchcount
      type(c_ptr),    value  :: dA_array    !! double_complex**
      type(c_ptr),    value  :: ipiv_array  !! int**
      type(c_ptr),    value  :: info_array  !! int*
      type(c_ptr),    value  :: queue
    end subroutine magma_zgetrf_batched

!! -------------------------------------------------------------------------
!! BLAS (matrices in GPU memory)
    subroutine magma_zaxpy( &
         n, &
         alpha, dx, incx, &
                dy, incy, &
         queue ) &
         bind(C, name="magma_zaxpy")
      use iso_c_binding
      integer(c_int),             value :: n, incx, incy
      complex(c_double_complex),  value :: alpha
      type(c_ptr),                value :: dx, dy
      type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_zgemv( &
         transA, m, n, &
         alpha, dA, lda, &
                dx, incx, &
         beta,  dy, incy, &
         queue ) &
         bind(C, name="magma_zgemv")
      use iso_c_binding
      integer(c_int),             value :: transA, m, n, lda, incx, incy
      complex(c_double_complex),  value :: alpha, beta
      type(c_ptr),                value :: dA, dx, dy
      type(c_ptr),                value :: queue  !! queue_t
    end subroutine magma_zgemv

    subroutine magma_zgemm( &
         transA, transB, m, n, k, &
         alpha, dA, lda, &
                dB, ldb, &
         beta,  dC, ldc, &
         queue ) &
         bind(C, name="magma_zgemm")
      use iso_c_binding
      integer(c_int),             value :: transA, transB, m, n, k, lda, ldb, ldc
      complex(c_double_complex),  value :: alpha, beta
      type(c_ptr),                value :: dA, dB, dC
      type(c_ptr),                value :: queue  !! queue_t
    end subroutine magma_zgemm

    SUBROUTINE magma_zgemm_batched( &
         transA, transB, m, n, k, &
         alpha, dA_array, ldda, &
                dB_array, lddb, &
         beta,  dC_array, lddc, &
         batchCount, queue ) &
         bind(C, name="magma_zgemm_batched")
      use iso_c_binding
      integer(c_int),                     value :: transA, transB, m, n, k, &
                                                   ldda, lddb, lddc, &
                                                   batchCount
      complex(c_double_complex),          value :: alpha, beta
      type(c_ptr), DIMENSION(batchCount)        :: dA_array, dB_array, dC_array
      type(c_ptr),                        value :: queue
    END SUBROUTINE magma_zgemm_batched

end interface
#endif /* _MAGMA */

!==============================================================================
! For helper function magma_char
  INTEGER, PARAMETER :: &
       trans = 0, &
       uplo  = 1, &
       diag  = 2, &
       side  = 3

!! =============================================================================
!! Fortran routines & functions
contains

!! -------------------------------------------------------------------------
!! queue support
  subroutine magma_queue_create( dev, queue_ptr )
    use iso_c_binding
    integer(c_int), value :: dev
    type(c_ptr), target :: queue_ptr  !! queue_t*

#ifdef _MAGMA_
    call magma_queue_create_internal( &
         dev, queue_ptr, &
         "magma_queue_create" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_queue_create

  subroutine magma_queue_destroy( queue )
    use iso_c_binding
    type(c_ptr), value :: queue  !! queue_t

#ifdef _MAGMA_
    call magma_queue_destroy_internal( &
         queue, &
         "magma_queue_destroy" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_queue_destroy

  subroutine magma_queue_sync( queue )
    use iso_c_binding
    type(c_ptr), value :: queue  !! queue_t

#ifdef _MAGMA_
    call magma_queue_sync_internal( &
         queue, &
         "magma_queue_sync" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_queue_sync

!! -------------------------------------------------------------------------
!! malloc wrappers
  integer(c_int) function magma_imalloc( ptr, n )
    use iso_c_binding
    type(c_ptr),       target :: ptr  !! void**
    integer(c_size_t), value  :: n

#ifdef _MAGMA_
    magma_imalloc = magma_malloc( ptr, n*sizeof_int )
#endif /* _MAGMA_ */

  end function magma_imalloc

  integer(c_int) function magma_imalloc_cpu( ptr, n )
    use iso_c_binding
    type(c_ptr),       target :: ptr  !! void**
    integer(c_size_t), value  :: n

#ifdef _MAGMA_
    magma_imalloc_cpu = magma_malloc_cpu( ptr, n*sizeof_int )
#endif /* _MAGMA_ */

  end function magma_imalloc_cpu

  integer(c_int) function magma_imalloc_pinned( ptr, n )
    use iso_c_binding
    type(c_ptr),       target :: ptr  !! void**
    integer(c_size_t), value  :: n

#ifdef _MAGMA_
    magma_imalloc_pinned = magma_malloc_pinned( ptr, n*sizeof_int )
#endif /* _MAGMA_ */

  end function magma_imalloc_pinned

!! -------------------------------------------------------------------------
!! magma_free wrappers
  integer(c_int) function magma_free( ptr )
    type(c_ptr) :: ptr

#ifdef _MAGMA_
    magma_free = magma_free_internal( &
         ptr, &
         "magma_free" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end function magma_free

  integer(c_int) function magma_free_pinned( ptr )
    type(c_ptr) :: ptr

#ifdef _MAGMA_
    magma_free_pinned = magma_free_internal( &
         ptr, &
         "magma_free_pinned" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end function magma_free_pinned

!! -------------------------------------------------------------------------
!! set/get wrappers
  subroutine magma_setmatrix( &
       m, n, elemsize, hA_src, lda, dB_dst, ldb, queue )
    use iso_c_binding
    integer(c_int),    value  :: m, n, elemsize, lda, ldb
    type(c_ptr),       value  :: hA_src
    type(c_ptr),       value  :: dB_dst
    type(c_ptr),       value  :: queue

#ifdef _MAGMA_
    call magma_setmatrix_internal( &
         m, n, elemsize, hA_src, lda, dB_dst, ldb, queue, &
         "magma_setmatrix" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_setmatrix

  subroutine magma_getmatrix( &
       m, n, elemsize, dA_src, lda, hB_dst, ldb, queue )
    use iso_c_binding
    integer(c_int),    value  :: m, n, elemsize, lda, ldb
    type(c_ptr),       value  :: dA_src
    type(c_ptr),       value  :: hB_dst
    type(c_ptr),       value  :: queue

#ifdef _MAGMA_
    call magma_getmatrix_internal( &
         m, n, elemsize, dA_src, lda, hB_dst, ldb, queue, &
         "magma_getmatrix" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_getmatrix

  subroutine magma_setvector( &
       n, elemsize, hx_src, incx, dy_dst, incy, queue )
    use iso_c_binding
    integer(c_int),    value  :: n, elemsize, incx, incy
    type(c_ptr),       value  :: hx_src
    type(c_ptr),       value  :: dy_dst
    type(c_ptr),       value  :: queue

#ifdef _MAGMA_
    call magma_setvector_internal( &
         n, elemsize, hx_src, incx, dy_dst, incy, queue, &
         "magma_setvector" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_setvector

  subroutine magma_getvector( &
       n, elemsize, dx_src, incx, hy_dst, incy, queue )
    use iso_c_binding
    integer(c_int),    value  :: n, elemsize, incx, incy
    type(c_ptr),       value  :: dx_src
    type(c_ptr),       value  :: hy_dst
    type(c_ptr),       value  :: queue

#ifdef _MAGMA_
    call magma_getvector_internal( &
         n, elemsize, dx_src, incx, hy_dst, incy, queue, &
         "magma_getvector" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_getvector

!! -------------------------------------------------------------------------
!! set/get wrappers
!! matrices & vectors of integers
  subroutine magma_isetmatrix( &
       m, n, hA_src, lda, dB_dst, ldb, queue )
    use iso_c_binding
    integer(c_int), value  :: m, n, lda, ldb
    integer(c_int), target :: hA_src(lda,*)
    type(c_ptr),    value  :: dB_dst
    type(c_ptr),    value  :: queue

#ifdef _MAGMA_
    call magma_setmatrix_internal( &
         m, n, int(sizeof_int), c_loc(hA_src), lda, dB_dst, ldb, queue, &
         "magma_isetmatrix" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_isetmatrix

  subroutine magma_igetmatrix( &
       m, n, dA_src, lda, hB_dst, ldb, queue )
    use iso_c_binding
    integer(c_int), value  :: m, n, lda, ldb
    type(c_ptr),    value  :: dA_src
    integer(c_int), target :: hB_dst(ldb,*)
    type(c_ptr),    value  :: queue

#ifdef _MAGMA_
    call magma_getmatrix_internal( &
         m, n, int(sizeof_int), dA_src, lda, c_loc(hB_dst), ldb, queue, &
         "magma_igetmatrix" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_igetmatrix

  subroutine magma_isetvector( &
       n, hx_src, incx, dy_dst, incy, queue )
    use iso_c_binding
    integer(c_int), value  :: n, incx, incy
    integer(c_int), target :: hx_src(*)
    type(c_ptr),    value  :: dy_dst
    type(c_ptr),    value  :: queue

#ifdef _MAGMA_
    call magma_setvector_internal( &
         n, int(sizeof_int), c_loc(hx_src), incx, dy_dst, incy, queue, &
         "magma_isetvector" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_isetvector

  subroutine magma_igetvector( &
       n, dx_src, incx, hy_dst, incy, queue )
    use iso_c_binding
    integer(c_int), value  :: n, incx, incy
    type(c_ptr),    value  :: dx_src
    integer(c_int), target :: hy_dst(*)
    type(c_ptr),    value  :: queue

#ifdef _MAGMA_
    call magma_getvector_internal( &
         n, int(sizeof_int), dx_src, incx, c_loc(hy_dst), incy, queue, &
         "magma_igetvector" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_igetvector

!! -------------------------------------------------------------------------
!! set/get wrappers
!! matrices & vectors of c_ptr pointers
  subroutine magma_psetmatrix( &
       m, n, hA_src, lda, dB_dst, ldb, queue )
    use iso_c_binding
    integer(c_int), value  :: m, n, lda, ldb
    type(c_ptr),    target :: hA_src(lda,*)
    type(c_ptr),    value  :: dB_dst
    type(c_ptr),    value  :: queue

#ifdef _MAGMA_
    call magma_setmatrix_internal( &
         m, n, int(sizeof_ptr), c_loc(hA_src), lda, dB_dst, ldb, queue, &
         "magma_psetmatrix" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_psetmatrix

  subroutine magma_pgetmatrix( &
       m, n, dA_src, lda, hB_dst, ldb, queue )
    use iso_c_binding
    integer(c_int), value  :: m, n, lda, ldb
    type(c_ptr),    value  :: dA_src
    type(c_ptr),    target :: hB_dst(ldb,*)
    type(c_ptr),    value  :: queue

#ifdef _MAGMA_
    call magma_getmatrix_internal( &
         m, n, int(sizeof_ptr), dA_src, lda, c_loc(hB_dst), ldb, queue, &
         "magma_pgetmatrix" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_pgetmatrix

  subroutine magma_psetvector( &
       n, hx_src, incx, dy_dst, incy, queue )
    use iso_c_binding
    integer(c_int), value  :: n, incx, incy
    type(c_ptr),    target :: hx_src(*)
    type(c_ptr),    value  :: dy_dst
    type(c_ptr),    value  :: queue

#ifdef _MAGMA_
    call magma_setvector_internal( &
         n, int(sizeof_ptr), c_loc(hx_src), incx, dy_dst, incy, queue, &
         "magma_psetvector" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_psetvector

  subroutine magma_pgetvector( &
       n, dx_src, incx, hy_dst, incy, queue )
    use iso_c_binding
    integer(c_int), value  :: n, incx, incy
    type(c_ptr),    value  :: dx_src
    type(c_ptr),    target :: hy_dst(*)
    type(c_ptr),    value  :: queue

#ifdef _MAGMA_
    call magma_getvector_internal( &
         n, int(sizeof_ptr), dx_src, incx, c_loc(hy_dst), incy, queue, &
         "magma_pgetvector" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_pgetvector

!! -------------------------------------------------------------------------
!! malloc wrappers
  integer(c_int) function magma_zmalloc( ptr, n )
    use iso_c_binding
    type(c_ptr),       target :: ptr  !! void**
    integer(c_size_t), value  :: n

#ifdef _MAGMA_
    magma_zmalloc = magma_malloc( ptr, n*sizeof_complex16 )
#endif /* _MAGMA_ */

  end function magma_zmalloc

  integer(c_int) function magma_zmalloc_cpu( ptr, n )
    use iso_c_binding
    type(c_ptr),       target :: ptr  !! void**
    integer(c_size_t), value  :: n

#ifdef _MAGMA_
    magma_zmalloc_cpu = magma_malloc_cpu( ptr, n*sizeof_complex16 )
#endif /* _MAGMA_ */

  end function magma_zmalloc_cpu

  integer(c_int) function magma_zmalloc_pinned( ptr, n )
    use iso_c_binding
    type(c_ptr),       target :: ptr  !! void**
    integer(c_size_t), value  :: n

#ifdef _MAGMA_
    magma_zmalloc_pinned = magma_malloc_pinned( ptr, n*sizeof_complex16 )
#endif /* _MAGMA_ */

  end function magma_zmalloc_pinned

!! -------------------------------------------------------------------------
!! set/get wrappers
  subroutine magma_zsetmatrix( &
       m, n, hA_src, lda, dB_dst, ldb, queue )
    use iso_c_binding
    integer(c_int),            value  :: m, n, lda, ldb
    complex(c_double_complex), target :: hA_src(lda,*)
    type(c_ptr),               value  :: dB_dst
    type(c_ptr),               value  :: queue

#ifdef _MAGMA_
    call magma_setmatrix_internal( &
         m, n, int(sizeof_complex16), c_loc(hA_src), lda, dB_dst, ldb, queue, &
         "magma_zsetmatrix" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_zsetmatrix

  subroutine magma_zgetmatrix( &
       m, n, dA_src, lda, hB_dst, ldb, queue )
    use iso_c_binding
    integer(c_int),            value  :: m, n, lda, ldb
    type(c_ptr),               value  :: dA_src
    complex(c_double_complex), target :: hB_dst(ldb,*)
    type(c_ptr),               value  :: queue

#ifdef _MAGMA_
    call magma_getmatrix_internal( &
         m, n, int(sizeof_complex16), dA_src, lda, c_loc(hB_dst), ldb, queue, &
         "magma_zgetmatrix" // c_null_char, &
         __FILE__ // c_null_char, &
         __LINE__ )
#endif /* _MAGMA_ */

  end subroutine magma_zgetmatrix

!==============================================================================
! Initialize MAGMA and create queue
  SUBROUTINE magma_init_f

    USE mod_mpi_grid, ONLY: iproc
    IMPLICIT NONE

#ifdef _MAGMA_

    ! Local variables
    INTEGER :: devnum

    ! Only the master thread has access to MAGMA
    !$OMP MASTER

    ! Initialize MAGMA
    CALL magma_init

    ! Get number of GPU devices and assign each rank to a different device
    ngpus = magma_num_gpus()
    devnum = MOD(iproc, ngpus)

    ! Create MAGMA queue on device
    CALL magma_queue_create( devnum, queue )
    WRITE(*,*) 'Rank ',  iproc, ' created MAGMA queue on GPU ', devnum

    !$OMP END MASTER

#else

    IF( iproc == 0 ) WRITE(*,*) 'Please recompile with _MAGMA_ flag enabled'

#endif /* _MAGMA_ */

    RETURN
  END SUBROUTINE magma_init_f

!==============================================================================
! Destroy MAGMA queue and finalize MAGMA
  SUBROUTINE magma_finalize_f
    IMPLICIT NONE

#ifdef _MAGMA_

    ! Only the master thread has access to MAGMA
    !$OMP MASTER

    ! Destroy MAGMA queue
    CALL magma_queue_destroy( queue )

    ! Finalize MAGMA
    CALL magma_finalize

    !$OMP END MASTER

#endif /* _MAGMA_ */

    RETURN
  END SUBROUTINE magma_finalize_f

!==============================================================================
! BLAS-like interface to ZGEMM

  SUBROUTINE magma_zgemm_f( transA, transB, m, n, k, &
                            alpha, dA, lda, &
                                   dB, ldb, &
                            beta,  dC, ldc )
    USE ISO_C_BINDING
    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: transA, transB
    INTEGER(KIND=C_INT), VALUE :: m, n, k, lda, ldb, ldc
    COMPLEX(KIND=C_DOUBLE_COMPLEX), VALUE :: alpha, beta
#ifdef _OPENACC
    ! The device pointers will be extracted
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(lda,*), TARGET :: dA
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(ldb,*), TARGET :: dB
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(ldc,*), TARGET :: dC
#else
    ! These are already device pointers
    TYPE(C_PTR), INTENT(IN) :: dA, dB
    TYPE(C_PTR), INTENT(INOUT) :: dC
#endif /* _OPENACC */

#ifdef _MAGMA_

    ! Internal variables
    INTEGER :: ierr
    INTEGER(KIND=C_INT) :: op_a, op_b
    TYPE(C_PTR) :: dptr_a, dptr_b, dptr_c

    ! Only the master thread has access to MAGMA
    !$OMP MASTER

    ! Map transA and transB to enum using helper function
    op_a = magma_char( 'magma_zgemm_f', transA, trans )
    op_b = magma_char( 'magma_zgemm_f', transB, trans )

#ifdef _OPENACC

    ! Expose device pointers
    !$ACC HOST_DATA USE_DEVICE( dA, dB, dC )

    ! Extract device pointers
    dptr_a = C_LOC( dA )
    dptr_b = C_LOC( dB )
    dptr_c = C_LOC( dC )

    ! Call MAGMA with extracted device pointers
    CALL magma_zgemm( op_a, op_b, m, n, k, alpha, dptr_a, lda, dptr_b, ldb, beta, dptr_c, ldc, queue )

    !$ACC END HOST_DATA

#else

    ! Call MAGMA with device pointers passed directly
    CALL magma_zgemm( op_a, op_b, m, n, k, alpha, dA, lda, dB, ldb, beta, dC, ldc, queue )

#endif /* _OPENACC */

    !$OMP END MASTER

#endif /* _MAGMA_ */

    RETURN
  END SUBROUTINE magma_zgemm_f

!==============================================================================
! BLAS-like interface to batched ZGEMM

  SUBROUTINE magma_zgemm_batched_f( transA, transB, m, n, k, &
                                    alpha, dA, lda, &
                                           dB, ldb, &
                                    beta,  dC, ldc, &
                                    batchCount )
    USE ISO_C_BINDING
    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: transA, transB
    INTEGER(KIND=C_INT), VALUE :: m, n, k, lda, ldb, ldc, batchCount
    COMPLEX(KIND=C_DOUBLE_COMPLEX), VALUE :: alpha, beta
#ifdef _OPENACC
    ! The device pointers will be extracted
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:,:,:), TARGET :: dA
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:,:,:), TARGET :: dB
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:,:,:), TARGET :: dC
#else
    ! These are already device pointers
    TYPE(C_PTR), DIMENSION(batchCount), INTENT(IN) :: dA, dB
    TYPE(C_PTR), DIMENSION(batchCount), INTENT(INOUT) :: dC
#endif /* _OPENACC */

#ifdef _MAGMA_

    ! Internal variables
    INTEGER :: ierr, batch
    INTEGER(KIND=C_INT) :: op_a, op_b
    TYPE(C_PTR), DIMENSION(batchCount) :: dptr_a, dptr_b, dptr_c

    ! Only the master thread has access to MAGMA
    !$OMP MASTER

    ! Map transA and transB to enum using helper function
    op_a = magma_char( 'magma_zgemm_f', transA, trans )
    op_b = magma_char( 'magma_zgemm_f', transB, trans )

#ifdef _OPENACC

    ! Expose device pointers
    !$ACC HOST_DATA USE_DEVICE( dA, dB, dC )

    ! Extract device pointers
    DO batch = 1, batchCount
       dptr_a(batch) = C_LOC( dA( LBOUND(dA,1), LBOUND(dA,2), batch) )
       dptr_b(batch) = C_LOC( dB( LBOUND(dB,1), LBOUND(dB,2), batch) )
       dptr_c(batch) = C_LOC( dC( LBOUND(dC,1), LBOUND(dC,2), batch) )
    END DO

    ! Call MAGMA with extracted device pointers
    CALL magma_zgemm_batched( op_a, op_b, m, n, k, &
                              alpha, dptr_a, lda, &
                                     dptr_b, ldb, &
                               beta, dptr_c, ldc, &
                              batchCount, queue )

    !$ACC END HOST_DATA

#else

    ! Call MAGMA with device pointers passed directly
    CALL magma_zgemm_batched( op_a, op_b, m, n, k, &
                              alpha, dA, lda, &
                                     dB, ldb, &
                               beta, dC, ldc, &
                              batchCount, queue )

#endif /* _OPENACC */

    !$OMP END MASTER

#endif /* _MAGMA_ */

    RETURN
  END SUBROUTINE magma_zgemm_batched_f

!==============================================================================
! Helper function to translate chars to enums
!==============================================================================
  INTEGER(C_INT) FUNCTION magma_char( fname, char, type ) RESULT(num)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: fname
    CHARACTER(LEN=1), INTENT(IN) :: char
    INTEGER, INTENT(IN) :: type

    SELECT CASE( type )

    CASE( trans )
       SELECT CASE( char )
       CASE( 'N', 'n' )
          num = MagmaNoTrans
       CASE( 'T', 't' )
          num = MagmaTrans
       CASE( 'C', 'c' )
          num = MagmaConjTrans
       CASE DEFAULT
          WRITE(*,*) fname, ': unrecognized option for trans (N/T/C): ', char
       END SELECT

    CASE( uplo )
       SELECT CASE( char )
       CASE( 'L', 'l' )
          num = MagmaLower
       CASE( 'U', 'u' )
          num = MagmaUpper
       CASE( 'F', 'f', 'G', 'g' )
          num = MagmaGeneral
       CASE DEFAULT
          WRITE(*,*) fname, ': unrecognized option for uplo (L/U/G): ', char
       END SELECT

    CASE( diag )
       SELECT CASE( char )
       CASE( 'N', 'n' )
          num = MagmaNonUnit
       CASE( 'U', 'u' )
          num = MagmaUnit
       CASE DEFAULT
          WRITE(*,*) fname, ': unrecognized option for diag (N/U): ', char
       END SELECT

    CASE( side )
       SELECT CASE( char )
       CASE( 'L', 'l' )
          num = MagmaLeft
       CASE( 'R', 'r' )
          num = MagmaRight
       CASE( 'B', 'b' )
          num = MagmaBothSides
       CASE DEFAULT
          WRITE(*,*) fname, ': unrecognized option for side (L/R/B): ', char
       END SELECT

    CASE DEFAULT
       WRITE(*,*) fname, ': unrecognized option ', char
    END SELECT

    RETURN
  END FUNCTION magma_char

!==============================================================================

END MODULE mod_magma
