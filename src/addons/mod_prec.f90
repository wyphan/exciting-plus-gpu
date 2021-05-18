MODULE mod_prec

#ifdef F2008
  USE ISO_FORTRAN_ENV
#endif

  IMPLICIT NONE

  ! Portable precision constants
  INTEGER, PARAMETER :: dl = SELECTED_INT_KIND(18)       ! INT64
  integer, parameter :: ds = selected_real_kind(6,30)    ! FP32
  integer, parameter :: dd = selected_real_kind(14,100)  ! FP64
  INTEGER, PARAMETER :: dc = KIND( (0.0_ds, 1.0_ds) )
  INTEGER, PARAMETER :: dz = KIND( (0.0_dd, 1.0_dd) )

  ! Test vars
  INTEGER(KIND=dl), PRIVATE :: dummy_l(1)
  REAL(KIND=ds),    PRIVATE :: dummy_s(1)
  REAL(KIND=dd),    PRIVATE :: dummy_d(1)
  COMPLEX(KIND=dc), PRIVATE :: dummy_c(1)
  COMPLEX(KIND=dz), PRIVATE :: dummy_z(1)

#if defined(__GFORTRAN__) || defined(__INTEL_COMPILER) || defined(__IBMC__)

  ! GCC, Intel, IBM include SIZEOF() as an extension
  CHARACTER(LEN=*), PARAMETER :: sz_method = 'intrinsic'
  INTEGER, PARAMETER :: sz_l = SIZEOF(dummy_l(1))
  INTEGER, PARAMETER :: sz_s = SIZEOF(dummy_s(1))
  INTEGER, PARAMETER :: sz_d = SIZEOF(dummy_d(1))
  INTEGER, PARAMETER :: sz_c = SIZEOF(dummy_c(1))
  INTEGER, PARAMETER :: sz_z = SIZEOF(dummy_z(1))

#elif defined(F2008)
  ! TODO: Find the correct universal preprocessor macro, if it exists

  ! Use F2008 STORAGE_SIZE() intrinsic from ISO_FORTRAN_ENV module
  CHARACTER(LEN=*), PARAMETER :: sz_method = 'iso_fortran_env'
  INTEGER, PARAMETER :: sz_l = STORAGE_SIZE(dummy_l) / CHARACTER_STORAGE_SIZE
  INTEGER, PARAMETER :: sz_s = STORAGE_SIZE(dummy_s) / CHARACTER_STORAGE_SIZE
  INTEGER, PARAMETER :: sz_d = STORAGE_SIZE(dummy_d) / CHARACTER_STORAGE_SIZE
  INTEGER, PARAMETER :: sz_c = STORAGE_SIZE(dummy_c) / CHARACTER_STORAGE_SIZE
  INTEGER, PARAMETER :: sz_z = STORAGE_SIZE(dummy_z) / CHARACTER_STORAGE_SIZE

#else

  ! Commonly used numbers in bytes
  CHARACTER(LEN=*), PARAMETER :: sz_method = 'constant'
  INTEGER, PARAMETER :: sz_l = 8
  INTEGER, PARAMETER :: sz_s = 4
  INTEGER, PARAMETER :: sz_d = 8
  INTEGER, PARAMETER :: sz_c = 8
  INTEGER, PARAMETER :: sz_z = 16

#endif /* method */

END MODULE mod_prec
