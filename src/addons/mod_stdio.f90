MODULE mod_stdio

#ifdef F2003
  ! Fetch units for stdin, stdout, and stderr
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: stdin  => input_unit, &
                                           stdout => output_unit, &
                                           stderr => error_unit
  IMPLICIT NONE
#else
  IMPLICIT NONE

  ! Define units for stdin, stdout, and stderr
  INTEGER, PARAMETER :: stdin  = 5
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: stderr = 0
#endif /* F2003 */

END MODULE mod_stdio
