MODULE mod_prec

  IMPLICIT NONE

  ! Portable precision constants
  ! TODO: use SELECTED_REAL_KIND() for true portability
  INTEGER, PARAMETER :: ds = KIND(1.E0)
  INTEGER, PARAMETER :: dd = KIND(1.D0)
  INTEGER, PARAMETER :: dc = KIND((0.E0,1.E0))
  INTEGER, PARAMETER :: dz = KIND((0.D0,1.D0))

END MODULE mod_prec
