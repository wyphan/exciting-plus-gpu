MODULE mod_axes3d

  ! Precision constants (from mod_prec)
  integer, parameter :: dd = selected_real_kind(14,100)  ! FP64

  ! Pi (from modmain)
  real(KIND=dd), parameter :: pi=3.1415926535897932385_dd

  ! Zero
  REAL(KIND=dd), parameter :: dzero = 0._dd

  ! Maximum number of iterations
  INTEGER, PARAMETER :: maxiter = 500

  TYPE params
     REAL(KIND=dd), DIMENSION(3) :: vecA, vecB
  END TYPE params

CONTAINS

!-------------------------------------------------------------------------------
! Calculates the cross product
! _   _   _
! C = A x B

  FUNCTION CROSS_PRODUCT( vecA, vecB )
    IMPLICIT NONE

    ! Arguments
    REAL(KIND=dd), DIMENSION(3), INTENT(IN) :: vecA, vecB
    REAL(KIND=dd), DIMENSION(3) :: CROSS_PRODUCT

    CROSS_PRODUCT(1) = vecA(2)*vecB(3) - vecA(3)*vecB(2)
    CROSS_PRODUCT(2) = vecA(3)*vecB(1) - vecA(1)*vecB(3)
    CROSS_PRODUCT(3) = vecA(1)*vecB(2) - vecA(2)*vecB(1)

    RETURN
  END FUNCTION CROSS_PRODUCT

!-------------------------------------------------------------------------------
! Calculates the angle between two vectors theta using the following formula:
!                   _   _        _   _  
! theta = arccos( ( A . B ) / ( |A| |B| ) )

  REAL(KIND=dd) FUNCTION vecang( vecA, vecB )
    IMPLICIT NONE

    ! Arguments
    REAL(KIND=dd), DIMENSION(3), INTENT(IN) :: vecA, vecB

    ! Internal variables
    REAL(KIND=dd) :: magA, magB

    ! Calculate magnitude of each vector
    magA = SQRT( DOT_PRODUCT( vecA, vecA ) )
    magB = SQRT( DOT_PRODUCT( vecB, vecB ) )
    
    vecang = ACOS( DOT_PRODUCT( vecA, vecB ) / ( magA * magB ) )

    RETURN
  END FUNCTION vecang

!-------------------------------------------------------------------------------
! Rotates a vector about a given axis n and a given angle phi (in radians,
! going counterclockwise) according to the following formula, taken from
! https://mathworld.wolfram.com/RotationFormula.html
! _    _            ^   ^   _                        _   ^
! r' = r cos(phi) + n ( n . r ) ( 1 - cos(phi) ) + ( r x n ) sin(phi)

  FUNCTION rotvec( vec, axis, phi )
    IMPLICIT NONE

    ! Arguments
    REAL(KIND=dd), DIMENSION(3), INTENT(IN) :: vec, axis
    REAL(KIND=dd), INTENT(IN) :: phi
    REAL(KIND=dd), DIMENSION(3) :: rotvec

    ! Internal variables
    REAL(KIND=dd), DIMENSION(3) :: uvec
    REAL(KIND=dd) :: cf, sf

    ! Calculate unit vector of axis
    uvec(:) = axis(:) / SQRT( DOT_PRODUCT( axis, axis ) )

    ! Calculate sin(phi) and cos(phi)
    cf = COS(phi)
    sf = SIN(phi)

    ! Finally, plug into the formula
    rotvec(:) = cf * vec(:) &
                + ( 1._dd - cf ) * DOT_PRODUCT( uvec, vec ) * uvec(:) &
                + sf * CROSS_PRODUCT( vec, uvec )

    RETURN
  END FUNCTION rotvec

!-------------------------------------------------------------------------------
! False position method
! - fx is the function to optimize to 0 (to find the root)
! - x is the variable to adjust. On input, the initial guess, and on output,
!     the optimized value that causes fx to converge
! - tol is the tolerance
! - iter is the number of iterations needed to converge within tolerance

  SUBROUTINE iterref( fx, x, param, tol, iter )
    IMPLICIT NONE

    ! Arguments
    REAL(KIND=dd), EXTERNAL :: fx, dfx
    REAL(KIND=dd), INTENT(INOUT) :: x
    TYPE(params), INTENT(IN) :: param
    REAL(KIND=dd), INTENT(IN) :: tol
    INTEGER, INTENT(OUT) :: iter

    ! Internal variables
    REAL(KIND=dd) :: x0, x1, y0, y1, dx
    INTEGER :: i, sgn

    ! Increment starts at 1 degree
    dx = pi/180._dd

    DO i = 1, maxiter

       ! Calculate old values
       x0 = x
       y0 = fx( x0, param )

!--DEBUG
!       WRITE(*,*) 'i=', i, 'x0=', x0, 'y0=', y0
!--DEBUG
         
       IF( ABS(y0) < tol ) THEN
          ! Tolerance achieved, we're done!
          iter = i
          EXIT
       END IF

       ! Calculate new values
       x1 = x0 + dx
       y1 = fx( x1, param )
       
       ! Use false position method
       x = ( x0 * y1 - x1 * y0 ) / ( y1 - y0 )

       ! Halve dx once per iteration
       dx = dx/2._dd

    END DO ! i

    RETURN
  END SUBROUTINE iterref

!-------------------------------------------------------------------------------
! The function to optimize by iterative refinement: the absolute value of the
! difference between pi/2 and theta, the angle between the two vectors (passed
! through vecs).
! Phi is the angle to rotate vecB, in case in needs to be adjusted.  

  REAL(KIND=dd) FUNCTION abstheta( phi, vecs )

    ! Arguments
    REAL(KIND=dd), INTENT(IN) :: phi
    TYPE(params), INTENT(IN) :: vecs

    ! Internal variables
    REAL(KIND=dd), DIMENSION(3) :: uvec, vecB1

    !                           _   _  
    ! Calculate rotation axis ( A x B )
    uvec(:) = CROSS_PRODUCT( vecs%vecA, vecs%vecB )
    !        _                   _   _            _
    ! Rotate B by phi about the (A x B) axis into B'
    vecB1(:) = rotvec( vecs%vecB, uvec, phi )
    !                       _     _
    ! Compute angle between A and B', and subtract it from pi/2
    abstheta = pi/2._dd - vecang( vecs%vecA, vecB1 )

    RETURN
  END FUNCTION abstheta

!-------------------------------------------------------------------------------
! Helper function to calculate orthogonality matrix to stdout

  FUNCTION orthmat( x, y, z ) RESULT( mat )
    IMPLICIT NONE

    ! Arguments
    REAL(KIND=dd), DIMENSION(3), INTENT(IN) :: x, y, z
    REAL(KIND=dd), DIMENSION(3,3) :: mat

    ! Internal variables
    REAL(KIND=dd) :: xdotx, ydoty, zdotz, xdoty, xdotz, ydotz

    ! Calculate the dot products for the six different combinations
    xdotx = DOT_PRODUCT( x, x )
    ydoty = DOT_PRODUCT( y, y )
    zdotz = DOT_PRODUCT( z, z )
    xdoty = DOT_PRODUCT( x, y )
    xdotz = DOT_PRODUCT( x, z )
    ydotz = DOT_PRODUCT( y, z )

    mat(:,:) = RESHAPE( (/ xdotx, xdoty, xdotz, xdoty, ydoty, ydotz, &
                           xdotz, ydotz, zdotz /), (/ 3, 3 /) )

    RETURN
  END FUNCTION orthmat

!-------------------------------------------------------------------------------
! Helper function to display matrix to stdout

  SUBROUTINE matdisp( mat, msg )
    IMPLICIT NONE

    ! Arguments
    REAL(KIND=dd), DIMENSION(:,:), INTENT(IN) :: mat
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: msg

    ! Internal variables
    INTEGER :: i, j

    IF( PRESENT(msg) ) THEN
       WRITE(*,*)
       WRITE(*,*) msg
    END IF

    DO i = 1, SIZE(mat,1)
       WRITE(*,*) ( mat(i,j), j = 1, SIZE(mat,2) )
    END DO

    RETURN
  END SUBROUTINE matdisp

!-------------------------------------------------------------------------------
END MODULE mod_axes3d
