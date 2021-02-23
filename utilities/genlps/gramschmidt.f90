PROGRAM gramschmidt

  USE mod_axes3d, ONLY: dd, CROSS_PRODUCT
  IMPLICIT NONE

  ! Zero
  REAL(KIND=dd), parameter :: dzero = 0._dd

  ! Input variables
  REAL(KIND=dd), DIMENSION(3) :: v1, v2

  ! Output variable
  REAL(KIND=dd), DIMENSION(3,3) :: u

  ! Internal variables
  REAL(KIND=dd), DIMENSION(3,3) :: v
  INTEGER :: i, j

  ! Interface to DNRM2() from BLAS
  INTERFACE
     FUNCTION DNRM2( N, X, INCX )
       USE mod_axes3d, ONLY: dd
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N, INCX
       REAL(KIND=dd), DIMENSION(*), INTENT(IN) :: X
       REAL(KIND=dd) :: DNRM2
     END FUNCTION DNRM2
  END INTERFACE

  ! Read input from stdin
  READ(*,*) v1(:)
  READ(*,*) v2(:)

  ! Echo inputs
  WRITE(*,*) 'v1= ', v1(:)
  WRITE(*,*) 'v2= ', v2(:)

  ! Set first and second vectors
  v(:,1) = v1(:)
  v(:,2) = v2(:)

  ! Set third vector to v1 x v2
  v(:,3) = CROSS_PRODUCT( v1, v2 )

  ! Begin main loop
  DO i = 1, 3
     u(:,i) = v(:,i) / DNRM2( 3, v(:,i), 1 )

     ! Orthogonalization loop
     DO j = i+1, 3
        v(:,j) = v(:,j) - DOT_PRODUCT( u(:,i), v(:,j) ) * u(:,i)
     END DO ! j
  END DO ! i

  ! Write output
  DO i = 1, 3
     WRITE(*,*) ( u(i,j), j = 1, 3 )
  END DO

  STOP
END PROGRAM gramschmidt
