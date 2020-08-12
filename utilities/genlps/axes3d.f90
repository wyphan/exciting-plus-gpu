PROGRAM axes3d

  USE mod_axes3d
  IMPLICIT NONE

  ! Tolerance for iterative optimization
  REAL(KIND=dd), PARAMETER :: tol = 1.e-14_dd

  ! Input and output matrix
  ! On input, the 1st, 2nd, and 3rd column contain the x', y', and z' axes,
  ! respectively. Only the z' axis will remain unchanged.
  ! On output, the orthogonalized set of axes x", y" and z' to feed into genlps.
  REAL(KIND=dd), DIMENSION(3,3) :: mat, orth, newmat

  ! Internal variables
  REAL(KIND=dd), DIMENSION(3) :: x1, y1, z1, x2, y2, z2
  REAL(KIND=dd) :: alpha, beta, gamma
  TYPE(params) :: pair
  INTEGER :: lmax, i, j, iter

  ! Read the input file (filename hardcoded to be 'genlps.in')
  OPEN( UNIT=10, FILE='genlps.in', ACTION='READ', FORM='FORMATTED' )
  READ(10,*) lmax ! Won't be changed
  DO i = 1, 3
     READ(10,*) ( mat(i,j), j = 1, 3 )
  END DO
  CLOSE(10)

  ! Echo input matrix
  CALL matdisp( mat, 'Input matrix:' )

  ! Extract the three axes
  x1(:) = mat(:,1)
  y1(:) = mat(:,2)
  z1(:) = mat(:,3)

  ! Scale these vectors to unity
  x1(:) = x1(:) / SQRT( DOT_PRODUCT( x1, x1 ) )
  y1(:) = y1(:) / SQRT( DOT_PRODUCT( y1, y1 ) )
  z1(:) = z1(:) / SQRT( DOT_PRODUCT( z1, z1 ) )

  ! x' is expected to be already close to the desired value
  ! Rotate x' by beta to optimize angle between x' and z'
  pair%vecA(:) = z1(:)
  pair%vecB(:) = x1(:)
  beta = pi/180._dd
  CALL iterref( abstheta, beta, pair, tol, iter )
  y2(:) = CROSS_PRODUCT( z1, x1 )
  WRITE(*,*)
  WRITE(*,*) "Rotate x' by ", beta*180._dd/pi, "degrees about ", y2(:)
  WRITE(*,*) '  (after ', iter, ' iterations)'
  x2(:) = rotvec( x1, y2, beta )
  WRITE(*,*) "x'=", x1(:)
  WRITE(*,*) 'x"=', x2(:)
  x1(:) = x2(:) / SQRT( DOT_PRODUCT( x2, x2 ) )

  DO i = 1, maxiter

     ! Check orthogonality
     orth(:,:) = orthmat( x1, y1, z1 )
     CALL matdisp( orth, 'Orthogonality matrix:')

     ! First step: check cross product
     y2(:) = CROSS_PRODUCT( z1, x1 )
     y2(:) = y2(:) / SQRT( DOT_PRODUCT( y2, y2 ) )
     IF( DOT_PRODUCT( y1, y2 ) < DOT_PRODUCT( y2, y2 ) ) THEN
        ! Left-handed coordinate
        WRITE(*,*)
        WRITE(*,*) "Rotate y' by 180 degrees about ", z1(:)
        y2(:) = rotvec( y1, z1, pi )
        WRITE(*,*) "y'=", y1(:)
        WRITE(*,*) 'y"=', y2(:)
        y1(:) = y2(:) / SQRT( DOT_PRODUCT( y2, y2 ) )
     END IF ! Left-handed coordinate
     
     IF( ( ABS(orth(1,2)) < tol ) .AND. ( ABS(orth(1,3)) < tol ) .AND. &
         ( ABS(orth(2,3)) < tol ) ) THEN
        WRITE(*,*)
        WRITE(*,*) 'Axes optimization complete after ', i, ' iterations'
        EXIT
     END IF

     ! Second step: rotate y' by gamma to optimize angle between y' and z'
     pair%vecA(:) = z1(:)
     pair%vecB(:) = y1(:)
     gamma = pi/180._dd
     CALL iterref( abstheta, gamma, pair, tol, iter )
     x2(:) = CROSS_PRODUCT( y1, z1 )
     WRITE(*,*)
     WRITE(*,*) "Rotate y' by ", gamma*180._dd/pi, "degrees about ", x2(:)
     WRITE(*,*) '  (after ', iter, ' iterations)'
     y2(:) = rotvec( y1, x2, gamma )
     WRITE(*,*) "y'=", y1(:)
     WRITE(*,*) 'y"=', y2(:)
     y1(:) = y2(:) / SQRT( DOT_PRODUCT( y2, y2 ) )

     ! Check orthogonality
     orth(:,:) = orthmat( x1, y1, z1 )
     CALL matdisp( orth, 'Orthogonality matrix:')

     ! Third step: rotate y' by gamma to optimize angle between x' and y'
     pair%vecA(:) = x1(:)
     pair%vecB(:) = y1(:)
     alpha = pi/180._dd
     CALL iterref( abstheta, alpha, pair, tol, iter )
     z2(:) = CROSS_PRODUCT( x1, y1 )
     WRITE(*,*)
     WRITE(*,*) "Rotate y' by ", alpha*180._dd/pi, "degrees about ", z2(:)
     WRITE(*,*) '  (after ', iter, ' iterations)'
     y2(:) = rotvec( y1, z2, alpha )
     WRITE(*,*) "y'=", y1(:)
     WRITE(*,*) 'y"=', y2(:)
     y1(:) = y2(:) / SQRT( DOT_PRODUCT( y2, y2 ) )

     ! Print warning message
     IF( i == maxiter ) THEN
        WRITE(*,*) 'Warning: reached max number of iterations', maxiter
     END IF

  END DO ! i

  ! Form new matrix
  newmat(:,1) = x1(:)
  newmat(:,2) = y1(:)
  newmat(:,3) = z1(:)

  ! Write output to stdout
  CALL matdisp( newmat, 'New matrix:' )

  ! Write output to disk
  OPEN( UNIT=10, FILE='genlps.in', ACTION='WRITE', FORM='FORMATTED' )
  WRITE(10,*) lmax
  DO i = 1, 3
     WRITE(10,*) ( newmat(i,j), j = 1, 3 )
  END DO
  WRITE(10,*)
  WRITE(10,*) '! Original matrix:'
  DO i = 1, 3
     WRITE(10,*) ( mat(i,j), j = 1, 3 )
  END DO  
  CLOSE(10)

  STOP
END PROGRAM axes3d
