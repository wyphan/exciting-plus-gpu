MODULE mod_fft_acc
! Implements the Fast Fourier Transform algorithm from
! W.M. Gentleman and G. Sande,
! "Fast Fourier Transforms -- for Fun and Profit"
! Proceedings of the 1966 AFIPS Fall Joint Computer Conference
! Pages 563-578
! doi:10.1145/1464291.1464352
!===============================================================================

  USE mod_prec, ONLY: dd, dz
  USE modmain, ONLY: zzero, zhalf, zone, zi, pi, twopi, sqtwo
  IMPLICIT NONE

  ! Mathematical constants
  ! For size 3, 6, 9, 12
  REAL(KIND=dd), PARAMETER :: cos10 = 0.984807753012208_dd
  REAL(KIND=dd), PARAMETER :: cos20 = 0.939692620785908_dd
  REAL(KIND=dd), PARAMETER :: cos30 = 0.866025403784439_dd
  REAL(KIND=dd), PARAMETER :: cos40 = 0.766044443118978_dd
  REAL(KIND=dd), PARAMETER :: sin10 = 0.173648177666930_dd
  REAL(KIND=dd), PARAMETER :: sin20 = 0.342020143325669_dd
  REAL(KIND=dd), PARAMETER :: sin30 = 0.5_dd
  REAL(KIND=dd), PARAMETER :: sin40 = 0.642787609686539_dd
  COMPLEX(KIND=dz), PARAMETER :: ei30  = CMPLX(  cos30,  sin30, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei40  = CMPLX(  cos40,  sin40, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei60  = CMPLX(  sin30,  cos30, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei80  = CMPLX(  sin10,  cos10, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei120 = CMPLX( -cos60,  sin60, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei160 = CMPLX( -cos20,  sin20, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei240 = CMPLX( -sin30, -cos30, KIND=dz )
  ! For size 4, 8, 16
  REAL(KIND=dd), PARAMETER :: cos23 = 0.923879532511287_dd ! cos(pi/8)
  REAL(KIND=dd), PARAMETER :: cos45 = sqtwo/2.0_dd
  REAL(KIND=dd), PARAMETER :: sin23 = 0.382683432365090_dd ! sin(pi/8)
  COMPLEX(KIND=dz), PARAMETER :: ei23  = CMPLX(  cos23,  sin23, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei45  = CMPLX(  cos45,  cos45, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei68  = CMPLX(  sin23,  cos23, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei135 = CMPLX( -cos45,  cos45, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei203 = CMPLX( -cos23, -sin23, KIND=dz )
  ! For size 5, 10, 15, 20, 25
  REAL(KIND=dd), PARAMETER :: cos4  = 0.998026728428272_dd ! cos(pi/50)
  REAL(KIND=dd), PARAMETER :: cos6  = 0.994521895368273_dd ! cos(pi/30)
  REAL(KIND=dd), PARAMETER :: cos7  = 0.992114701314478_dd ! cos(pi/25)
  REAL(KIND=dd), PARAMETER :: cos12 = 0.978147600733806_dd ! cos(pi/15)
  REAL(KIND=dd), PARAMETER :: cos14 = 0.968583161128631_dd ! cos(2pi/25)
  REAL(KIND=dd), PARAMETER :: cos18 = 0.951056516295154_dd ! cos(pi/10)
  REAL(KIND=dd), PARAMETER :: cos24 = 0.913545457642601_dd ! cos(2pi/15)
  REAL(KIND=dd), PARAMETER :: cos25 = 0.904827052466020_dd ! cos(7pi/50)
  REAL(KIND=dd), PARAMETER :: cos29 = 0.876306680043864_dd ! cos(4pi/25)
  REAL(KIND=dd), PARAMETER :: cos32 = 0.844327925502015_dd ! cos(9pi/50)
  REAL(KIND=dd), PARAMETER :: cos36 = 0.809016994374947_dd ! cos(pi/5)
  REAL(KIND=dd), PARAMETER :: cos40 = 0.770513242775789_dd ! sin(11pi/50)
  REAL(KIND=dd), PARAMETER :: cos43 = 0.728968627421412_dd ! cos(6pi/25)
  REAL(KIND=dd), PARAMETER :: sin4  = 0.0627905195293134_dd ! sin(pi/50)
  REAL(KIND=dd), PARAMETER :: sin6  = 0.104528463267653_dd ! sin(pi/30)
  REAL(KIND=dd), PARAMETER :: sin7  = 0.125333233564304_dd ! sin(pi/25)
  REAL(KIND=dd), PARAMETER :: sin12 = 0.207911690817759_dd ! sin(pi/15)
  REAL(KIND=dd), PARAMETER :: sin14 = 0.248689887164855_dd ! sin(2pi/25)
  REAL(KIND=dd), PARAMETER :: sin18 = 0.309016994374947_dd ! sin(pi/10)
  REAL(KIND=dd), PARAMETER :: sin24 = 0.406736643075800_dd ! sin(2pi/15)
  REAL(KIND=dd), PARAMETER :: sin25 = 0.425779291565073_dd ! sin(7pi/50)
  REAL(KIND=dd), PARAMETER :: sin29 = 0.481753674101715_dd ! sin(4pi/25)
  REAL(KIND=dd), PARAMETER :: sin32 = 0.535826794978997_dd ! sin(9pi/50)
  REAL(KIND=dd), PARAMETER :: sin36 = 0.587785252292473_dd ! sin(pi/5)
  REAL(KIND=dd), PARAMETER :: sin40 = 0.637423989748690_dd ! sin(11pi/50)
  REAL(KIND=dd), PARAMETER :: sin43 = 0.684547105928689_dd ! sin(6pi/25)
  COMPLEX(KIND=dz), PARAMETER :: ei14  = CMPLX(  cos14,  sin14, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei18  = CMPLX(  cos18,  sin18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei24  = CMPLX(  cos24,  sin24, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei29  = CMPLX(  cos29,  sin29, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei36  = CMPLX(  cos36,  sin36, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei43  = CMPLX(  cos43,  sin43, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei54  = CMPLX(  sin36,  cos36, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei58  = CMPLX(  sin32,  cos32, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei72  = CMPLX(  sin18,  cos18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei86  = CMPLX(  sin4,   cos4,  KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei96  = CMPLX( -sin6,   cos6,  KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei108 = CMPLX( -sin18,  cos18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei115 = CMPLX( -sin25,  cos25, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei130 = CMPLX( -sin40,  cos40, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei144 = CMPLX( -cos36,  sin36, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei162 = CMPLX( -cos18,  sin18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei173 = CMPLX( -cos7,   sin7,  KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei192 = CMPLX( -cos12, -sin12, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei216 = CMPLX( -cos36, -sin36, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei230 = CMPLX( -cos40, -sin40, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei288 = CMPLX(  sin18, -cos18, KIND=dz )

  COMPLEX(KIND=dz), DIMENSION(0:1,0:1) :: twiddle2
  DATA twiddle2 / zone, zone, &
                  zone, -zone /
  
  COMPLEX(KIND=dz), DIMENSION(0:2,0:2) :: twiddle3
  DATA twiddle3 / zone, zone,  zone,  &
                  zone, ei120, ei240, &
                  zone, ei240, ei120 /

  COMPLEX(KIND=dz), DIMENSION(0:3,0:3) :: twiddle4
  DATA twiddle4 / zone, zone, zone, zone, &
                  zone, zi,  -zone,-zi,   &
                  zone,-zone, zone,-zone, &
                  zone,-zi,  -zone, zi    /

  COMPLEX(KIND=dz), DIMENSION(0:4,0:4) :: twiddle5
  DATA twiddle5 / zone, zone,  zone,  zone,  zone,  &
                  zone, ei72,  ei144, ei216, ei288, &
                  zone, ei144, ei288, ei72,  ei216, &
                  zone, ei216, ei72,  ei288, ei144, &
                  zone, ei288, ei216, ei144, ei72 /

  COMPLEX(KIND=dz), DIMENSION(0:1,0:1) :: twiddle22
  DATA twiddle22 / zone, zone, &
                   zone, zi    /

  COMPLEX(KIND=dz), DIMENSION(0:1,0:2) :: twiddle23
  DATA twiddle23 / zone, zone, &
                   zone, ei60, &
                   zone, ei120 /

  COMPLEX(KIND=dz), DIMENSION(0:1,0:3) :: twiddle24
  DATA twiddle24 / zone, zone, &
                   zone, ei45, &
                   zone, zi,   &
                   zone, ei135 /
  
  COMPLEX(KIND=dz), DIMENSION(0:1,0:4) :: twiddle25
  DATA twiddle25 / zone, zone, &
                   zone, ei36, &
                   zone, ei72, &
                   zone, ei108, &
                   zone, ei144 /

  COMPLEX(KIND=dz), DIMENSION(0:2,0:2) :: twiddle33
  DATA twiddle33 / zone, zone, zone, &
                   zone, ei40, ei80, &
                   zone, ei80, ei160 /
  
  COMPLEX(KIND=dz), DIMENSION(0:2,0:3) :: twiddle34
  DATA twiddle34 / zone, zone, zone, &
                   zone, ei30, ei60, &
                   zone, ei60, ei120, &
                   zone, zi,  -zone /

  COMPLEX(KIND=dz), DIMENSION(0:2,0:4) :: twiddle35
  DATA twiddle35 / zone, zone, zone, &
                   zone, ei24, ei48, &
                   zone, ei48, ei96, &
                   zone, ei72, ei144, &
                   zone, ei96, ei192 /

  COMPLEX(KIND=dz), DIMENSION(0:3,0:3) :: twiddle44
  DATA twiddle44 / zone, zone, zone,  zone, &
                   zone, ei23, ei45,  ei68, &
                   zone, ei45, zi,    ei135, &
                   zone, ei68, ei135, ei203 /

  COMPLEX(KIND=dz), DIMENSION(0:3,0:4) :: twiddle45
  DATA twiddle45 / zone, zone, zone, zone, &
                   zone, ei18, ei36, ei54, &
                   zone, ei36, ei72, ei108, &
                   zone, ei54, ei108, ei162, &
                   zone, ei72, ei144, ei216 /

  COMPLEX(KIND=dz), DIMENSION(0:4,0:4) :: twiddle55
  DATA twiddle55 / zone, zone, zone,  zone,  zone, &
                   zone, ei14, ei29,  ei43,  ei58, &
                   zone, ei29, ei43,  ei86,  ei115, &
                   zone, ei43, ei86,  ei130, ei173, &
                   zone, ei58, ei115, ei173, ei230 /

CONTAINS

!===============================================================================
! Transposes a double complex matrix
  FUNCTION ZT( Z )
    USE mod_prec, ONLY: dz

    ! Input arguments
    COMPLEX(KIND=dz), DIMENSION(:,:), INTENT(IN) :: Z

    ! Output arguments
    COMPLEX(KIND=dz), DIMENSION(SIZE(Z,2),SIZE(Z,1)) :: ZT

    ! Internal variables
    INTEGER :: i, j, ncols, nrows

    ncols = SIZE(Z,1)
    nrows = SIZE(Z,2)

    DO i = 1, ncols
       DO j = 1, nrows
          ZT(j,i) = Z(i,j)
       END DO ! j
    END DO ! i

    RETURN
  END FUNCTION ZT
    
!===============================================================================
! Generates the twiddle factor e( a \hat{a} / A )
! (Notably, for N = 6 or 7)  
  FUNCTION twiddleN( N )

    USE mod_prec, ONLY: dd, dz
    IMPLICIT NONE

    ! Input argument
    INTEGER, INTENT(IN) :: N

    ! Output argument
    COMPLEX(KIND=dz), DIMENSION(0:(N-1),0:(N-1)) :: twiddleN

    ! Internal arguments
    INTEGER :: i, j
    REAL(KIND=dd) :: x, per, cx, sx

    ! Quick exit
    IF( N == 2 ) THEN
       twiddleN = twiddle2
    ELSE IF( N == 3 ) THEN
       twiddleN = twiddle3
    ELSE IF( N == 4 ) THEN
       twiddleN = twiddle4
    ELSE IF( N == 5 ) THEN
       twiddleN = twiddle5
    ELSE
       
       per = twopi / REAL( N, KIND=dd )
       DO j = 0, N-1
          DO i = 0, N-1
             x = per * REAL( i*j, KIND=dd )
             cx = COS(x)
             sx = SIN(x)
             twiddleN(i,j) = CMPLX( cx, sx, KIND=dz )
          END DO ! i
       END DO ! j

    END IF ! N

    RETURN
  END FUNCTION twiddleN

!===============================================================================
! Generates the twiddle factor e( \hat{a} b / A B )
  FUNCTION twiddleAB( A, B )
    USE mod_prec, ONLY: dd, dz
    IMPLICIT NONE

    ! Input argument
    INTEGER, INTENT(IN) :: A, B

    ! Output argument
    COMPLEX(KIND=dz), DIMENSION(0:(A-1),0:(B-1)) :: twiddleAB

    ! Internal arguments
    INTEGER :: i, j
    REAL(KIND=dd) :: x, per, cx, sx

    ! Quick exit
    IF( ( A == 2 ) .AND. ( B == 2 ) ) THEN
       twiddleAB = twiddle22
    ELSE IF( ( A == 2 ) .AND. ( B == 3 ) ) THEN
       twiddleAB = twiddle23
    ELSE IF( ( A == 2 ) .AND. ( B == 4 ) ) THEN
       twiddleAB = twiddle24
    ELSE IF( ( A == 2 ) .AND. ( B == 5 ) ) THEN
       twiddleAB = twiddle25
    ELSE IF( ( A == 3 ) .AND. ( B == 2 ) ) THEN
       twiddleAB = ZT( twiddle23 )
    ELSE IF( ( A == 3 ) .AND. ( B == 3 ) ) THEN
       twiddleAB = twiddle33
    ELSE IF( ( A == 3 ) .AND. ( B == 4 ) ) THEN
       twiddleAB = twiddle34
    ELSE IF( ( A == 3 ) .AND. ( B == 5 ) ) THEN
       twiddleAB = twiddle35
    ELSE IF( ( A == 4 ) .AND. ( B == 2 ) ) THEN
       twiddleAB = ZT( twiddle24 )
    ELSE IF( ( A == 4 ) .AND. ( B == 3 ) ) THEN
       twiddleAB = ZT( twiddle34 )
    ELSE IF( ( A == 4 ) .AND. ( B == 4 ) ) THEN
       twiddleAB = twiddle44
    ELSE IF( ( A == 4 ) .AND. ( B == 5 ) ) THEN
       twiddleAB = twiddle45
    ELSE IF( ( A == 5 ) .AND. ( B == 2 ) ) THEN
       twiddleAB = ZT( twiddle25 )
    ELSE IF( ( A == 5 ) .AND. ( B == 3 ) ) THEN
       twiddleAB = ZT( twiddle35 )
    ELSE IF( ( A == 5 ) .AND. ( B == 4 ) ) THEN
       twiddleAB = ZT( twiddle45 )
    ELSE IF( ( A == 5 ) .AND. ( B == 5 ) ) THEN
       twiddleAB = twiddle55
    ELSE

       per = twopi / REAL( A*B, KIND=dd )
       DO j = 0, B-1
          DO i = 0, A-1
             x = per * REAL( i*j, KIND=dd )
             cx = COS(x)
             sx = SIN(x)
             twiddleAB(i,j) = CMPLX( cx, sx, KIND=dz )
          END DO ! i
       END DO ! j
       
    END IF ! A, B

    RETURN
  END FUNCTION twiddleAB
       
!===============================================================================
! Factorizes a number into products of 2, 3, 4, 5 and 7

  SUBROUTINE factorize23457( N, nfactors, factors )

    USE mod_prec, ONLY: dd
    IMPLICIT NONE

    ! Input argument
    INTEGER, INTENT(IN) :: N

    ! Output arguments
    INTEGER, INTENT(OUT) :: nfactors
    INTEGER, DIMENSION(nfactors), INTENT(OUT) :: factors

    ! Internal variables
    INTEGER :: i, m, nf
    INTEGER, DIMENSION(N/2) :: f ! Temporary array to hold factors
    INTEGER, DIMENSION(5) :: r
    DATA r /7,5,4,3,2/ ! Note they're in descending order ("greedy" algorithm)

    ! Initialize
    nf = 0
    f(:) = 0

    DO i = 1, SIZE(r,1)
       modloop: DO
          m = MOD( N, r(i) )
          IF( m == 0 ) THEN
             nf = nf + 1
             f(nf) = r(i)
             m = m / r(i)
          ELSE
             EXIT modloop
          END IF ! m
       END DO modloop
    END DO ! i

    ! Assign output
    nfactors = nf
    IF( ALLOCATED(factors) ) DEALLOCATE( factors )
    ALLOCATE( factors(nfactors) )
    factors(1:nfactors) = f(1:nfactors)

    RETURN
  END SUBROUTINE factorize23457
  
!===============================================================================
! Performs the 1-D Fourier transform
! \hat{X} ( \hat(t) ) = \sum_t=0^N-1 X(t) e( t \hat{t} / N )
! Only for N <= 7!
  SUBROUTINE fftXt( X, N, Xhat )

    USE modmain, ONLY: zzero, zone
    USE mod_prec, ONLY: dd, dz
    USE mod_lapack, ONLY: ZGEMV
    IMPLICIT NONE

    ! Input arguments
    INTEGER, INTENT(IN) :: N
    COMPLEX(KIND=dz), DIMENSION(0:N-1), INTENT(IN) :: X

    ! Output argument
    COMPLEX(KIND=dz), DIMENSION(0:N-1), INTENT(OUT) :: Xhat

    ! Dependencies
    INTERFACE
       FUNCTION twiddleN( N )
         USE mod_prec, ONLY: dd, dz
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: N
         COMPLEX(KIND=dz), DIMENSION(0:(N-1),0:(N-1)) :: twiddleN
       END FUNCTION twiddleN
    END INTERFACE

    ! Internal variables
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: mat

    ! Allocate twiddle matrix
    ALLOCATE( mat( 0:(N-1), 0:(N-1) ))
    mat = twiddleN( N )

    ! Calculate Xhat using ZGEMV
    ! \hat{X} = twiddle * X
    CALL ZGEMV( 'N', &
                N, N, &
                zone,  mat,  N, &
                       X,    1, &
                zzero, Xhat, 1 )

    ! Clean up
    DEALLOCATE( mat )
    
    RETURN
  END SUBROUTINE fftXt

!===============================================================================
! Performs the 1-D Fourier transform
! \hat{X} ( \hat(a) + \hat(b) B ) = \sum_b=0^B-1 e( b \hat{b} / B ) Z_hat{a} (b)
! Z_hat{a} (b) = e ( \hat{a} b / A B ) \sum_a=0^A-1 e( a \hat{a} / A ) W_b (a)
  SUBROUTINE fftXab( X, A, B, Xhat )

    USE modmain, ONLY: zzero, zone
    USE mod_prec, ONLY: dd, dz
    USE mod_lapack, ONLY: ZGEMM
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: A, B
    COMPLEX(KIND=dz), DIMENSION(0:(A*B)-1), INTENT(IN) :: X 

    ! Output variable
    COMPLEX(KIND=dz), DIMENSION(0:(A*B)-1), INTENT(IN) :: Xhat

    ! Dependencies
    INTERFACE
       FUNCTION twiddleN( N )
         USE mod_prec, ONLY: dd, dz
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: N
         COMPLEX(KIND=dz), DIMENSION(0:(N-1),0:(N-1)) :: twiddleN
       END FUNCTION twiddleN
       FUNCTION twiddleAB( A, B )
         USE mod_prec, ONLY: dd, dz
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: A, B
         COMPLEX(KIND=dz), DIMENSION(0:(A-1),0:(B-1)) :: twiddleAB
       END FUNCTION twiddleAB
    END INTERFACE
 
    ! Internal variables
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: matA, matB, matAB
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: matWb, matWhat
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: matZa, matZhat

    ! Twiddle matrices
    ALLOCATE( matA(  0:(A-1), 0:(A-1) ))
    ALLOCATE( matB(  0:(B-1), 0:(B-1) ))
    ALLOCATE( matAB( 0:(A-1), 0:(B-1) ))
    matA = twiddleN( A )
    matB = twiddleN( B )
    matAB = twiddleAB( A, B )

    ! Decimate by B
    ALLOCATE( matWb(   0:(A-1), 0:(B-1) ))
    ALLOCATE( matWhat( 0:(A-1), 0:(B-1) ))
    DO j = 0, B-1
       DO i = 0, A-1
          matWb(i,j) = X(j + i*B)
       END DO ! a
    END DO ! b

    ! Twiddle matrix for A
    matA = twiddleN( A )

    ! "Do B different A point Fourier transforms" all at the same time
    ! using ZGEMM:
    ! \hat{W} = twiddleA * Wb(a)
    CALL ZGEMM( 'N', 'N', &
                A, B, A, &
                zone,  matA,    A, &
                       matWb,   A, &
                zzero, matWhat, A )

    ! Multiply by twiddle factor, point by point
    ALLOCATE( matZa( 0:(A-1), 0:(B-1) ))
    DO j = 0, B-1
       DO i = 0, A-1
          matZa(i,j) = matAB(i,j) * matWhat(i,j)
       END DO ! a
    END DO ! b

    ! "Do the A different B point Fourier transforms" all at the same time
    ! using ZGEMM:
    ! \hat{Z} = twiddleB * Z^T (b)
    ALLOCATE( matZhat( 0:(B-1), 0:(A-1) ))
    CALL ZGEMM( 'N', 'T', &
                B, A, B, &
                zone,  matB,    B, &
                       matZa,   A, &
                zzero, matZhat, B )

    ! Finally, transpose and reshape to get the final result
    Xhat = RESHAPE( ZT( matZhat ), (/ A*B /) )

    ! Clean up
    DEALLOCATE( matA )
    DEALLOCATE( matB )
    DEALLOCATE( matAB )
    DEALLOCATE( matWb )
    DEALLOCATE( matWhat )
    DEALLOCATE( matZa )
    DEALLOCATE( matZhat )

    RETURN
  END SUBROUTINE fftXab

!===============================================================================
  SUBROUTINE fft1d_kern_acc( zin, zout, ngrid, sgn )
    IMPLICIT NONE

    ! Arguments
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(IN) :: zin
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(OUT) :: zout
    INTEGER, DIMENSION(1), INTENT(IN) :: ngrid
    INTEGER, INTENT(IN) :: sgn

    ! Internal variables
    INTEGER :: nx, ny, nz, a, b, c, ka, kb, kc
    COMPLEX(KIND=dz) :: z

    IF( ( sgn /= 1 ) .OR. ( sgn /= -1 ) ) THEN
       WRITE(*,*) 'Error(fft_kern_acc): Invalid direction (should be 1 or -1): ', sgn
    END IF

    nx = ngrid(1)
    ny = ngrid(2)
    nz = ngrid(3)

    DO ka = 0, nx
       DO kb = 0, ny
          DO kc = 0, nz

             z = zzero

             IF( sgn == 1 ) THEN

                ! Forward transform
                DO a = 0, nx-1
                   DO b = 0, ny-1
                      DO c = 0, nz-1
                         z = z + tblC( c, kc ) * zin( a + b*nx + c*nx*ny )
                      END DO ! sum_c
                      z = z + z * tblB( b, kb ) * tblABC( a+b*nb, kc )
                   END DO ! sum_b
                   z = z + z * tblA( a, ka ) * tblAB( a, kb )
                END DO ! sum_a

             ELSE

                ! Backward transform
                DO a = 0, nx-1
                   DO b = 0, ny-1
                      DO c = 0, nz-1
                         z = z + CONJG( tblC( c, kc ) ) * zin( a + b*nx + c*nx*ny )
                      END DO ! sum_c
                      z = z + z * CONJG( tblB( b, kb ) * tblABC( a+b*nb, kc ) )
                   END DO ! sum_b
                   z = z + z * CONJG( tblA( a, ka ) * tblAB( a, kb ) )
                END DO ! sum_a

             END IF ! sgn

             zout( kc + kb*nz + ka*nx*ny ) = z/nx/ny/nz

          END DO ! kc
       END DO ! kb
    END DO ! ka

    RETURN
  END SUBROUTINE fft_kern_acc

!===============================================================================

END MODULE mod_fft_acc
