MODULE mod_fft_acc

  USE mod_prec, ONLY: dd, dz
  USE modmain, ONLY: zzero, zone, zi, pi, twopi, sqtwo
  IMPLICIT NONE

  ! FFT methods

  INTEGER, PARAMETER :: fft_1d_lib           = 10
  INTEGER, PARAMETER :: fft_1d_zgemm_direct  = 11
  INTEGER, PARAMETER :: fft_1d_zgemm_reshape = 12
  INTEGER, PARAMETER :: fft_1d_zgemm_factor  = 12
  INTEGER, PARAMETER :: fft_1d_kron_direct   = 15

  INTEGER, PARAMETER :: fft_2d_lib           = 20
  INTEGER, PARAMETER :: fft_2d_zgemm_direct  = 21
  INTEGER, PARAMETER :: fft_2d_kron_direct   = 25

  INTEGER, PARAMETER :: fft_3d_lib           = 30
!  INTEGER, PARAMETER :: fft_3d_zgemm_direct  = 31
  INTEGER, PARAMETER :: fft_3d_kron_direct   = 35

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
  COMPLEX(KIND=dz), PARAMETER :: ei120 = CMPLX( -sin30,  cos30, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei160 = CMPLX( -cos20,  sin20, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: zmone = -zone
  COMPLEX(KIND=dz), PARAMETER :: ei200 = CMPLX( -cos20, -sin20, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei240 = CMPLX( -sin30, -cos30, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: zmi   = -zi
  COMPLEX(KIND=dz), PARAMETER :: ei280 = CMPLX( -sin10, -cos10, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei300 = CMPLX(  sin30, -cos30, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei320 = CMPLX(  cos40, -sin40, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei330 = CMPLX(  cos30, -sin30, KIND=dz )

  ! For size 4, 8, 16
  REAL(KIND=dd), PARAMETER :: cos23 = 0.923879532511287_dd ! cos(pi/8)
  REAL(KIND=dd), PARAMETER :: cos45 = sqtwo/2.0_dd
  REAL(KIND=dd), PARAMETER :: sin23 = 0.382683432365090_dd ! sin(pi/8)
  COMPLEX(KIND=dz), PARAMETER :: ei23  = CMPLX(  cos23,  sin23, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei45  = CMPLX(  cos45,  cos45, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei68  = CMPLX(  sin23,  cos23, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei135 = CMPLX( -cos45,  cos45, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei203 = CMPLX( -cos23, -sin23, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei215 = CMPLX(  cos45, -cos45, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei315 = CMPLX(  cos45, -cos45, KIND=dz )
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
  REAL(KIND=dd), PARAMETER :: cos39 = 0.770513242775789_dd ! cos(11pi/50)
  REAL(KIND=dd), PARAMETER :: cos42 = 0.743144825477394_dd ! cos(7pi/30)
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
  REAL(KIND=dd), PARAMETER :: sin39 = 0.637423989748690_dd ! sin(11pi/50)
  REAL(KIND=dd), PARAMETER :: sin42 = 0.669130606358858_dd ! sin(7pi/30)
  REAL(KIND=dd), PARAMETER :: sin43 = 0.684547105928689_dd ! sin(6pi/25)
  COMPLEX(KIND=dz), PARAMETER :: ei14  = CMPLX(  cos14,  sin14, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei18  = CMPLX(  cos18,  sin18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei24  = CMPLX(  cos24,  sin24, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei29  = CMPLX(  cos29,  sin29, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei36  = CMPLX(  cos36,  sin36, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei43  = CMPLX(  cos43,  sin43, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei48  = CMPLX(  sin42,  cos42, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei54  = CMPLX(  sin36,  cos36, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei58  = CMPLX(  sin32,  cos32, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei72  = CMPLX(  sin18,  cos18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei86  = CMPLX(  sin4,   cos4,  KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei96  = CMPLX( -sin6,   cos6,  KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei108 = CMPLX( -sin18,  cos18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei115 = CMPLX( -sin25,  cos25, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei130 = CMPLX( -sin39,  cos39, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei144 = CMPLX( -cos36,  sin36, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei162 = CMPLX( -cos18,  sin18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei173 = CMPLX( -cos7,   sin7,  KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei192 = CMPLX( -cos12, -sin12, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei216 = CMPLX( -cos36, -sin36, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei230 = CMPLX( -cos40, -sin40, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei252 = CMPLX( -sin18, -cos18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei288 = CMPLX(  sin18, -cos18, KIND=dz )
  COMPLEX(KIND=dz), PARAMETER :: ei324 = CMPLX(  cos36, -sin36, KIND=dz )

  COMPLEX(KIND=dz), DIMENSION(0:1,0:1) :: twiddle2
  DATA twiddle2 / zone, zone, &
                  zone, zmone /
  
  COMPLEX(KIND=dz), DIMENSION(0:2,0:2) :: twiddle3
  DATA twiddle3 / zone, zone,  zone,  &
                  zone, ei120, ei240, &
                  zone, ei240, ei120 /

  COMPLEX(KIND=dz), DIMENSION(0:3,0:3) :: twiddle4
  DATA twiddle4 / zone, zone, zone, zone, &
                  zone, zi,   zmone,zmi,   &
                  zone, zmone,zone, zmone, &
                  zone, zmi,  zmone,zi /

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
                   zone, zi,   zmone /

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

#ifdef _CUFFT_
  LOGICAL :: lplanned = .FALSE.
  INTEGER :: plan
#endif /* _CUFFT_ */

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
! Transposes and conjugates a double complex matrix
  FUNCTION ZH( Z )
    USE mod_prec, ONLY: dz

    ! Input arguments
    COMPLEX(KIND=dz), DIMENSION(:,:), INTENT(IN) :: Z

    ! Output arguments
    COMPLEX(KIND=dz), DIMENSION(SIZE(Z,2),SIZE(Z,1)) :: ZH

    ! Dependencies
    INTRINSIC :: CONJG

    ! Internal variables
    INTEGER :: i, j, ncols, nrows

    ncols = SIZE(Z,1)
    nrows = SIZE(Z,2)

    DO i = 1, ncols
       DO j = 1, nrows
          ZH(j,i) = CONJG(Z(i,j))
       END DO ! j
    END DO ! i

    RETURN
  END FUNCTION ZH

!===============================================================================
! Generates the twiddle factor e( a \hat{a} / A )
! (Notably, for N = 6 or 8)
! dir is direction of FFT (+1 = forward, -1 = backward)
  FUNCTION twiddleN( N, dir )

    USE mod_prec, ONLY: dd, dz
    IMPLICIT NONE

    ! Input argument
    INTEGER, INTENT(IN) :: N, dir
    
    ! Output argument
    COMPLEX(KIND=dz), DIMENSION(0:(N-1),0:(N-1)) :: twiddleN

    ! Internal arguments
    INTEGER :: i, j
    REAL(KIND=dd) :: x, per, cx, sx

    ! Quick exit
    IF( (dir /= 1) .OR. (dir /= +1) ) THEN
       WRITE(*,*) 'Error[twiddleN]: direction should be +1 (forward) &
                  &or -1 (backward)'
       STOP
    END IF
    IF( N == 2 ) THEN
       twiddleN = twiddle2
    ELSE IF( N == 3 ) THEN
       IF( dir == +1 ) THEN
          twiddleN = twiddle3
       ELSE ! dir == -1
          twiddleN = CONJG( twiddle3 )
       END IF ! dir
    ELSE IF( N == 4 ) THEN
       IF( dir == +1 ) THEN
          twiddleN = twiddle4
       ELSE ! dir == -1
          twiddleN = CONJG( twiddle4 )
       END IF ! dir
    ELSE IF( N == 5 ) THEN
       IF( dir == +1 ) THEN
          twiddleN = twiddle5
       ELSE ! dir == -1
          twiddleN = CONJG( twiddle5 )
       END IF ! dir
    ELSE
       
       per = dir * twopi / REAL( N, KIND=dd )
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
! dir is direction of FFT (+1 = forward, -1 = backward)
  FUNCTION twiddleAB( A, B, dir )
    USE mod_prec, ONLY: dd, dz
    IMPLICIT NONE

    ! Input argument
    INTEGER, INTENT(IN) :: A, B, dir

    ! Output argument
    COMPLEX(KIND=dz), DIMENSION(0:(A-1),0:(B-1)) :: twiddleAB

    ! Internal arguments
    INTEGER :: i, j
    REAL(KIND=dd) :: x, per, cx, sx

    ! Quick exit
    IF( (dir /= 1) .OR. (dir /= +1) ) THEN
       WRITE(*,*) 'Error[twiddleAB]: direction should be +1 (forward) &
                  &or -1 (backward)'
       STOP
    END IF
    IF( ( A == 2 ) .AND. ( B == 2 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle22
       ELSE
          twiddleAB = CONJG( twiddle22 )
       END IF ! dir
    ELSE IF( ( A == 2 ) .AND. ( B == 3 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle23
       ELSE
          twiddleAB = CONJG( twiddle23 )
       END IF ! dir
    ELSE IF( ( A == 2 ) .AND. ( B == 4 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle24
       ELSE
          twiddleAB = CONJG( twiddle24 )
       END IF ! dir
    ELSE IF( ( A == 2 ) .AND. ( B == 5 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle25
       ELSE
          twiddleAB = CONJG( twiddle25 )
       END IF ! dir
    ELSE IF( ( A == 3 ) .AND. ( B == 2 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = ZT( twiddle23 )
       ELSE
          twiddleAB = ZH( twiddle23 )
       END IF ! dir
    ELSE IF( ( A == 3 ) .AND. ( B == 3 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle33
       ELSE
          twiddleAB = CONJG( twiddle33 )
       END IF ! dir       
    ELSE IF( ( A == 3 ) .AND. ( B == 4 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle34
       ELSE
          twiddleAB = CONJG( twiddle34 )
       END IF ! dir
    ELSE IF( ( A == 3 ) .AND. ( B == 5 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle35
       ELSE
          twiddleAB = CONJG( twiddle35 )
       END IF ! dir
    ELSE IF( ( A == 4 ) .AND. ( B == 2 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = ZT( twiddle24 )
       ELSE
          twiddleAB = ZH( twiddle24 )
       END IF ! dir
    ELSE IF( ( A == 4 ) .AND. ( B == 3 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = ZT( twiddle34 )
       ELSE
          twiddleAB = ZH( twiddle34 )
       END IF ! dir
    ELSE IF( ( A == 4 ) .AND. ( B == 4 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle44
       ELSE
          twiddleAB = CONJG( twiddle44 )
       END IF ! dir
    ELSE IF( ( A == 4 ) .AND. ( B == 5 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle45
       ELSE
          twiddleAB = CONJG( twiddle45 )
       END IF ! dir
    ELSE IF( ( A == 5 ) .AND. ( B == 2 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = ZT( twiddle25 )
       ELSE
          twiddleAB = ZH( twiddle25 )
       END IF ! dir
    ELSE IF( ( A == 5 ) .AND. ( B == 3 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = ZT( twiddle35 )
       ELSE
          twiddleAB = ZH( twiddle35 )
       END IF ! dir
    ELSE IF( ( A == 5 ) .AND. ( B == 4 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = ZT( twiddle45 )
       ELSE
          twiddleAB = ZH( twiddle45 )
       END IF ! dir
    ELSE IF( ( A == 5 ) .AND. ( B == 5 ) ) THEN
       IF( dir == +1 ) THEN
          twiddleAB = twiddle55
       ELSE
          twiddleAB = CONJG( twiddle55 )
       END IF ! dir
    ELSE

       per = dir * twopi / REAL( A*B, KIND=dd )
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

  SUBROUTINE factorize( N, nfactors, factors )

    USE mod_prec, ONLY: dd
    IMPLICIT NONE

    ! Input argument
    INTEGER, INTENT(IN) :: N

    ! Output arguments
    INTEGER, INTENT(OUT) :: nfactors
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: factors

    ! Internal variables
    INTEGER :: i, m, nf
    INTEGER, DIMENSION(N/2) :: f ! Temporary array to hold factors
    INTEGER, DIMENSION(6) :: r
    DATA r /8,6,5,4,3,2/ ! Note they're in descending order ("greedy" algorithm)

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
  END SUBROUTINE factorize
  
!===============================================================================
! Performs the 1-D fast Fourier transform from
! W.M. Gentleman and G. Sande,
! "Fast Fourier Transforms -- for Fun and Profit"
! Proceedings of the 1966 AFIPS Fall Joint Computer Conference
! Pages 563-578
! doi:10.1145/1464291.1464352
! \hat{X} ( \hat(t) ) = \sum_t=0^N-1 X(t) e( t \hat{t} / N )
! Only for N <= 8!
  SUBROUTINE fftXt( X, N, Nvec, dir, Xhat )

    USE modmain, ONLY: zzero, zone
    USE mod_prec, ONLY: dd, dz
    USE mod_gpu, ONLY: ZGEMM_acc
    IMPLICIT NONE

    ! Input arguments
    INTEGER, INTENT(IN) :: N, Nvec, dir
    COMPLEX(KIND=dz), DIMENSION(0:N-1,Nvec), INTENT(IN) :: X

    ! Output argument
    COMPLEX(KIND=dz), DIMENSION(0:N-1,Nvec), INTENT(OUT) :: Xhat

    ! Internal variables
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: mat
    COMPLEX(KIND=dz) :: perN

    ! Quick exit
    IF( (dir /= 1) .OR. (dir /= +1) ) THEN
       WRITE(*,*) 'Error[fftXt]: direction should be +1 (forward) &
                  &or -1 (backward)'
       STOP
    END IF

    ! Allocate twiddle matrix
    ALLOCATE( mat( 0:(N-1), 0:(N-1) ))
    mat = twiddleN( N, dir )

    ! Calculate multiple Xhat vectors using ZGEMM
    ! \hat{X} = twiddle * X
    IF( dir == +1 ) THEN
       CALL ZGEMM_acc( 'N', 'N', &
                        N, Nvec, N, &
                        zone,  mat,  N, &
                               X,    N, &
                        zzero, Xhat, N )
    ELSE
       perN = CMPLX( 1._dd/REAL( N, KIND=dd ), 0._dd, KIND=dz )
       CALL ZGEMM_acc( 'N', 'N', &
                       N, Nvec, N, &
                       perN,  mat,  N, &
                              X,    N, &
                       zzero, Xhat, N )
    END IF ! dir

    ! Clean up
    DEALLOCATE( mat )
    
    RETURN
  END SUBROUTINE fftXt

!===============================================================================
! Performs the 2-D fast Fourier transform from Gentleman & Sande
! \hat{X} ( \hat(a) + \hat(b) B ) = \sum_b=0^B-1 e( b \hat{b} / B ) Z_hat{a} (b)
! Z_hat{a} (b) = e ( \hat{a} b / A B ) \sum_a=0^A-1 e( a \hat{a} / A ) W_b (a)
  SUBROUTINE fftXab( X, A, B, dir, Xhat )

    USE modmain, ONLY: zzero, zone
    USE mod_prec, ONLY: dd, dz
    USE mod_gpu, ONLY: ZGEMM_acc
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: A, B, dir
    COMPLEX(KIND=dz), DIMENSION(0:(A*B)-1), INTENT(IN) :: X 

    ! Output variable
    COMPLEX(KIND=dz), DIMENSION(0:(A*B)-1), INTENT(OUT) :: Xhat

    ! Dependencies
    INTRINSIC :: RESHAPE

    ! Internal variables
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: matA, matB, matAB
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: matWb, matWhat
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: matZa, matZhat
    COMPLEX(KIND=dz) :: perN
    INTEGER :: i, j

    ! Twiddle matrices
    ALLOCATE( matA(  0:(A-1), 0:(A-1) ))
    ALLOCATE( matB(  0:(B-1), 0:(B-1) ))
    ALLOCATE( matAB( 0:(A-1), 0:(B-1) ))
    matA = twiddleN( A, dir )
    matB = twiddleN( B, dir )
    matAB = twiddleAB( A, B, dir )

    ! Decimate by B
    ALLOCATE( matWb(   0:(A-1), 0:(B-1) ))
    ALLOCATE( matWhat( 0:(A-1), 0:(B-1) ))
    DO j = 0, B-1
       DO i = 0, A-1
          matWb(i,j) = X(j + i*B)
       END DO ! a
    END DO ! b

    ! "Do B different A point Fourier transforms" all at the same time
    ! using ZGEMM:
    ! \hat{W} = twiddleA * Wb(a)
    CALL ZGEMM_acc( 'N', 'N', &
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
    IF( dir == +1 ) THEN
       CALL ZGEMM_acc( 'N', 'T', &
                       B, A, B, &
                       zone,  matB,    B, &
                              matZa,   A, &
                       zzero, matZhat, B )
    ELSE
       perN = CMPLX( 1._dd/REAL( A*B, KIND=dd ), 0._dd, KIND=dz )
       CALL ZGEMM_acc( 'N', 'T', &
                       B, A, B, &
                       perN,  matB,    B, &
                              matZa,   A, &
                       zzero, matZhat, B )
    END IF ! dir

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
! Perform 1-D FFT using the selected method
! If no method is given, defaults to fft_1d_lib (call library)
! For fft_1d_zgemm_reshape, param = radix
! For fft_1d_kron_direct,   param = Nvec 

  SUBROUTINE fft1d_kern_acc( zin, zout, N, dir, method, param )
    USE mod_kron, ONLY: zkronmult1
    IMPLICIT NONE

    ! Arguments
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(IN), TARGET :: zin
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(OUT) :: zout
    INTEGER, INTENT(IN) :: N, dir
    INTEGER, INTENT(IN), OPTIONAL :: method, param

    ! Internal variables
    COMPLEX(KIND=dz), DIMENSION(:), ALLOCATABLE :: tmpvec
    COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: tmpmat, twiddle
    INTEGER :: i

    IF( (.NOT. PRESENT(method)) .OR. (method == fft_1d_lib) ) THEN

       ! Call whatever FFT library interface is enabled
       CALL zfftifc_gpu_exec( 1, (/ N /), dir, tmpvec )

    ELSE IF( method == fft_1d_zgemm_direct ) THEN

       ! Call Gentleman-Sande 1-D FFT routine
       CALL fftXt( RESHAPE( zin(1:N), (/ N, 1 /)), N, 1, dir, zout )

    ELSE IF( method == fft_1d_zgemm_reshape ) THEN
       ! For this method, param = FFT radix

       ! Allocate temporary matrix
       ALLOCATE( tmpmat(param,N/param) )

       ! Call Gentleman-Sande 1-D FFT routine
       CALL fftXt( RESHAPE( zin(1:N), (/ param, N/param /)), &
                   param, N/param, dir, tmpmat )

       ! Copy temporary matrix to output
       zout(1:N) = RESHAPE( tmpmat(1:param,1:N/param), (/ N /) )

       ! Clean up
       DEALLOCATE( tmpmat )

    ELSE IF( method == fft_1d_zgemm_factor ) THEN

       WRITE(*,*) 'fft_1d_zgemm_factor method not yet implemented'
       STOP

    ELSE IF( method == fft_1d_kron_direct ) THEN
       ! For this method, param = Nvec
       
       ! Allocate temporary matrix
       ALLOCATE( tmpmat(N,param) )

       ! Form twiddle matrix
       ALLOCATE( twiddle(N,N) )
       twiddle = twiddleN( N, dir )

       ! Perform Kronecker product
       CALL zkronmult1( N, N, twiddle, N, param, &
                        RESHAPE( zin(1:N), (/ N, param /)), &
                        tmpmat )

       ! Copy temporary matrix to output
       zout(1:N*param) = RESHAPE( tmpmat(1:N,1:param), (/ N*param /) )

       ! Clean up
       DEALLOCATE( tmpmat )
       DEALLOCATE( twiddle )

    END IF ! method

    RETURN
  END SUBROUTINE fft1d_kern_acc

!===============================================================================

  SUBROUTINE zfftifc_gpu_init( nd, ngrid )

#if defined(_CUFFT_)
    USE cufft
#elif defined(_ROCFFT_)
    USE rocfft
#endif /* _CUFFT_ || _ROCFFT_ */

    USE mod_prec, ONLY: dz
    USE ISO_FORTRAN_ENV, ONLY: u => error_unit
    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: nd
    INTEGER, DIMENSION(nd), INTENT(IN) :: ngrid

    ! Internal variables
    INTEGER :: ierr

#ifdef _CUFFT_

    SELECT CASE(nd)
    CASE(1)
       ierr = cufftPlan1D( plan, ngrid(1), CUFFT_Z2Z, 1 )
       IF( ierr /= 0 ) THEN
          WRITE(u,*) 'Error[zfftifc_gpu_init]: cufftPlan1D returned ', ierr
          STOP
       END IF
    CASE(2)
       ierr = cufftPlan2D( plan, ngrid(1), ngrid(2), CUFFT_Z2Z )
       IF( ierr /= 0 ) THEN
          WRITE(u,*) 'Error[zfftifc_gpu_init]: cufftPlan2D returned ', ierr
          STOP
       END IF
    CASE(3)
       ierr = cufftPlan3D( plan, ngrid(1), ngrid(2), ngrid(3), CUFFT_Z2Z )
       IF( ierr /= 0 ) THEN
          WRITE(u,*) 'Error[zfftifc_gpu_init]: cufftPlan3D returned ', ierr
          STOP
       END IF
    CASE DEFAULT
       WRITE(u,*) 'Error[zfftifc_gpu_init]: unsupported number of dimensions nd = ', nd
       STOP
    END SELECT
    lplanned = .TRUE.

#elif defined(_ROCFFT_)

#endif /* _CUFFT_ || _ROCFFT_ */

    RETURN
  END SUBROUTINE zfftifc_gpu_init

!===============================================================================

  SUBROUTINE zfftifc_gpu_exec( nd, ngrid, dir, z )

#if defined(_CUFFT_)
    USE cufft
#elif defined(_ROCFFT_)
    USE rocfft
#endif /* _CUFFT_ || _ROCFFT_ */

    USE mod_prec, ONLY: dz
    USE ISO_FORTRAN_ENV, ONLY: u => error_unit
    IMPLICIT NONE
    
    ! Arguments
    INTEGER, INTENT(IN) :: nd, dir
    INTEGER, DIMENSION(nd), INTENT(IN) :: ngrid
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(INOUT) :: z

    ! Internal variables
    INTEGER :: ierr
    LOGICAL :: toggle = .FALSE.

#ifdef _CUFFT_

    ! cuFFT interface
    IF( .NOT. lplanned ) THEN
       WRITE(u,*), 'Warning[zfftifc_gpu]: please create a plan first using zfftifc_gpu_init()'
       CALL zfftifc_gpu_init( nd, ngrid )
       toggle = .TRUE.
    END IF
    IF( dir == 1 ) THEN
       ierr = cufftExecZ2Z( plan, z, z, CUFFT_FORWARD )
    ELSE IF( dir == -1 ) THEN
       ierr = cufftExecZ2Z( plan, z, z, CUFFT_INVERSE )
       !$ACC KERNELS
       z = z / PRODUCT(ngrid)
       !$ACC END KERNELS
    ELSE
       WRITE(u,*), 'Error[zfftifc_gpu]: unknown direction dir=', dir
       STOP
    END IF
    IF( ierr /= 0 ) THEN
       WRITE(u,*) 'Error[zfftifc_gpu]: cufftExecZ2Z returned ', ierr
       STOP
    END IF
    IF( toggle ) THEN
       CALL zfftifc_gpu_fin()
       toggle = .FALSE.
    END IF

#elif defined(_ROCFFT_)

    ! rocFFT interface
    WRITE(*,*) 'zfftifc_gpu: rocFFT interface not implemented yet'

#endif /* _CUFFT_ || _ROCFFT_ */

    RETURN
  END SUBROUTINE zfftifc_gpu_exec

!===============================================================================

  SUBROUTINE zfftifc_gpu_fin()

#if defined(_CUFFT_)
    USE cufft
#elif defined(_ROCFFT_)
    USE rocfft
#endif /* _CUFFT_ || _ROCFFT_ */

    USE mod_prec, ONLY: dz
    USE ISO_FORTRAN_ENV, ONLY: u => error_unit
    IMPLICIT NONE    

    ! Internal variables
    INTEGER :: ierr

#if defined(_CUFFT_)

    ierr = cufftDestroy( plan )
    IF( ierr /= 0 ) THEN
       WRITE(u,*) 'Error[zfftifc_gpu_fin]: cufftDestroy returned ', ierr
       STOP
    END IF
    lplanned = .FALSE.

#elif defined(_ROCFFT_)

#endif /* _CUFFT_ || _ROCFFT_ */

    RETURN
  END SUBROUTINE zfftifc_gpu_fin

!===============================================================================

END MODULE mod_fft_acc
