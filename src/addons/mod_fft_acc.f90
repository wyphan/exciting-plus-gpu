MODULE mod_fft_acc
! Implements the Fast Fourier Transform algorithm from
! W.M. Gentleman and G. Sande,
! "Fast Fourier Transforms -- for Fun and Profit"
! Proceedings of the 1966 AFIPS Fall Joint Computer Conference
! Pages 563-578
! doi:10.1145/1464291.1464352
!===============================================================================

  USE mod_prec, ONLY: dd, dz
  USE modmain, ONLY: zzero
  IMPLICIT NONE

  ! Mathematical constants
  REAL(KIND=dd), PARAMETER    :: pi = 2._dd * ASIN(1._dd)
  COMPLEX(kind=dz), PARAMETER :: ii = (0._dd, 1._dd)

  ! Block size for expi()
  INTEGER, PARAMETER :: nb = 2

  ! Twiddle factor tables
  COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: tblA, tblB, tblC
  COMPLEX(KIND=dz), DIMENSION(:,:), ALLOCATABLE :: tblAB, tblABC

CONTAINS

!===============================================================================
! Allocates and fills in the arrays for the twiddle factor tables
  SUBROUTINE inittable(ngrid)
    IMPLICIT NONE

    ! Input argument
    INTEGER, DIMENSION(3), INTENT(IN) :: ngrid

    ! Internal variables
    INTEGER :: nx, ny, nz, nAB
    INTEGER :: a, b, c, ka, kb, kc, ib
    REAL(KIND=dd) :: dx, dy, dz, dkx, dky, dkz
    REAL(KIND=dd), DIMENSION(0:nb-1) :: da, db, dc, kx, ky, kz, kdotx

    nx = ngrid(1)
    ny = ngrid(2)
    nz = ngrid(3)

    nAB = ngrid(1)*ngrid(2)

    ALLOCATE( tblA(   0:(nx-1),  0:(nx-1) ))
    ALLOCATE( tblB(   0:(ny-1),  0:(ny-1) ))
    ALLOCATE( tblC(   0:(nz-1),  0:(nz-1) ))
    ALLOCATE( tblAB(  0:(nx-1),  0:(ny-1) ))
    ALLOCATE( tblABC( 0:(nAB-1), 0:(nz-1) ))

    ! Prepare increments
    dx = 1._dd / REAL(nx, KIND=dd)
    dy = 1._dd / REAL(ny, KIND=dd)
    dz = 1._dd / REAL(nz, KIND=dd)
    dkx = 1._dd / REAL(nx, KIND=dd)
    dky = 1._dd / REAL(ny, KIND=dd)
    dkz = 1._dd / REAL(nz, KIND=dd)

    ! Fill in twiddle factor table e( a \hat{a} / A )
    DO ib = 0, nb
       da(ib) = ib * dx ! a/A
    END DO
    DO ka = 0, nx
       kx(:) = REAL(ka, KIND=dd) ! \hat{a}
       DO a = 0, nx, nb
          kdotx(:) = kx(:) * ( a*dx + da(:) )
          tblA( a:(a+nb-1), ka ) = expi( kdotx(:) )
       END DO
    END DO

    ! Fill in twiddle factor table e( a \hat{b} /( A B ) )
    DO kb = 0, ny
       ky(:) = kb * dky ! \hat{b}/B
       DO a = 0, nx, nb
          kdotx(:) = ky(:) * ( a*dx + da(:) )
          tblAB( a:(a+nb-1), kb ) = expi( kdotx(:) )
       END DO
    END DO

    ! Fill in twiddle factor table e( \hat{c} (a + bA) / ( A B C ) )
    DO kc = 0, nz
       kz(:) = kc * dkz ! \hat{c}/C
       DO b = 0, ny
          DO a = 0, nx, nb
             kdotx(:) = kz(:) * dxb*( a*dx + da(:) + REAL(b,kind=dd) )
             tblABC( (b*nx+a):(b*nx+a+nb-1), kc ) = expi( kdotx(:) )
          END DO
       END DO
    END DO

    ! Fill in twiddle factor table e( b \hat{b} / B )
    DO ib = 0, nb
       db(ib) = ib * dy ! b/B
    END DO
    DO kb = 0, ny
       ky(:) = REAL(kb, KIND=dd) ! \hat{b}
       DO b = 0, ny, nb
          kdotx(:) = ky(:) * ( b*dy + db(:) )
          tblB( b:(b+nb-1), kb ) = expi( kdotx(:) )
       END DO
    END DO

    ! Fill in twiddle factor table e( c \hat{c} / C )
    DO ib = 0, nb
       dc(ib) = ib * dz ! c/C
    END DO
    DO kc = 0, nz
       kz(:) = REAL(kc, KIND=dd) ! \hat{c}
       DO c = 0, nz, nb
          kdotx(:) = kz(:) * ( c*dz + dc(:) )
          tblC( c:(c+nb-1), kc ) = expi( kdotx(:) )
       END DO
    END DO

    RETURN
  END SUBROUTINE inittable

!===============================================================================
! Deallocates the arrays for the twiddle factor tables
  SUBROUTINE freetable
    IMPLICIT NONE

    DEALLOCATE( tblA )
    DEALLOCATE( tblB )
    DEALLOCATE( tblC )
    DEALLOCATE( tblAB )
    DEALLOCATE( tblABC )

    RETURN
  END SUBROUTINE freetable

!===============================================================================
! Euler's formula
! exp(ix) = cos(x) + i sin(x)
  FUNCTION expi(x)
    IMPLICIT NONE

    ! Input argument
    REAL(KIND=dd), INTENT(IN), DIMENSION(0:(nb-1)) :: x

    ! Output argument
    COMPLEX(KIND=dz) :: expi(0:(nb-1))

    ! Internal variables
    REAL(KIND=dd), DIMENSION(nb) :: c, s

    c(:) = COS( x(:) )
    s(:) = SIN( x(:) )
    expi(:) = CMPLX( c(:), s(:), KIND=dz)

    END DO

    RETURN
  END FUNCTION expi

!===============================================================================
! Performs the 3-D Fourier transform
  SUBROUTINE fft_kern( zin, zout, ngrid, sgn )
    IMPLICIT NONE

    ! Arguments
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(IN) :: zin
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(OUT) :: zout
    INTEGER, DIMENSION(3), INTENT(IN) :: ngrid
    INTEGER :: INTENT(IN) :: sgn

    ! Internal variables
    INTEGER :: nx, ny, nz, a, b, c, ka, kb, kc
    COMPLEX(KIND=dz) :: z

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
  END SUBROUTINE fft_kern

!===============================================================================

END MODULE mod_fft_acc
