MODULE mod_sparse
  IMPLICIT NONE

CONTAINS

!===============================================================================
! Finds the nonzero rows and columns of a double complex matrix
! mat(1:nrows,1:ncols), and returns them as logical arrays
! lrownz(1:nrownz) and lcolnz(1:lcolnz), respectively.  
  SUBROUTINE zge2sp_findnnz( nrows, ncols, mat, nrownz, lrownz, ncolnz, lcolnz)

    USE mod_prec, ONLY: dd, dz
    USE mod_mpi_grid, ONLY: iproc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: nrows, ncols
    COMPLEX(KIND=dz), DIMENSION(nrows,ncols), INTENT(IN) :: mat
    INTEGER, INTENT(OUT) :: nrownz, ncolnz
    LOGICAL, DIMENSION(nrownz), INTENT(OUT) :: lrownz
    LOGICAL, DIMENSION(ncolnz), INTENT(OUT) :: lcolnz

    ! Internal variables
    LOGICAL, DIMENSION(nrows,ncols) :: tblnz
    INTEGER, DIMENSION(nrows) :: col
    INTEGER, DIMENSION(ncols) :: row
    INTEGER :: i, j

    ! Initialize vars
    nrownz = 0
    ncolnz = 0
    lrownz(:) = .FALSE.
    lcolnz(:) = .FALSE.
    tblnz(:,:) = .FALSE.
    col(:) = 0
    row(:) = 0

    ! Fill in tblnz
    ! TODO: Parallelize
    DO j = 1, ncols
       DO i = 1, nrows
          IF( ABS(mat(i,j)) /= 0._dd ) tblnz(i,j) = .TRUE.
       END DO
       col(j) = COUNT( tblnz(:,j) )
    END DO
    DO i = 1, nrows
       row(i) = COUNT( tblnz(i,:) )
    END DO

    ! Count nrownz and ncolnz
    nrownz = COUNT( row(:) )
    ncolnz = COUNT( col(:) )

    ! Fill in lrownz
    ! TODO: Parallelize
    DO i = 1, nrownz
       IF( row(i) /= 0 ) lrownz(i) = .TRUE.
    END DO

    ! Fill in lcolnz
    ! TODO: Parallelize
    DO j = 1, ncolnz
       IF( col(j) /= 0 ) lcolnz(j) = .TRUE.
    END DO

    RETURN
  END SUBROUTINE zge2sp_findnnz

!===============================================================================
! Packs a sparse matrix mat(1:nrows,1:ncols) into matnz(1:nrownz,1:ncolnz)
! Call zge2sp_findnnz() first before calling this subroutine!
  SUBROUTINE zge2sp_pack( nrows, ncols, mat,
                          nrownz, lrownz, ncolnz, lcolnz, matnz )

    USE mod_prec, ONLY: dd, dz
    USE mod_mpi_grid, ONLY: iproc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: nrows, ncols, nrownz, ncolnz
    COMPLEX(KIND=dz), DIMENSION(nrows,ncols), INTENT(IN) :: mat
    LOGICAL, DIMENSION(nrownz), INTENT(IN) :: lrownz
    LOGICAL, DIMENSION(ncolnz), INTENT(IN) :: lcolnz
    COMPLEX(KIND=dz), DIMENSION(nrownz,ncolnz), INTENT(OUT) :: matnz

    ! Internal variables
    INTEGER :: i, j, irow, icol

    ! Permute the non-zero data into matnz
    ! TODO: Parallelize
    DO icol = 1, ncolnz
       j = icolnz(icol)
       DO irow = 1, nrownz
          i = icolnz(irow)
          matnz(irow,icol) = mat(i,j)
       END DO ! irow
    END DO ! icol

    RETURN
  END SUBROUTINE zge2sp_pack

!===============================================================================

END MODULE mod_sparse
