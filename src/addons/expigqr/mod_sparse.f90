MODULE mod_sparse

  USE mod_prec, ONLY: dz
  IMPLICIT NONE

  COMPLEX(KIND=dz), ALLOCATABLE :: gntuju_packed(:,:,:,:)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nmt, nareanz
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: igntujunz, tblgntujunz

CONTAINS

!-------------------------------------------------------------------------------
! Finds the non-zero rows of a symmetric double complex matrix mat(nrows,nrows)
! as nrownz, and fills in the translation table icolnz(nrownz)
! TODO: OpenACC port
SUBROUTINE zsy2sp_findnnz( nrows, mat, nrownz, icolnz )

  USE mod_prec, ONLY: dd, dz
  USE mod_mpi_grid, ONLY: iproc

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: nrows
  COMPLEX(KIND=dz), DIMENSION(nrows,nrows), INTENT(IN) :: mat
  INTEGER, INTENT(OUT) :: nrownz
  INTEGER, DIMENSION(nrownz), INTENT(OUT) :: icolnz

  ! Internal variables
  INTEGER, DIMENSION(nrows,nrows) :: tblnz
  INTEGER, DIMENSION(nrows) :: col
  INTEGER :: i, j, mask

  ! Find nonzero rows on the lower triangular matrix
  ! TODO: Parallelize
  tblnz(:,:) = 0
  col(:) = 0
  DO j = 1, nrows

     icolnz(j) = 0

     DO i = 1, j

        tblnz(i,j) = 0
        IF( mat(i,j) /= 0._dd ) THEN
           tblnz(i,j) = 1
        END IF

     END DO ! i

     col(j) = SUM( tblnz(:,j) )

  END DO ! j

  ! Finally, count nrownz and
  ! fill in icolnz, the translation table from symmetric to sparse
  ! TODO: Parallelize
  nrownz = 0
  DO j = 1, nrows
     IF( col(j) /= 0 ) THEN
        nrownz = nrownz + 1
        icolnz(nrownz) = j
     END IF

#if EBUG >= 3
     WRITE(*,*) 'zsy2sp_findnnz: iproc=', iproc, 'nrownz=', nrownz, ' icolnz=', icolnz(nrownz)
#endif /* DEBUG */

  END DO ! j

  RETURN
END SUBROUTINE zsy2sp_findnnz

!-------------------------------------------------------------------------------
! Permutes a sparse symmetric double complex matrix mat(nrows,nrows) such that
! only the first nrownz rows and columns are filled in matnz(nrownz,nrownz).
! The number of distinct non-zero areas is output as nareanz,
! and the start indices of each area is output to tblcolnz
SUBROUTINE zsy2sp_pack( nrows, mat, nrownz, icolnz, nareanz, tblcolnz, matnz )

  USE mod_prec, ONLY: dz
  USE mod_mpi_grid, ONLY: iproc

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: nrows, nrownz
  COMPLEX(KIND=dz), DIMENSION(nrows,nrows), INTENT(IN) :: mat
  INTEGER, DIMENSION(nrownz), INTENT(IN) :: icolnz
  INTEGER, INTENT(OUT) :: nareanz
  INTEGER, DIMENSION(0:nareanz), INTENT(OUT) :: tblcolnz
  COMPLEX(KIND=dz), DIMENSION(nrownz,nrownz), INTENT(OUT) :: matnz

  ! Internal variables
  INTEGER :: i, j, irow, icol, iarea, jarea
  LOGICAL :: toggle

  ! Count number of distinct non-zero areas
  ! Note: this part stays sequential
  nareanz = 1
  toggle = .FALSE.
  tblcolnz(0) = 1
  DO i = 1, nrownz

     IF( icolnz(i) /= i ) THEN
        ! Found a new area
        toggle = .TRUE.
     ELSE
        ! Inside an area, indices are consecutive
        irow = icolnz(i)
     END IF ! icolnz

#if EBUG >= 3
     WRITE(*,*) 'zsy2sp_pack: iproc=', iproc, ' nareanz=', nareanz, ' i=', i
#endif /* DEBUG */

     IF( toggle ) THEN
        ! Mark end of area
        tblcolnz(nareanz) = irow

        ! Move on to the next area
        nareanz = nareanz + 1
        toggle = .FALSE.
     END IF ! toggle

  END DO ! i

  ! Permute the non-zero data into matnz
  IF( nareanz > 1 ) THEN
     DO jarea = 1, nareanz
        DO icol = tblcolnz(jarea-1), tblcolnz(jarea)-1
           j = icolnz(icol)
           DO iarea = 1, nareanz
              DO irow = tblcolnz(iarea-1), tblcolnz(iarea)-1
                 i = icolnz(irow)
                 matnz(i,j) = mat(irow,icol)
              END DO ! irow
           END DO ! iarea
        END DO ! icol
     END DO ! jarea
  END IF ! nareanz

  RETURN
END SUBROUTINE zsy2sp_pack

!-------------------------------------------------------------------------------

END MODULE mod_sparse
