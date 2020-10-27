MODULE mod_sparse

  USE mod_prec, ONLY: dz
  IMPLICIT NONE

  COMPLEX(KIND=dz), ALLOCATABLE :: gntuju_packed(:,:,:,:)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nmt
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: nareanz, igntujunz, tblgntujunz

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
  icolnz(:) = 0
  DO j = 1, nrows
     IF( col(j) /= 0 ) THEN
        nrownz = nrownz + 1
        icolnz(nrownz) = j
     END IF
  END DO

#if EBUG >= 3
  WRITE(*,*) 'zsy2sp_findnnz: iproc=', iproc, 'nrownz=', nrownz, ' icolnz=', icolnz
#endif /* DEBUG */

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
  irow = 0
  toggle = .FALSE.
  tblcolnz(0) = 1
  DO i = 1, nrownz

     IF( icolnz(i) /= irow ) THEN

        toggle = .TRUE.
        tblcolnz(nareanz) = icolnz(i)
        IF( toggle ) THEN

#if EBUG >= 3
           WRITE(*,*) 'zsy2sp_pack: iproc=', iproc, ' nareanz=', nareanz, ' irow=', irow
#endif /* DEBUG */

           irow = icolnz(i)
           nareanz = nareanz + 1
           toggle = .FALSE.

        END IF ! toggle

     ELSE

        toggle = .FALSE.
        irow = irow + 1

     END IF ! icolnz
  END DO ! i

  ! Permute the data
  IF( nareanz > 1 ) THEN
     DO iarea = 1, nareanz
        DO jarea = 1, nareanz

           DO icol = tblcolnz(jarea-1), tblcolnz(jarea)-1

              j = icolnz(icol)
              DO irow = tblcolnz(iarea-1), tblcolnz(iarea)-1

                 i = icolnz(irow)
                 matnz(irow,icol) = mat(i,j)

              END DO ! irow
           END DO ! icol
        END DO ! jarea
     END DO ! iarea
  END IF ! nareanz

  RETURN
END SUBROUTINE zsy2sp_pack

!-------------------------------------------------------------------------------

END MODULE mod_sparse
