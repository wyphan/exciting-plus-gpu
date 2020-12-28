MODULE mod_sparse
  IMPLICIT NONE

CONTAINS

!===============================================================================
! Simple check for valid values in iperm_row(:) and iperm_col(:)
  subroutine check_iperm( nrow_small, ncol_small, nrow_big, ncol_big, &
                          iperm_row, iperm_col )

    implicit none

    ! Input arguments
    integer, intent(in) :: nrow_small, ncol_small
    integer, intent(in) :: nrow_big, ncol_big
    integer, intent(in) :: iperm_row(nrow_small)
    integer, intent(in) :: iperm_col(ncol_small)

    ! Internal variables
    integer :: nerrors
    integer :: irow_small, jcol_small
    integer :: irow_big, jcol_big
    logical :: isok

    nerrors = 0
    do irow_small=1,nrow_small
       irow_big = iperm_row(irow_small)
       isok = (1 <= irow_big).and.(irow_big <= nrow_big)
       if (.not.isok) then
          nerrors = nerrors + 1
          print *, 'irow_small, irow_big, nrow_big ', &
                   irow_small, irow_big, nrow_big
       endif
    enddo

    do jcol_small=1,ncol_small
       jcol_big = iperm_col(jcol_small)
       isok = (1 <= jcol_big).and.(jcol_big <= ncol_big)
       if (.not.isok) then
          nerrors = nerrors + 1
          print *, 'jcol_small, jcol_big, ncol_big ', &
                   jcol_small, jcol_big, ncol_big
       endif
    enddo

    if (nerrors.ne.0) then
       stop 'error in check_iperm '
    endif

    return
  end subroutine check_iperm

!===============================================================================
! Finds the nonzero rows and columns of a double complex matrix
! mat(1:nrows,1:ncols), and returns them as permutation vectors
! irownz(1:nrows) and icolnz(1:ncols), respectively.  
  SUBROUTINE zge2sp_findnnz( nrows, ncols, mat, lda, &
                             nrownz, irownz, ncolnz, icolnz )

    USE mod_prec, ONLY: dd, dz
    USE mod_mpi_grid, ONLY: iproc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: nrows, ncols, lda
    COMPLEX(KIND=dz), DIMENSION(lda,ncols), INTENT(IN) :: mat
    INTEGER, INTENT(OUT) :: nrownz, ncolnz
    INTEGER, DIMENSION(nrows), INTENT(OUT) :: irownz
    INTEGER, DIMENSION(ncols), INTENT(OUT) :: icolnz

    ! Internal variables
    REAL(KIND=dd) :: tol = 1.D-10
    LOGICAL, DIMENSION(nrows,ncols) :: tblnz
    INTEGER, DIMENSION(nrows) :: col
    INTEGER, DIMENSION(ncols) :: row
    INTEGER :: i, j

    ! Initialize vars
    nrownz = 0
    ncolnz = 0
    irownz(:) = 0
    icolnz(:) = 0
    tblnz(:,:) = .FALSE.
    col(:) = 0
    row(:) = 0

    ! Fill in tblnz
    ! TODO: Parallelize
    DO j = 1, ncols
       DO i = 1, nrows
          IF( ABS(mat(i,j)) >= tol ) tblnz(i,j) = .TRUE.
       END DO
       col(j) = COUNT( tblnz(:,j) )
    END DO
    DO i = 1, nrows
       row(i) = COUNT( tblnz(i,:) )
    END DO

    ! Count nrownz and ncolnz
    DO i = 1, nrows
       IF( row(i) /= 0 ) nrownz = nrownz + 1
    END DO
    DO j = 1, ncols
       IF( col(j) /= 0 ) ncolnz = ncolnz + 1
    END DO

#if EBUG >= 3
    WRITE(*,*) 'zsy2sp_findnnz: iproc=', iproc, 'nrownz=', nrownz, ' ncolnz=', ncolnz
#endif /* DEBUG */

    ! Fill in irownz
    ! TODO: Parallelize
    DO i = 1, nrownz
       IF( row(i) /= 0 ) irownz(i) = i

#if EBUG >= 3
       WRITE(*,*) 'zsy2sp_findnnz: iproc=', iproc, ' irownz=', irownz(i)
#endif /* DEBUG */

    END DO

    ! Fill in icolnz
    ! TODO: Parallelize
    DO j = 1, ncolnz
       IF( col(j) /= 0 ) icolnz(j) = j

#if EBUG >= 3
       WRITE(*,*) 'zsy2sp_findnnz: iproc=', iproc, ' icolnz=', icolnz(j)
#endif /* DEBUG */

    END DO

    RETURN
  END SUBROUTINE zge2sp_findnnz

!===============================================================================
! Packs a sparse matrix mat(1:nrows,1:ncols) into matnz(1:nrownz,1:ncolnz)
! Call zge2sp_findnnz() first before calling this subroutine!
  SUBROUTINE zge2sp_pack( nrows, ncols, mat, lda, &
                          nrownz, irownz, ncolnz, icolnz, matnz, ldanz )

    USE mod_prec, ONLY: dd, dz
    USE mod_mpi_grid, ONLY: iproc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: nrows, ncols, lda, nrownz, ncolnz, ldanz
    INTEGER, DIMENSION(nrows), INTENT(IN) :: irownz
    INTEGER, DIMENSION(ncols), INTENT(IN) :: icolnz
    COMPLEX(KIND=dz), DIMENSION(lda,ncols), INTENT(IN) :: mat
    COMPLEX(KIND=dz), DIMENSION(ldanz,ncolnz), INTENT(OUT) :: matnz

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
! Unpacks a sparse matrix matnz(1:nrownz,1:ncolnz) into mat(1:nrows,1:ncols)
  SUBROUTINE zsp2ge_unpack( nrownz, irownz, ncolnz, icolnz, matnz, ldanz, &
                            nrows, ncols, mat, lda )
  
    USE mod_prec, ONLY: dd, dz
    USE mod_mpi_grid, ONLY: iproc

    IMPLICIT NONE

    ! Arguments
    INTEGER, INTENT(IN) :: nrows, ncols, lda, nrownz, ncolnz, ldanz
    INTEGER, DIMENSION(nrows), INTENT(IN) :: irownz
    INTEGER, DIMENSION(ncols), INTENT(IN) :: icolnz
    COMPLEX(KIND=dz), DIMENSION(ldanz,ncolnz), INTENT(IN) :: matnz
    COMPLEX(KIND=dz), DIMENSION(lda,ncols), INTENT(OUT) :: mat

    ! Internal variables
    INTEGER, PARAMETER :: idebug = 0
    INTEGER :: irow_small, jcol_small, irow_big, jcol_big
    COMPLEX(KIND=dz) :: aij

    if (idebug >= 1) then
       call check_iperm( nrownz, ncolnz, nrows, ncols, &
                         irownz, icolnz )
    endif

    ! Set default value to zero
    do jcol_big = 1, ncols
       do irow_big = 1, nrows
          mat(irow_big,jcol_big) = 0
       enddo
    enddo

    ! Set non-zero values
    ! equivalent to
    ! A_big(iperm_row(:),iperm_col(:)) = A_small(:,:)
    do jcol_small = 1, ncolnz
       do irow_small = 1, nrownz

          irow_big = irownz(irow_small)
          jcol_big = icolnz(jcol_small)

          aij = matnz(irow_small,jcol_small)
          mat(irow_big,jcol_big) = aij

       enddo
    enddo

    RETURN
  END SUBROUTINE zsp2ge_unpack

!===============================================================================

END MODULE mod_sparse
