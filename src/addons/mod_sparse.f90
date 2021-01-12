MODULE mod_sparse

  USE mod_prec, ONLY: dd

  IMPLICIT NONE

  ! Packing tolerance
  REAL(KIND=dd) :: packtol ! Default is 1.D-10, set in readinput()

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
! Helper function to check uplo for zsy2sp_* subroutines
  LOGICAL FUNCTION check_uplo( uplo, func )
    IMPLICIT NONE

    ! Input argument
    CHARACTER(LEN=1), INTENT(IN) :: uplo
    CHARACTER(LEN=*), INTENT(IN) :: func

    SELECT CASE( uplo )
    CASE('U')
       check_uplo = .TRUE.
    CASE('L')
       check_uplo = .FALSE.
    CASE DEFAULT
       WRITE(*,*) 'Error(', TRIM(func), '): invalid uplo argument ', uplo
       STOP 'Error in check_uplo'
    END SELECT

    RETURN
  END FUNCTION check_uplo

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
    REAL(KIND=dd), DIMENSION(nrows) :: rownorm
    REAL(KIND=dd), DIMENSION(ncols) :: colnorm
    REAL(KIND=dd) :: aij, maxrownorm, maxcolnorm
    LOGICAL, DIMENSION(nrows) :: keeprow
    LOGICAL, DIMENSION(ncols) :: keepcol
    LOGICAL :: lnz
    INTEGER :: i, j, nrowsmall, ncolsmall

    ! Initialize vars
    nrowsmall = 0
    ncolsmall = 0
    ! irownz(:) = 0 ! for some reason Valgrind flags these two lines as
    ! icolnz(:) = 0 ! out-of-bounds memory access... disabling for now
    keeprow(:) = .FALSE.
    keepcol(:) = .FALSE.
    rownorm(:) = 0._dd
    colnorm(:) = 0._dd

    ! Find nonzeroes using absolute norm magnitudes
    DO j = 1, ncols
       DO i = 1, nrows
          lnz = ( ABS(mat(i,j)) >= packtol )
          IF( lnz ) THEN
             keeprow(i) = .TRUE.
             keepcol(j) = .TRUE.
          END IF
       END DO ! i
    END DO ! j

    ! Find nonzeroes using relative norm magnitudes
    !DO j = 1, ncols
    !   DO i = 1, nrows
    !      aij = ABS(mat(i,j))
    !      rownorm(i) = rownorm(i) + aij
    !      colnorm(j) = colnorm(j) + aij
    !   END DO ! i
    !END DO ! j
    !maxrownorm = MAXVAL(rownorm)
    !maxcolnorm = MAXVAL(colnorm)
    !keeprow(1:nrows) = ( rownorm(1:nrows) > packtol*maxrownorm )
    !keepcol(1:ncols) = ( colnorm(1:ncols) > packtol*maxcolnorm )

!#if EBUG >= 3
    !WRITE(*,*) 'zge2sp_findnnz: iproc=', iproc, &
    !           ' maxrownorm=', maxrownorm, ' maxcolnorm=', maxcolnorm
!#endif /* DEBUG */

    ! Count nrownz and fill in irownz
    DO i = 1, nrows
       IF( keeprow(i) ) THEN
          nrowsmall = nrowsmall + 1
          irownz(nrowsmall) = i

#if EBUG >= 3
          WRITE(*,*) 'zge2sp_findnnz: iproc=', iproc, ' irownz(', nrowsmall, ')=', irownz(nrowsmall)
#endif /* DEBUG */

       END IF
    END DO ! i
    nrownz = nrowsmall
    
    ! Count ncolnz and fill in icolnz
    DO j = 1, ncols
       IF( keepcol(j) ) THEN
          ncolsmall = ncolsmall + 1
          icolnz(ncolsmall) = j

#if EBUG >= 3
          WRITE(*,*) 'zge2sp_findnnz: iproc=', iproc, ' icolnz(', ncolsmall, ')=', icolnz(ncolsmall)
#endif /* DEBUG */

       END IF
    END DO ! j
    ncolnz = ncolsmall


#if EBUG >= 3
    WRITE(*,*) 'zge2sp_findnnz: iproc=', iproc, 'nrownz=', nrowsmall, ' ncolnz=', ncolsmall
#endif /* DEBUG */

    RETURN
  END SUBROUTINE zge2sp_findnnz

!===============================================================================
! Finds the nonzero rows of a double complex symmetric matrix
! mat(1:nrows,1:ncols), and returns it as a permutation vector
! irownz(1:nrows).
  SUBROUTINE zsy2sp_findnnz( uplo, nrows, mat, lda, &
                             nrownz, irownz )

    USE mod_prec, ONLY: dd, dz
    USE mod_mpi_grid, ONLY: iproc

    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: uplo
    INTEGER, INTENT(IN) :: nrows, lda
    COMPLEX(KIND=dz), DIMENSION(lda,nrows), INTENT(IN) :: mat
    INTEGER, INTENT(OUT) :: nrownz
    INTEGER, DIMENSION(nrows), INTENT(OUT) :: irownz

    ! Internal variables
    REAL(KIND=dd), DIMENSION(nrows) :: rownorm
    REAL(KIND=dd) :: aii, aij, maxrownorm
    LOGICAL, DIMENSION(nrows) :: keeprow
    LOGICAL :: lup, lnz
    INTEGER :: i, j, rowstart, rowend, nrowsmall

    ! Check uplo
    lup = check_uplo(uplo,'zsy2sp_findnnz')

    ! Initialize vars
    nrowsmall = 0
    ! irownz(:) = 0 ! See comment in zge2sp_findnnz() above
    keeprow(:) = .FALSE.
    rownorm(:) = 0._dd

    ! Initialize rownorm with diagonal values
    DO i = 1, nrows
       aii = ABS(mat(i,i))
       rownorm(i) = aii
    END DO ! i

    ! Find nonzeroes using absolute norm magnitudes
    DO j = 1, nrows ! technically ncols
       rowstart = MERGE( 1,   j+1,   lup )
       rowend   = MERGE( j-1, nrows, lup )
       DO i = rowstart, rowend
          lnz = ( ABS(mat(i,j)) >= packtol )
          IF( lnz ) THEN
             keeprow(i) = .TRUE.
             keeprow(j) = .TRUE.
          END IF ! lnz
       END DO ! i
    END DO ! j
    
    ! Find nonzeroes using relative norm magnitudes
    !DO j = 1, nrows ! technically ncols
    !   rowstart = MERGE( 1,   j+1,   lup )
    !   rowend   = MERGE( j-1, nrows, lup )
    !   DO i = rowstart, rowend
    !      aij = ABS(mat(i,j))
    !      rownorm(i) = rownorm(i) + aij
    !      rownorm(j) = rownorm(j) + aij
    !   END DO ! i
    !END DO ! j             
    !maxrownorm = MAXVAL(rownorm)
    !keeprow(1:nrows) = ( rownorm(1:nrows) > packtol*maxrownorm )

!#if EBUG >= 3
    !WRITE(*,*) 'zsy2sp_findnnz: iproc=', iproc, &
    !           ' maxrownorm=', maxrownorm
!#endif /* DEBUG */

    ! Count nrownz and fill in irownz
    DO i = 1, nrows
       IF( keeprow(i) ) THEN
          nrowsmall = nrowsmall + 1
          irownz(nrowsmall) = i

#if EBUG >= 3
          WRITE(*,*) 'zsy2sp_findnnz: iproc=', iproc, ' irownz(', nrowsmall, ')=', irownz(nrowsmall)
#endif /* DEBUG */

       END IF
    END DO ! i
    nrownz = nrowsmall
    
#if EBUG >= 3
    WRITE(*,*) 'zsy2sp_findnnz: iproc=', iproc, 'nrownz=', nrowsmall
#endif /* DEBUG */

    RETURN
  END SUBROUTINE zsy2sp_findnnz

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
          IF( i /= 0 .AND. j /= 0 ) matnz(irow,icol) = mat(i,j)
       END DO ! irow
    END DO ! icol

    RETURN
  END SUBROUTINE zge2sp_pack

!===============================================================================
! Packs a symmetric sparse matrix mat(1:nrows,1:nrows) into
! matnz(1:nrownz,1:nrownz)
! Call zsy2sp_findnnz() first before calling this subroutine!
  SUBROUTINE zsy2sp_pack( uplo, nrows, mat, lda, &
                          nrownz, irownz, matnz, ldanz )

    USE mod_prec, ONLY: dd, dz
    USE mod_mpi_grid, ONLY: iproc

    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: uplo
    INTEGER, INTENT(IN) :: nrows, lda, nrownz, ldanz
    INTEGER, DIMENSION(nrows), INTENT(IN) :: irownz
    COMPLEX(KIND=dz), DIMENSION(lda,nrows), INTENT(IN) :: mat
    COMPLEX(KIND=dz), DIMENSION(ldanz,nrownz), INTENT(OUT) :: matnz

    ! Internal variables
    INTEGER :: i, j, irow, icol, rowstart, rowend
    LOGICAL :: lup

    ! Check uplo
    lup = check_uplo(uplo,'zsy2sp_pack')

    ! Permute the non-zero data into matnz
    ! TODO: Parallelize
    DO icol = 1, nrownz ! technically ncolnz
       j = irownz(icol)
       rowstart = MERGE( 1,      icol+1, lup )
       rowend   = MERGE( icol-1, nrownz, lup )
       DO irow = rowstart, rowend
          i = irownz(irow)
          IF( i /= 0 .AND. j /= 0 ) matnz(irow,icol) = mat(i,j)
       END DO ! irow
    END DO ! icol

    RETURN
  END SUBROUTINE zsy2sp_pack

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
    INTEGER, DIMENSION(nrows) :: map_row
    INTEGER, DIMENSION(ncols) :: map_col
    LOGICAL :: is_nonzero
    COMPLEX(KIND=dz) :: aij

    if (idebug >= 1) then
       call check_iperm( nrownz, ncolnz, nrows, ncols, &
                         irownz, icolnz )
    endif

    map_row(:) = 0
    map_col(:) = 0

    do irow_small = 1, nrownz
       irow_big = irownz(irow_small)
       map_row(irow_big) = irow_small
    enddo
    
    do jcol_small = 1, ncolnz
       jcol_big = icolnz(jcol_small)
       map_col(jcol_big) = jcol_small
    enddo

    
    ! Initialize A_big(:,:) in a single pass
    ! equivalent to
    ! A_big(iperm_row(:),iperm_col(:)) = A_small(:,:)
    do jcol_big = 1, ncols

       jcol_small = map_col(jcol_big)
       
       do irow_big = 1, nrows

          irow_small = map_row(irow_big)
          aij = (0._dd,0._dd)

          is_nonzero = ( irow_small >= 1 ) .and. ( jcol_small >= 1 )
          if( is_nonzero ) aij = matnz(irow_small,jcol_small)

          mat(irow_big,jcol_big) = aij

       enddo
    enddo

    RETURN
  END SUBROUTINE zsp2ge_unpack

!===============================================================================
! Unpacks a symmetric sparse matrix matnz(1:nrownz,1:nrownz) into
! mat(1:nrows,1:nrows)
  SUBROUTINE zsp2sy_unpack( uplo, nrownz, irownz, matnz, ldanz, &
                            nrows, mat, lda )
  
    USE mod_prec, ONLY: dd, dz
    USE mod_mpi_grid, ONLY: iproc

    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: uplo
    INTEGER, INTENT(IN) :: nrows, lda, nrownz, ldanz
    INTEGER, DIMENSION(nrows), INTENT(IN) :: irownz
    COMPLEX(KIND=dz), DIMENSION(ldanz,nrownz), INTENT(IN) :: matnz
    COMPLEX(KIND=dz), DIMENSION(lda,nrows), INTENT(OUT) :: mat

    ! Internal variables
    INTEGER, PARAMETER :: idebug = 0
    INTEGER :: irow_small, jcol_small, irow_big, jcol_big, rowstart, rowend
    INTEGER, DIMENSION(nrows) :: map_row
    LOGICAL :: lup, is_nonzero
    COMPLEX(KIND=dz) :: aij

    ! Check uplo
    lup = check_uplo(uplo,'zsp2sy_unpack')

    ! Technically this checks the same thing twice...
    if (idebug >= 1) then
       call check_iperm( nrownz, nrownz, nrows, nrows, &
                         irownz, irownz )
    endif

    map_row(:) = 0

    do irow_small = 1, nrownz
       irow_big = irownz(irow_small)
       map_row(irow_big) = irow_small
    enddo
    
    ! Initialize A_big(:,:) in a single pass
    ! equivalent to
    ! A_big(iperm_row(:),iperm_row(:)) = A_small(:,:)
    do jcol_big = 1, nrows

       jcol_small = map_row(jcol_big)

       rowstart = MERGE( 1,          jcol_big+1,   lup )
       rowend   = MERGE( jcol_big-1, nrows,        lup )
       DO irow_big = rowstart, rowend

          irow_small = map_row(irow_big)
          aij = (0._dd,0._dd)

          is_nonzero = ( irow_small >= 1 ) .and. ( jcol_small >= 1 )
          if( is_nonzero ) aij = matnz(irow_small,jcol_small)

          mat(irow_big,jcol_big) = aij

       enddo ! irow_big
    enddo ! jcol_big

    RETURN
  END SUBROUTINE zsp2sy_unpack

!===============================================================================

END MODULE mod_sparse
