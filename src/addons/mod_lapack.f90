MODULE mod_lapack

  USE mod_prec

  IMPLICIT NONE

  INTERFACE

     ! BLAS subroutine
     ! Double precision matrix-matrix multiply C = alpha * A x B + beta * C
     SUBROUTINE DGEMM( transA, transB, m, n, k, &
                       alpha, matA, lda, &
                              matB, ldb, &
                       beta,  matC, ldc )
       USE mod_prec, ONLY: dd
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: transA, transB
       INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc
       REAL(KIND=dd), INTENT(IN) :: alpha, beta
       REAL(KIND=dd), INTENT(IN), DIMENSION(lda,*) :: matA
       REAL(KIND=dd), INTENT(IN), DIMENSION(ldb,*) :: matB
       REAL(KIND=dd), INTENT(INOUT), DIMENSION(ldc,*) :: matC
     END SUBROUTINE DGEMM

     ! BLAS subroutine
     ! Double complex matrix-matrix multiply C = alpha * A x B + beta * C
     SUBROUTINE ZGEMM( transA, transB, m, n, k, &
                       alpha, matA, lda, &
                              matB, ldb, &
                       beta,  matC, ldc )
       USE mod_prec, ONLY: dz
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: transA, transB
       INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc
       COMPLEX(KIND=dz), INTENT(IN) :: alpha, beta
       COMPLEX(KIND=dz), INTENT(IN), DIMENSION(lda,*) :: matA
       COMPLEX(KIND=dz), INTENT(IN), DIMENSION(ldb,*) :: matB
       COMPLEX(KIND=dz), INTENT(INOUT), DIMENSION(ldc,*) :: matC
     END SUBROUTINE ZGEMM

     ! LAPACK subroutine
     ! Double precision LU decomposition
     SUBROUTINE DGETRF( m, n, matA, lda, ipiv, info )
       USE mod_prec, ONLY: dd
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, n, lda
       INTEGER, INTENT(OUT) :: info
       INTEGER, DIMENSION(*), INTENT(OUT) :: ipiv
       REAL(KIND=dd), INTENT(INOUT), DIMENSION(lda,*) :: matA
     END SUBROUTINE DGETRF

     ! LAPACK subroutine
     ! Double precision matrix inversion on LU-decomposed matrix
     SUBROUTINE DGETRI( n, matA, lda, ipiv, work, lwork, info )
       USE mod_prec, ONLY: dd
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, lda, lwork
       INTEGER, INTENT(OUT) :: info
       INTEGER, DIMENSION(n), INTENT(IN) :: ipiv
       REAL(KIND=dd), INTENT(INOUT), DIMENSION(lda,*) :: matA
       REAL(KIND=dd), INTENT(OUT), DIMENSION(*) :: work
     END SUBROUTINE DGETRI

     ! LAPACK subroutine
     ! Double complex LU decomposition
     SUBROUTINE ZGETRF( m, n, matA, lda, ipiv, info )
       USE mod_prec, ONLY: dz
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: m, n, lda
       INTEGER, INTENT(OUT) :: info
       INTEGER, DIMENSION(*), INTENT(OUT) :: ipiv
       COMPLEX(KIND=dz), INTENT(INOUT), DIMENSION(lda,*) :: matA
     END SUBROUTINE ZGETRF

     ! LAPACK subroutine
     ! Double complex matrix inversion on LU-decomposed matrix
     SUBROUTINE ZGETRI( n, matA, lda, ipiv, work, lwork, info )
       USE mod_prec, ONLY: dz
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: n, lda, lwork
       INTEGER, INTENT(OUT) :: info
       INTEGER, DIMENSION(n), INTENT(IN) :: ipiv
       COMPLEX(KIND=dz), INTENT(INOUT), DIMENSION(lda,*) :: matA
       COMPLEX(KIND=dz), INTENT(OUT), DIMENSION(*) :: work
     END SUBROUTINE ZGETRI

  END INTERFACE

END MODULE mod_lapack

!===============================================================================
! TODO: Include this into the module

  subroutine invzge(mtrx,ndim)
    USE mod_prec
    USE mod_lapack ! TODO: remove this once included in module
    implicit none

    ! passed var
    integer, intent(in) :: ndim
    COMPLEX(KIND=dz), INTENT(INOUT) :: mtrx(ndim,ndim)

    ! local var
    integer lwork,nb,info
    COMPLEX(KIND=dz), ALLOCATABLE :: work(:)

!integer, external :: ilaenv
!nb=ilaenv(1,'zgetri','U',ndim,-1,-1,-1)
!lwork=ndim*nb

!--begin IBM ESSL fix
INTEGER :: ipiv(ndim)
COMPLEX(KIND=dz), DIMENSION(1) :: query

! Query workspace
CALL zgetri( ndim, mtrx, ndim, ipiv, query, -1, info)
lwork = CEILING(REAL(query(1),KIND=dd))
!--end IBM ESSL fix

allocate(work(lwork))
call zgetrf(ndim,ndim,mtrx,ndim,ipiv,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(invzge) zgetrf returned ",I4)')info
  call pstop
endif
call zgetri(ndim,mtrx,ndim,ipiv,work,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(invzge): zgetri returned ",I4)')info
  call pstop
endif
deallocate(work)
end subroutine invzge

!===============================================================================
! TODO: Include this into the module

subroutine invdsy(n,mtrx)
  USE mod_lapack ! TODO: remove this once included in module
  implicit none
integer, intent(in) :: n
real(8), intent(inout) :: mtrx(n,n)

real(8) t1
integer lwork,info
real(8), allocatable :: work(:)

!--begin IBM ESSL fix
INTEGER :: ipiv(n)
REAL(KIND=dd), DIMENSION(1) :: query
!--end IBM ESSL fix

!allocate(ipiv(n))
!lwork=-1
!call dsytrf('U',n,mtrx,n,ipiv,t1,lwork,info)
!lwork=int(t1)+1
!allocate(work(lwork))
!call dsytrf('U',n,mtrx,n,ipiv,work,lwork,info)
!if (info.ne.0) then
!  write(*,*)
!  write(*,'("Error(invdsy) dsytrf returned ",I4)')info
!  call pstop
!endif
!call dsytri('U',n,mtrx,n,ipiv,work,info)
!if (info.ne.0) then
!  write(*,*)
!  write(*,'("Error(invdsy) dsytri returned ",I4)')info
!  call pstop
!endif
!deallocate(ipiv,work)

!--begin IBM ESSL fix

! For some reason, ESSL doesn't provide *sytri
! Use *getrf and *getri instead

! Query workspace
CALL DGETRI( n, mtrx, n, ipiv, query, -1, info )
lwork = CEILING(query(1))

! LU factorization using DGETRF
CALL DGETRF( n, n, mtrx, n, ipiv, info )
IF ( info /= 0 ) THEN
  WRITE(*,*)
  WRITE(*,'("Error(invdsy) dgetrf returned ",I4)') info
  CALL pstop
END IF

! Allocate workspace
ALLOCATE( work(lwork) )

! Invert matrix using DGETRI
CALL DGETRI( n, mtrx, n, ipiv, work, lwork, info )
! Regardless of results, ipiv and work are no longer needed at this point
DEALLOCATE( work )
IF ( info /= 0 ) THEN
  WRITE(*,*)
  WRITE(*,'("Error(invdsy) dgetri returned ",I4)') info
  CALL pstop
END IF

!--end IBM ESSL fix
return
end subroutine invdsy

!===============================================================================
! TODO: Write interfaces to LAPACK for the subroutines called from the following

subroutine diagzhe(ndim,mtrx,evalue)
implicit none
integer ,intent(in) :: ndim
complex(8), intent(out) :: mtrx(ndim,ndim)
real(8), intent(out) :: evalue(ndim)
integer nb,lwork,inf, dummy
real*8, allocatable :: work(:),rwork(:)

!integer, external :: ilaenv
!nb=ilaenv(1,'ZHETRD','U',ndim,-1,-1,-1)
!lwork=(nb+1)*ndim

!--begin IBM ESSL fix
INTEGER :: dummy2
REAL(KIND=KIND(1.D0)) :: query, dummy1

! Query workspace
CALL zheev( 'V', 'U', ndim, mtrx, ndim, dummy1, query, -1, dummy2, inf )
lwork = CEILING(query)
!--end IBM ESSL fix

allocate(work(lwork*2))
allocate(rwork(3*ndim+2))
call zheev('V','U',ndim,mtrx,ndim,evalue,work,lwork,rwork,inf)
if (inf.ne.0) then
  write(*,*)
  write(*,'("Error(diagzhe) zheev returned ",I4)')inf
  call pstop
endif
deallocate(work,rwork)
end

subroutine diagdsy(n,mtrx,eval)
implicit none
integer, intent(in) :: n
real(8), intent(inout) :: mtrx(n,n)
real(8), intent(out) :: eval(n)
integer lwork,info
real(8), allocatable :: work(:)
real(8) t1
lwork=-1
call dsyev('V','U',n,mtrx,n,eval,t1,lwork,info)
lwork=int(t1)+1
allocate(work(lwork))
call dsyev('V','U',n,mtrx,n,eval,work,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Warning(diagdsy) dsyev returned ",I4)')info
  write(*,*)
endif
deallocate(work)
return
end subroutine diagdsy

subroutine isqrtzhe(ndim,mtrx,ierr)
use mod_mpi_grid
implicit   none
! arguments
integer, intent(in) :: ndim
integer, intent(out) :: ierr
complex(8), intent(inout) :: mtrx(ndim,ndim)

integer nb,lwork,info,i,j,n, dummy
real(8), allocatable :: work(:),rwork(:),ev(:),ev1(:)
complex(8), allocatable :: z1(:,:)

!integer, external :: ilaenv
!nb=ilaenv(1,'ZHETRD','U',ndim,-1,-1,-1)
!lwork=(nb+1)*ndim

!--begin IBM ESSL fix
INTEGER :: dummy2
REAL(KIND=KIND(1.D0)) :: query, dummy1

! Query workspace
CALL zheev( 'V', 'U', ndim, mtrx, ndim, dummy1, query, -1, dummy2, info )
lwork = CEILING(query)
!--end IBM ESSL fix

allocate(z1(ndim,ndim))
allocate(work(lwork*2))
allocate(rwork(3*ndim+2))
allocate(ev(ndim))
allocate(ev1(ndim))

!if (mpi_grid_root()) write(*,*)mtrx

call zheev('V','U',ndim,mtrx,ndim,ev1,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(isqrtzhe): zheev returned ",I4)')info
  write(*,*)
  call pstop
endif

!if (mpi_grid_root()) write(*,*)ev1
!if (mpi_grid_root()) write(*,*)mtrx

ierr=0
do i=1,ndim
  if (ev1(i).lt.1.d-12) then
    ierr=i
    ev(i)=0.d0
  else
    ev(i)=1.d0/dsqrt(ev1(i))
  endif
enddo

z1(:,:)=mtrx(:,:)
mtrx=dcmplx(0.d0,0.d0)
do i=1,ndim
  do j=1,ndim
    do n=1,ndim
      mtrx(i,j)=mtrx(i,j)+z1(i,n)*dconjg(z1(j,n))*ev(n) 
    enddo
  enddo
enddo
deallocate(z1,work,rwork,ev,ev1)
return
end subroutine isqrtzhe

subroutine isqrtdsy(n,mtrx,ierr)
use mod_mpi_grid
implicit   none
! arguments
integer, intent(in) :: n
integer, intent(out) :: ierr
real(8), intent(inout) :: mtrx(n,n)
!
integer i,j,k
real(8), allocatable :: eval(:),eval_isq(:)
real(8), allocatable :: evec(:,:)

allocate(evec(n,n))
allocate(eval(n),eval_isq(n))
call diagdsy(n,mtrx,eval)

ierr=0
do i=1,n
  if (eval(i).lt.1.d-12) then
    ierr=i
    eval_isq(i)=0.d0
  else
    eval_isq(i)=1.d0/dsqrt(eval(i))
  endif
enddo
evec(:,:)=mtrx(:,:)
mtrx=0.d0
do i=1,n
  do j=1,n
    do k=1,n
      mtrx(i,j)=mtrx(i,j)+evec(i,k)*evec(j,k)*eval_isq(k) 
    enddo
  enddo
enddo
deallocate(evec,eval,eval_isq)
return
end subroutine isqrtdsy

subroutine diagzge(ndim,mtrx,evalue)
implicit   none
integer, intent(in) :: ndim
complex*16, intent(inout) :: mtrx(ndim,ndim)
complex*16, intent(out) :: evalue(ndim)
integer lwork,inf
complex(8) zt1
real*8, allocatable :: rwork(:)
complex(8), allocatable :: work(:)
complex(8), allocatable :: evec(:,:)

!--begin IBM ESSL fix
!integer, external :: ilaenv
!--end IBM ESSL fix

lwork=-1
call zgeev('N','V',ndim,mtrx,ndim,evalue,evec,ndim,evec,ndim,zt1,lwork,rwork,inf)
lwork=dble(zt1)+1
allocate(work(lwork))
allocate(rwork(2*ndim))
allocate(evec(ndim,ndim))
call zgeev('N','V',ndim,mtrx,ndim,evalue,evec,ndim,evec,ndim,work,lwork,rwork,inf)
if (inf.ne.0) then
  write(*,'("Error(diagzge) zgeev returned ",I4)')inf
  call pstop
endif
mtrx=evec
deallocate(work,rwork,evec)
end subroutine diagzge

subroutine diagzheg(n,nv,ld,etol,a,b,eval,evec)
implicit none
integer, intent(in) :: n
integer, intent(in) :: nv
integer, intent(in) :: ld
real(8), intent(in) :: etol
complex(8), intent(inout) :: a(n,n)
complex(8), intent(inout) :: b(n,n)
real(8), intent(out) :: eval(nv)
complex(8), intent(out) :: evec(ld,nv)
!
integer m,info,i,nb,lwork
real(8) vl,vu
integer, allocatable :: iwork(:)
integer, allocatable :: ifail(:)
real(8), allocatable :: w(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)

!integer, external :: ilaenv 
!nb=ilaenv(1,'ZHETRD','U',n,-1,-1,-1)
!lwork=(nb+1)*n

!--begin IBM ESSL fix
INTEGER :: dummy1, dummy2
REAL(KIND=KIND(1.D0)) :: query

! Query workspace
CALL zhegvx( 1, 'V', 'I', 'U', n, a, n, b, n, vl, vu, 1, nv, etol, m, w, evec,&
             ld, query, -1, dummy1, dummy2, ifail, info )
lwork = CEILING(query)
!--end IBM ESSL fix

allocate(iwork(5*n))
allocate(ifail(n))
allocate(w(n))
allocate(rwork(7*n))
allocate(work(lwork))
call zhegvx(1,'V','I','U',n,a,n,b,n,vl,vu,1,nv,etol,m,w,evec,ld,&
  work,lwork,rwork,iwork,ifail,info)
eval(1:nv)=w(1:nv)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(diagzheg): diagonalisation failed")')
  write(*,'("  zhegvx returned info = ",I8)') info
  if (info.gt.n) then
    i=info-n
    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
    write(*,'("  is not positive definite")')
    write(*,'(" Order of overlap matrix : ",I8)') n
    write(*,*)
  end if
  call pstop
end if
deallocate(iwork,ifail,w,rwork,work)
return
end subroutine diagzheg

