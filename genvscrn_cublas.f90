SUBROUTINE genvscrn( iq, chi0, krnl, vscrn, epsilon )

  USE modmain
  USE mod_addons_q
  USE mod_expigqr
  USE mod_linresp

  USE ISO_C_BINDING
  USE cublas_f

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: iq
  COMPLEX(8), INTENT(IN), TARGET :: chi0( ngq(iq), ngq(iq) )
  COMPLEX(8), INTENT(IN), TARGET :: krnl( ngq(iq), ngq(iq) )
  COMPLEX(8), INTENT(OUT), TARGET :: vscrn( ngq(iq), ngq(iq) )
  COMPLEX(8), INTENT(OUT), TARGET :: epsilon( ngq(iq), ngq(iq) )

  ! Local variables
  INTEGER :: ig, ig1, ig2, stage

  ! CUDA stuff
  INTEGER(C_INT) :: stat, lwork
  TYPE(C_PTR) :: handleblas, handlesolv
  INTEGER(C_INT), DIMENSION(:), TARGET :: ipiv
  INTEGER(C_INT), TARGET :: info
  TYPE(C_PTR) :: d_chi0, d_krnl, d_epsilon, d_work, d_ipiv, d_info
  COMPLEX(8), DIMENSION(:), POINTER :: h_chi0, h_krnl, h_vscrn, h_epsilon
  INTEGER(C_INT), DIMENSION(:), POINTER :: h_ipiv
  INTEGER(C_INT), POINTER :: h_info => info
  TYPE(C_PTR), POINTER :: offset, address, address2
  INTEGER(8) :: n, sz, bytes

  CALL papi_timer_start( pt_vscrn )

  ! Both cublas and cusolver has been initialized on main
  !stat = cublasCreate( handleblas )
  !stat = cublasSetStream( handleblas, stream )
  !stat = cusolverDnCreate( handlesolv )
  !stat = cusolverDnSetStream( handlesolv, stream )

  ! Matrix sizes
  n = ngq(iq)
  sz = n**2
  bytes = szmat * sizeof_complex
  WRITE(*,'("[genvscrn] Matrix sizes are ",F7.2," MB each")') DBLE(bytes)*toMB

  ! Set up host pointers
  ALLOCATE( h_chi0, sz )
  ALLOCATE( h_krnl, sz )
  ALLOCATE( h_vscrn, sz )
  ALLOCATE( h_epsilon, sz )
  ALLOCATE( ipiv, n*sizeof_int )
  h_chi0 => chi0
  h_krnl => krnl
  h_vscrn => vscrn
  h_epsilon => epsilon
  h_info => info
  h_ipiv => ipiv

  stage = 1

  ! Init device variables
  stat =        cudaMalloc( d_chi0, bytes )
  stat = stat + cudaMalloc( d_krnl, bytes )
  stat = stat + cudaMalloc( d_epsilon, bytes )
  IF( stat /= 0) THEN
     WRITE(*,'("[genvscrn] Failed to allocate device pointers")')
     CALL genvscrn_cleanup( stage ) ! stage 1
     CALL pstop
  END IF
  
  stage = 2

  ! Transfer input argument arrays to device
  stat = cudaMemcpy( d_chi0, C_LOC(h_chi0(1)), bytes, cudaMemcpyHostToDevice )
  stat = cudaMemcpy( d_krnl, C_LOC(h_krnl(1)), bytes, cudaMemcpyHostToDevice )

  ! Init epsilon on device
  ! TODO: There should be a better way, maybe using cudaMemset2D?
  stat = cudaMemset( d_epsilon, 0, bytes )
  address = C_LOC( d_epsilon )
  DO ig = 1, n
     stat = cudaMemset( address, 1, sizeof_complex )
     offset = (n+1)*sizeof_complex
     stat = addOffsetToPtr( address, offset )
  END DO ! ig

  CALL cudaStreamSynchronize( stream )

  ! zgemm: epsilon(ngq,ngq) = 1 * epsilon - 1 * chi0(ngq,ngq) x krnl(ngq,ngq)

  stat = cublasZgemm( handleblas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, -zone, &
                      d_chi0, n, d_krnl, n, zone, d_epsilon, n )

  IF( stat /= 0 ) THEN
     WRITE(*,'("[genvscrn] zgemm I - chi_0 x V failed")')
     CALL genvscrn_cleanup( stage ) ! stage 2
     CALL pstop
  END IF

  CALL cudaStreamSynchronize( stream )

  ! invzge: epsilon -> epsilon^-1
  ! zgemm: vscrn(ngq,ngq) = krnl(ngq,ngq) x epsilon^-1(ngq,ngq)

  ! Instead of actually inverting epsilon, let's just call zgetrf and ztrsm
  ! (zgetri doesn't exist yet in cuSOLVER, and zgetrs solves Ax=B, not xA=B)
  ! zgetrf: epsilon -> LU
  ! ztrsm: solve yU = B, where B = krnl and y is a temporary matrix
  ! ztrsm: solve xL = y, where x = vscrn

  ! First, query the workspace size for zgetrf
  stat = cusolverDnZgetrf_bufferSize( handlesolv, n, n, d_epsilon, n, lwork )
  WRITE(*,'("[genvscrn] Work matrix for zgetrf is ",F7.2," MB")') &
             DBLE(lwork)*sizeof_complex*toMB
  IF (stat /= 0) THEN
     WRITE(*,'("[genvscrn] cusolverDnZgetrf_bufferSize failed, stat = ",I2.2)')&
           stat
     CALL genvscrn_cleanup( stage ) ! stage 2
     CALL pstop
  END IF

  ! Allocate workspaces
  stat =        cudaMalloc( d_work, lwork*sizeof_complex )
  stat = stat + cudaMalloc( d_ipiv, n*sizeof_int )
  stat = stat + cudaMalloc( d_info, sizeof_int )

  IF (stat /= 0) THEN
     WRITE(*,'("[genvscrn] Failed to allocate device workspaces, stat = ",&
             &I2)') stat
     CALL genvscrn_cleanup( stage ) ! stage 2
     CALL pstop
  END IF

  stage = 3
  
  ! Perform LU decomposition
  ! zgetrf: epsilon -> LU
  stat = cusolverDnZgetrf( handlesolv, n, n, d_epsilon, n, d_work, d_ipiv, &
                           d_info)
  ! Send info to host
  stat = cudaMemcpy( h_info, d_info, sizeof_int, cudaMemcpyDeviceToHost )
  IF( info /= 0 ) THEN
     WRITE(*,'("[genvscrn] zgetrf failed, stat = ",I2.2," info = ",I2.2'), &
              stat, info
     CALL genvscrn_cleanup( stage ) ! stage 3
     CALL pstop
  END IF

  ! Send ipiv to host for pivoting
  ! TODO: Find a way to prevent sending ipiv to host
  stat = cudaMemcpy( h_ipiv, d_ipiv, n*sizeof_int, cudaMemcpyDeviceToHost )
  
  ! Perform necessary pivoting operations on B = krnl
  ! TODO: This is basically zlaswp, move it to its own subroutine?
  ! CALL zlaswp( nrhs, b, ldb, 1, n, ipiv, 1 )
  DO ig = 1, n

     IF( ipiv(ig) /= 0 ) THEN

        ! Obtain address for current row
        address = C_LOC( d_krnl )
        offset = (ig-1)*sizeof_complex
        stat = addOffsetToPtr( address, offset )

        ! Obtain address for destination row
        address2 = address
        offset = ( ipiv(ig) )*sizeof_complex
        stat = addOffsetToPtr( address2, offset )

        DO ig1 = 1, n
           ! Perform swap
           stat = cudaMemcpy( C_LOC(d_work), address, sizeof_complex, &
                              cudaMemcpyDevicetoDevice )
           stat = cudaMemcpy( address, address2, sizeof_complex, &
                              cudaMemcpyDevicetoDevice )
           stat = cudaMemcpy( address2, C_LOC(d_work), sizeof_complex, &
                              cudaMemcpyDevicetoDevice )
           ! Increment pointers to next element in row
           offset = n*sizeof_complex
           stat = addOffsetToPtr( address, offset )
           stat = addOffsetToPtr( address2, offset )
        END DO ! swap

     END IF ! ipiv /= 0

  END DO ! zlaswp

  CALL cudaStreamSynchronize( stream )

  ! Solve for temporary matrix y
  ! ztrsm: solve yU = B, where B = krnl and y is a temporary matrix
  stat = cublasZtrsm( handleblas, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, &
                      CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, zzero, &
                      d_epsilon, n, d_krnl, n )
  IF( stat /= 0 ) THEN
     WRITE(*,'("[genvscrn] ztrsm yU=B failed with stat = ",I2.2)') stat
     CALL genvscrn_cleanup( stage ) ! stage 3
     CALL pstop
  END IF
  ! If the previous step succeeded, d_krnl should now contain y

  ! Finally, solve for vscrn
  ! ztrsm: solve xL = y, where x = vscrn and y = temporary matrix
  stat = cublasZtrsm( handleblas, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, &
                      CUBLAS_OP_N, CUBLAS_DIAG_UNIT, n, n, zzero, &
                      d_epsilon, n, d_krnl, n )
  IF( stat /= 0 ) THEN
     WRITE(*,'("[genvscrn] ztrsm xL=y failed with stat = ",I2.2)') stat
     CALL genvscrn_cleanup( stage ) ! stage 3
     CALL pstop
  END IF
  ! If the previous step succeeded, d_krnl should now contain x = vscrn

  CALL cudaStreamSynchronize( stream )

  ! Init epsilon and vscrn on host
  epsilon(:,:) = zzero
  vscrn(:,:)   = zzero

  ! Transfer epsilon and vscrn to host
  stat = cudaMemcpy( C_LOC(h_epsilon), d_epsilon, bytes, cudaMemcpyDeviceToHost)
  stat = cudaMemcpy( C_LOC(h_vscrn), d_krnl, bytes, cudaMemcpyDeviceToHost)

  CALL cudaStreamSynchronize( stream )

  ! Clean up
  CALL genvscrn_cleanup( stage )

  RETURN

CONTAINS
  
  SUBROUTINE genvscrn_cleanup( stage )

    IMPLICIT NONE

    IF( stage > 0 ) THEN
       DEALLOCATE( h_chi0 )
       DEALLOCATE( h_krnl )
       DEALLOCATE( h_vscrn )
       DEALLOCATE( h_epsilon )
       DEALLOCATE( ipiv )
    END IF

    IF( stage > 1 ) THEN
       CALL cudaFree( d_chi0 )
       CALL cudaFree( d_krnl )
       CALL cudaFree( d_epsilon )
    END IF

    IF( stage > 2 ) THEN
       CALL cudaFree( d_work )
       CALL cudaFree( d_ipiv )
       CALL cudaFree( d_info )
    END IF

    RETURN
  END SUBROUTINE genvscrn_cleanup

END SUBROUTINE genvscrn
