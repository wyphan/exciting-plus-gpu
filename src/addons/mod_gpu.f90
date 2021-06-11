MODULE mod_gpu

  USE mod_prec
  USE ISO_C_BINDING

  IMPLICIT NONE

#ifdef NGPUS
  LOGICAL, PARAMETER :: usegpu = .TRUE.
  INTEGER(C_INT) :: ngpus  ! Number of GPUs detected at runtime
  INTEGER(C_INT) :: devnum ! Selected GPU device number
#else
  ! GPU offloading disabled
  LOGICAL, PARAMETER :: usegpu = .FALSE.
  INTEGER(C_INT) :: ngpus = 0
  INTEGER(C_INT) :: devnum = -1
#endif /* NGPUS */

!------------------------------------------------------------------------------
! GPU Programming API
!------------------------------------------------------------------------------

  ! OpenACC
#ifdef _OPENACC
  LOGICAL, PARAMETER :: useacc  = .TRUE.
#else
  LOGICAL, PARAMETER :: useacc  = .FALSE.
#endif /* _OPENACC */

  ! CUDA direct access
#ifdef _CUDA_
  LOGICAL, PARAMETER :: usecuda = .TRUE.
  TYPE(C_PTR) :: stream
#else
  LOGICAL, PARAMETER :: usecuda = .FALSE.
#endif /* _CUDA_ */

  ! HIP
#ifdef _HIP_
  LOGICAL, PARAMETER :: usehip  = .TRUE.
#else
  LOGICAL, PARAMETER :: usehip  = .FALSE.
#endif /* _HIP_ */

!------------------------------------------------------------------------------
! GPU Libraries
!------------------------------------------------------------------------------

  ! MAGMA
#ifdef _MAGMA_
  LOGICAL, PARAMETER :: usemagma = .TRUE.
#else
  LOGICAL, PARAMETER :: usemagma = .FALSE.
#endif /* _MAGMA_ */
  TYPE(C_PTR) :: queue ! MAGMA queue
  
  ! cuBLAS
#ifdef _CUBLAS_
  LOGICAL, PARAMETER :: usecublas = .TRUE.
#else
  LOGICAL, PARAMETER :: usecublas = .FALSE.
#endif /* _CUBLAS_ */
  TYPE(C_PTR) :: blashandle
  
  ! cuSPARSE
#ifdef _CUSPARSE_
  LOGICAL, PARAMETER :: usecusparse = .TRUE.
#else
  LOGICAL, PARAMETER :: usecusparse = .FALSE.
#endif /* _CUBLAS_ */
  TYPE(C_PTR) :: sparsehandle

  ! cuSolver
#ifdef _CUSOLVERDN_
  LOGICAL, PARAMETER :: usecusolverdn = .TRUE.
#else
  LOGICAL, PARAMETER :: usecusolverdn = .FALSE.
#endif /* _CUSOLVERDN_ */
  TYPE(C_PTR) :: densehandle

!------------------------------------------------------------------------------
! GPU array sizes in bytes
!------------------------------------------------------------------------------

  ! mod_addons_q
  INTEGER(KIND=dl) :: sz_sfacgq = 0

  ! mod_nrkp
  INTEGER(KIND=dl) :: sz_wfsvmtnrloc = 0

  ! mod_expigqr
#ifdef _PACK_gntuju_
  INTEGER(KIND=dl) :: sz_gntuju_packed = 0
#else
  INTEGER(KIND=dl) :: sz_gntuju = 0
#endif /* _PACK_gntuju_ */

  ! mod_expigqr::genmegq()
  INTEGER(KIND=dl) :: sz_wfsvmt_jk = 0

  ! mod_genmegqblh_gpu
  INTEGER(KIND=dl) :: sz_b1 = 0
  INTEGER(KIND=dl) :: sz_b2 = 0
  INTEGER(KIND=dl) :: sz_wftmp1mt = 0

!------------------------------------------------------------------------------
! GPU FLOP estimate
!------------------------------------------------------------------------------

  INTEGER(KIND=dl) :: flop_fillbatch = 0
  INTEGER(KIND=dl) :: flop_batchzgemm = 0

!==============================================================================

CONTAINS

!==============================================================================
! Counts the number of GPUs (ngpus) and selects one (devnum)

  SUBROUTINE gpu_enumerate
    USE ISO_C_BINDING

#if defined(_OPENACC)
    USE openacc
#endif /* _OPENACC */

#if defined(_CUDA_)
    !USE mod_cuda
#endif /* _CUDA_ */

#if defined(_ROCM_)
    !USE mod_rocm
#endif /* _ROCM_ */

#if defined(_HIP_)
    !USE mod_hip
#endif /* _HIP_ */

    IMPLICIT NONE

    ! Internal variables
    INTEGER :: devtype ! OpenACC device type
    INTEGER(C_INT) :: stat ! cudaError_t

    IF( usegpu ) THEN

#if defined(_OPENACC)

       ! Use OpenACC runtime procedures
       ! For now, assume NVIDIA GPUs
       ! TODO: generalize for AMD GPUs
       devtype = acc_device_nvidia
       ngpus = acc_get_num_devices( devtype )
       devnum = acc_get_device_num( devtype )

#elif defined(_CUDA_)

       ! stat = cudaGetDeviceCount( ngpus )
       ! stat = cudaGetDevice( devnum )
       ! stat = cudaStreamCreate( stream )

#elif defined(_ROCM_)

       ! Not implemented

#endif /* GPU API */

    END IF

  END SUBROUTINE gpu_enumerate

!==============================================================================
! Initialize GPU libraries

  SUBROUTINE gpu_init_libs

#ifdef _MAGMA_
    USE mod_magma
#endif

#ifdef _CUDA_
    !USE cuda
#endif

    IMPLICIT NONE

    ! Internal variables
    INTEGER(C_INT) :: stat ! cudaError_t

    !IF( usecuda ) THEN
    !   stat = cudaStreamCreate( stream )
    !   IF( usecublas ) THEN
    !      stat = cublasCreate( blashandle )
    !   END IF ! usecublas
    !   IF( usesparse ) THEN
    !      stat = cusparseCreate( sparsehandle )
    !   END IF ! usecusparse
    !   IF( usecusolverdn ) THEN
    !      stat = cusolverDnCreate( densehandle )
    !   END IF ! usecusolver
    !END IF ! usecuda

    IF( usemagma ) THEN
    !   IF( usecublas .AND. usecusparse ) THEN
    !      CALL magma_init()
    !      CALL magma_setdevice( devnum )
    !      CALL magma_queue_create_from_cuda( devnum, stream, &
    !                                         blashandle, sparsehandle, queue )
    !   ELSE
#ifdef _MAGMA_   
          ! Only the master thread has access to MAGMA
          ! TODO: test thread safety
          !$OMP MASTER
          CALL magma_init()
          CALL magma_set_device( devnum )
          CALL magma_queue_create( devnum, queue )
          !$OMP END MASTER
#endif /* _MAGMA_ */
    !   END IF
    END IF ! usemagma

    RETURN
  END SUBROUTINE gpu_init_libs

!==============================================================================
! Finalize GPU libraries

  SUBROUTINE gpu_fin_libs

#ifdef _MAGMA_
    USE mod_magma
#endif

#ifdef _CUDA_
    !USE cuda
#endif

    ! Internal variables
    INTEGER(C_INT) :: stat ! cudaError_t

    IF( usemagma ) THEN
#ifdef _MAGMA_
       ! Only the master thread has access to MAGMA
       ! TODO: test thread safety
       !$OMP MASTER
       CALL magma_queue_destroy( queue )
       CALL magma_finalize()
       !$OMP END MASTER
#endif /* _MAGMA_ */       
    END IF

    !IF( usecuda ) THEN
    !   IF( usecublas ) THEN
    !      stat = cublasDestroy( blashandle )
    !   END IF ! usecublas
    !   IF( usecusparse ) THEN
    !      stat = cusparseDestroy( sparsehandle )
    !   END IF ! usecublas
    !   IF( usecusolverdn ) THEN
    !      stat = cusolverDnDestroy( densehandle )
    !   END IF ! usecublas
    !   stat = cudaStreamDestroy( stream )
    !END IF ! usecuda

    RETURN
  END SUBROUTINE gpu_fin_libs

!==============================================================================
! ZGEMM using OpenACC or OpenMP target
  
      subroutine ZGEMM_acc(transA,transB, m,n,kk,                            &
     &   alpha,A,lda,B,ldb,beta,C,ldc)
      use mod_prec
      implicit none
#ifdef _OPENACC
!$acc routine vector
#else
!$omp declare target  
#endif
      integer, intent(in) :: m,n,kk,lda,ldb,ldc
      COMPLEX(KIND=dz), intent(in) :: alpha,beta
      character, intent(in) :: transA, transB

      COMPLEX(KIND=dz), intent(in) :: A(lda,*)
      COMPLEX(KIND=dz), intent(in) :: B(ldb,*)
      COMPLEX(KIND=dz), intent(inout) :: C(ldc,*)

      integer :: i,j,k
      COMPLEX(KIND=dz) :: cij, aik, bkj
      logical :: is_Aconj,is_Atrans,is_AN,no_transA
      logical :: is_Bconj,is_Btrans,is_BN,no_transB

      is_Aconj = (transA.eq.'C').or.(transA.eq.'c')
      is_Bconj = (transB.eq.'C').or.(transB.eq.'c')
      is_Atrans = (transA.eq.'T').or.(transA.eq.'t')
      is_Btrans = (transB.eq.'T').or.(transB.eq.'t')

      is_AN = (transA.eq.'N').or.(transA.eq.'n')
      is_BN = (transB.eq.'N').or.(transB.eq.'n')

      no_transA = is_AN .or. (.not.(is_Aconj.or.is_Atrans))
      no_transB = is_BN .or. (.not.(is_Bconj.or.is_Btrans))

#ifdef _OPENACC
!$acc loop vector collapse(2) private(cij,aik,bkj,k)
#else
!$omp parallel do simd collapse(2) private(cij,aik,bkj,k)
#endif
      do j=1,n
      do i=1,m

         cij = 0
         do k=1,kk
            if (no_transA) then
               aik = A(i,k)
            else
               if (is_Aconj) then
                  aik = conjg(A(k,i))
               else
                  aik = A(k,i)
               endif
            endif
            
            if (no_transB) then
               bkj = B(k,j)
            else
               if (is_Bconj) then
                  bkj = conjg(B(j,k))
               else
                  bkj = B(j,k)
               endif
            endif
            
            cij = cij + aik * bkj
         enddo
         
         if (beta.eq.0) then
            C(i,j) = alpha * cij
         else
            C(i,j) = beta*C(i,j) + alpha * cij
         endif

      enddo
      enddo
#ifdef _OPENACC
!$acc end loop
#else
!$omp end parallel do simd
#endif

      return
      end subroutine ZGEMM_acc

!==============================================================================
! Batched ZGEMM using OpenACC and MAGMA
! The arrays should already be in device memory;
! this subroutine will internally expose the device pointers

  SUBROUTINE zgemm_batched_gpu_acc_magma( transA, transB, m, n, k, &
                                          alpha, dA_r, ldda, &
                                                 dB_r, lddb, &
                                          beta,  dC_r, lddc, &
                                          batchCount )
#ifdef _MAGMA_
    ! Batched zgemm is not available in "magma" module
    USE mod_magma
#endif /* _MAGMA_ */

    USE ISO_C_BINDING
    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: transA, transB
    INTEGER, INTENT(IN) :: m, n, k, ldda, lddb, lddc, batchCount
    COMPLEX(KIND=dz), INTENT(IN) :: alpha, beta
    COMPLEX(KIND=dz), DIMENSION(:,:,:), INTENT(IN), TARGET :: dA_r, dB_r
    COMPLEX(KIND=dz), DIMENSION(:,:,:), INTENT(INOUT), TARGET :: dC_r

#ifdef _MAGMA_

    ! Internal variables
    INTEGER :: ierr, ibatch
    INTEGER(KIND=C_INT) :: op_a, op_b
    INTEGER(KIND=C_INT) :: h_m, h_n, h_k, h_ldda, h_lddb, h_lddc, d_batchCount
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: dptr_A, dptr_B, dptr_C

    ! TODO: test thread safety
    !$OMP MASTER
    
    ! Map transA and transB to enum
    op_a = magma_trans_const( transA )
    op_b = magma_trans_const( transB )

    ALLOCATE( dptr_A( batchCount ))
    ALLOCATE( dptr_B( batchCount ))
    ALLOCATE( dptr_C( batchCount ))

    ! Check arguments and initialize device pointer arrays
    !$ACC DATA PRESENT_OR_COPYIN( dA_r, dB_r, dC_r ) &
    !$ACC      CREATE( dptr_A, dptr_B, dptr_C )

    ! Convert integer arguments
    h_m = m
    h_n = n
    h_k = k
    h_ldda = ldda
    h_lddb = lddb
    h_lddc = lddc
    d_batchCount = batchCount

    ! Copy to device
    !$ACC DATA COPYIN( d_batchCount )
    
    ! Expose device pointers
    !$ACC HOST_DATA USE_DEVICE( dA_r, dB_r, dC_r )

    ! Extract device pointers
    !$ACC KERNELS LOOP PRIVATE(ibatch)
    DO ibatch = 1, d_batchCount
       dptr_A(ibatch) = C_LOC( dA_r( LBOUND(dA_r,1), LBOUND(dA_r,2), ibatch) )
       dptr_B(ibatch) = C_LOC( dB_r( LBOUND(dB_r,1), LBOUND(dB_r,2), ibatch) )
       dptr_C(ibatch) = C_LOC( dC_r( LBOUND(dC_r,1), LBOUND(dC_r,2), ibatch) )
    END DO
    !$ACC END KERNELS

    ! Expose both the device pointers and the array of pointers on device
    !$ACC END HOST_DATA
    !$ACC HOST_DATA USE_DEVICE( dA_r, dB_r, dC_r, dptr_A, dptr_B, dptr_C )

    ! Call MAGMA with extracted device pointers
    CALL magmablas_zgemm_batched( op_a, op_b, h_m, h_n, h_k, &
                                  alpha, C_LOC(dptr_A), h_ldda, &
                                         C_LOC(dptr_B), h_lddb, &
                                  beta,  C_LOC(dptr_C), h_lddc, &
                                  d_batchCount, queue )

    ! dA_r, dB_r, dC_r, dptr_A, dptr_B, dptr_C
    !$ACC END HOST_DATA

    ! d_batchCount
    !$ACC END DATA

    ! everything else
    !$ACC END DATA

    DEALLOCATE( dptr_A )
    DEALLOCATE( dptr_B )
    DEALLOCATE( dptr_C )

    !$OMP END MASTER

#endif /* _MAGMA_ */

    RETURN
  END SUBROUTINE zgemm_batched_gpu_acc_magma

!==============================================================================
! Batched ZGEMM using OpenACC and MAGMA
! The arrays and their pointer arrays should already be in device memory;
! this subroutine is simply calls MAGMA, passing along the pointer arrays
  SUBROUTINE zgemm_batched_gpu_acc_magma_ptr( transA, transB, m, n, k, &
                                              alpha, dptr_A, ldda, &
                                                     dptr_B, lddb, &
                                              beta,  dptr_C, lddc, &
                                              batchCount )
#ifdef _MAGMA_
    ! Batched zgemm is not available in "magma" module
    USE mod_magma
#endif /* _MAGMA_ */

    USE mod_prof

    USE ISO_C_BINDING
    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: transA, transB
    INTEGER, INTENT(IN) :: m, n, k, ldda, lddb, lddc, batchCount
    COMPLEX(KIND=dz), INTENT(IN) :: alpha, beta
    TYPE(C_PTR), DIMENSION(batchCount), INTENT(IN) :: dptr_A, dptr_B
    TYPE(C_PTR), DIMENSION(batchCount), INTENT(INOUT) :: dptr_C
#ifdef _CUDA
    ATTRIBUTES(DEVICE) :: dptr_A, dptr_B, dptr_C
#endif /* _CUDA */

#ifdef _MAGMA_

    ! Internal variables
    INTEGER :: ierr, ibatch
    INTEGER(KIND=C_INT) :: op_a, op_b
    INTEGER(KIND=C_INT) :: h_m, h_n, h_k, h_ldda, h_lddb, h_lddc, d_batchCount

    ! TODO: test thread safety
    !$OMP MASTER
    
    ! Map transA and transB to enum
    op_a = magma_trans_const( transA )
    op_b = magma_trans_const( transB )

    ! Convert integer arguments
    h_m = m
    h_n = n
    h_k = k
    h_ldda = ldda
    h_lddb = lddb
    h_lddc = lddc
    d_batchCount = batchCount

    ! Denote arguments that are already on device
    !$ACC DATA PRESENT( dptr_A, dptr_B, dptr_C )

    CALL profstart( "magmablas_zgemm_batched" )
    
    ! Call MAGMA with device pointer arrays
    CALL magmablas_zgemm_batched( op_a, op_b, &
                                  h_m, h_n, h_k, &
                                  alpha, C_LOC(dptr_A), h_ldda, &
                                         C_LOC(dptr_B), h_lddb, &
                                  beta,  C_LOC(dptr_C), h_lddc, &
                                  d_batchCount, queue )

    ! dptr_A, dptr_B, dptr_C
    !$ACC END DATA

    CALL magma_queue_sync( queue )

    CALL profend( "magmablas_zgemm_batched" )

    !$OMP END MASTER

#endif /* _MAGMA_ */

    RETURN
  END SUBROUTINE zgemm_batched_gpu_acc_magma_ptr

!==============================================================================
! Batched ZGEMM using OpenACC and MAGMA (with variable dimensions)
! The arrays and their pointer arrays should already be in device memory;
! this subroutine is simply calls MAGMA, passing along the pointer arrays
  SUBROUTINE zgemm_vbatched_gpu_acc_magma_ptr( transA, transB, m, n, k, &
                                              alpha, dptr_A, ldda, &
                                                     dptr_B, lddb, &
                                              beta,  dptr_C, lddc, &
                                              batchCount )
#ifdef _MAGMA_
    ! Batched zgemm is not available in "magma" module
    USE mod_magma
#endif /* _MAGMA_ */

    USE ISO_C_BINDING
    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: transA, transB
    INTEGER, DIMENSION(batchCount+1), INTENT(IN) :: m, n, k, ldda, lddb, lddc
    COMPLEX(KIND=dz), INTENT(IN), VALUE :: alpha, beta
    TYPE(C_PTR), DIMENSION(batchCount), INTENT(IN) :: dptr_A, dptr_B
    TYPE(C_PTR), DIMENSION(batchCount), INTENT(INOUT) :: dptr_C
    INTEGER, INTENT(IN), VALUE :: batchCount

#ifdef _MAGMA_

    ! Internal variables
    INTEGER :: ierr, ibatch
    INTEGER(KIND=C_INT) :: op_a, op_b
    INTEGER(KIND=C_INT) :: h_batchCount
    INTEGER(KIND=C_INT), DIMENSION(batchCount+1) :: d_m, d_n, d_k
    INTEGER(KIND=C_INT), DIMENSION(batchCount+1) :: d_ldda, d_lddb, d_lddc
    
    ! TODO: test thread safety
    !$OMP MASTER
    
    ! Map transA and transB to enum
    op_a = magma_trans_const( transA )
    op_b = magma_trans_const( transB )

    ! Convert integer arguments
    d_m(:) = m(:)
    d_n(:) = n(:)
    d_k(:) = k(:)
    d_ldda(:) = ldda(:)
    d_lddb(:) = lddb(:)
    d_lddc(:) = lddc(:)
    h_batchCount = batchCount
    
    ! Check arguments
    !$ACC DATA DEVICEPTR( dptr_A, dptr_B, dptr_C ) &
    !$ACC      COPYIN( d_m, d_n, d_k, d_ldda, d_lddb, d_lddc )

    ! Expose device pointers
    !$ACC HOST_DATA USE_DEVICE( dptr_A, dptr_B, dptr_C, &
    !$ACC                       d_m, d_n, d_k, d_ldda, d_lddb, d_lddc )

    ! Call MAGMA with device pointer arrays
    CALL magmablas_zgemm_vbatched( op_a, op_b, &
                                   d_m, d_n, d_k, &
                                   alpha, dptr_A, d_ldda, &
                                          dptr_B, d_lddb, &
                                   beta,  dptr_C, d_lddc, &
                                   h_batchCount, queue )

    ! dptr_A, dptr_B, dptr_C, d_m, d_n, d_k, d_ldda, d_lddb, d_lddc
    !$ACC END HOST_DATA
    !$ACC END DATA

    CALL magma_queue_sync( queue )

    !$OMP END MASTER

#endif /* _MAGMA_ */

    RETURN
  END SUBROUTINE zgemm_vbatched_gpu_acc_magma_ptr

!==============================================================================
! Fallback mechanism: Batched ZGEMM on CPU using OpenMP parallel do
! Each thread operates on a different batch

  SUBROUTINE zgemm_batched_omp( transA, transB, m, n, k, &
                                alpha, A_array, lda, &
                                       B_array, ldb, &
                                beta,  C_array, ldc, &
                                batchCount )
    USE mod_lapack, only: ZGEMM
#ifdef _OPENMP
    USE omp_lib
#endif /* _OPENMP */

    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: transA, transB
    INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc, batchCount
    COMPLEX(KIND=dz), INTENT(IN) :: alpha, beta
    COMPLEX(KIND=dz), DIMENSION(:,:,:), INTENT(IN) :: A_array, B_array
    COMPLEX(KIND=dz), DIMENSION(:,:,:), INTENT(INOUT) :: C_array

    ! Internal variables
    !COMPLEX(KIND=dz), DIMENSION(SIZE(A_array,1),SIZE(A_array,2)) :: matA
    !COMPLEX(KIND=dz), DIMENSION(SIZE(B_array,1),SIZE(B_array,2)) :: matB
    !COMPLEX(KIND=dz), DIMENSION(SIZE(C_array,1),SIZE(C_array,2)) :: matC
    INTEGER :: ld1, ld2, ld3
    INTEGER :: ibatch, tid

    !WRITE(*,*) 'zgemm_batched_omp: batchCount=', batchCount

    ld1 = size(A_array,1)
    ld2 = size(B_array,1)
    ld3 = size(C_array,1)
    
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP   PRIVATE( ibatch, tid )
!    !$OMP   PRIVATE( ibatch, tid, matA, matB, matC )
    DO ibatch = 1, batchCount

!--DEBUG
#if EBUG > 2 && defined(_OPENMP)
       tid = OMP_GET_THREAD_NUM()
       WRITE(*,*) 'zgemm_batched_omp: Thread ', tid, ' executes batch ', ibatch
#endif /* DEBUG */
!--DEBUG

       ! Fetch arrays
       !matA(:,:) = A_array(:,:,ibatch)
       !matB(:,:) = B_array(:,:,ibatch)
       !matC(:,:) = C_array(:,:,ibatch)

       ! Call ZGEMM (let BLAS check the arguments)
       CALL zgemm( transA, transB, m, n, k, &
                   alpha, A_array(:,:,ibatch), ld1, &
                          B_array(:,:,ibatch), ld2, &
                   beta,  C_array(:,:,ibatch), ld3 )

       ! Save result
!       !$OMP CRITICAL
       !C_array(:,:,ibatch) = matC(:,:)
!       !$OMP END CRITICAL

    END DO ! ibatch

    RETURN
  END SUBROUTINE zgemm_batched_omp

!==============================================================================
! Fallback mechanism: Batched & strided ZGEMM on CPU using OpenMP parallel do
! Each thread operates on a different batch

  SUBROUTINE zgemm_strided_batched_omp( transA, transB, m, n, k, &
                                        alpha, A_array, lda, strideA, &
                                               B_array, ldb, strideB, &
                                        beta,  C_array, ldc, strideC, &
                                        batchCount )                                
    USE mod_lapack, only: ZGEMM
#ifdef _OPENMP
    USE omp_lib
#endif /* _OPENMP */

    IMPLICIT NONE

    ! Arguments
    CHARACTER(LEN=1), INTENT(IN) :: transA, transB
    INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc, batchCount
    INTEGER, INTENT(IN) :: strideA, strideB, strideC
    COMPLEX(KIND=dz), INTENT(IN) :: alpha, beta
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(IN) :: A_array, B_array
    COMPLEX(KIND=dz), DIMENSION(*), INTENT(INOUT) :: C_array

    ! Internal variables
    INTEGER :: ld1, ld2, ld3
    INTEGER :: ipA, ipB, ipC, ibatch, tid
    INTEGER :: nrowsA, nrowsB, nrowsC
    LOGICAL :: ltransA, ltransB, lokA, lokB, lokC

    !WRITE(*,*) 'zgemm_batched_strided_omp: batchCount=', batchCount
    !WRITE(*,*) 'zgemm_batched_strided_omp: strideA=', strideA
    !WRITE(*,*) 'zgemm_batched_strided_omp: strideB=', strideB
    !WRITE(*,*) 'zgemm_batched_strided_omp: strideC=', strideC

    ! Set nrows
    ltransA = ( transA == 'T' ) .OR. ( transA == 't' )
    ltransB = ( transB == 'T' ) .OR. ( transB == 't' )
    nrowsA = MERGE( k, m, ltransA )
    nrowsB = MERGE( n, k, ltransB )
    nrowsC = m

    ! Check dimensions
    lokA = ( 1 <= nrowsA ) .AND. ( nrowsA <= lda )
    lokB = ( 1 <= nrowsB ) .AND. ( nrowsB <= ldb )
    lokC = ( 1 <= nrowsC ) .AND. ( nrowsC <= ldc )
    IF( .NOT. lokA ) THEN
       WRITE(*,*) 'Error[zgemm_strided_batched_omp]: # of rows of A ', nrowsA, &
                  'is out of bounds (ldA = ', lda, ' )'
       STOP
    END IF
    IF( .NOT. lokB ) THEN
       WRITE(*,*) 'Error[zgemm_strided_batched_omp]: # of rows of B ', nrowsB, &
                  'is out of bounds (ldB = ', ldb, ' )'
       STOP
    END IF
    IF( .NOT. lokC ) THEN
       WRITE(*,*) 'Error[zgemm_strided_batched_omp]: # of rows of C ', nrowsC, &
                  'is out of bounds (ldC = ', ldc, ' )'
       STOP
    END IF

    ld1 = lda
    ld2 = ldb
    ld3 = ldc

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP   PRIVATE( ibatch, ipA, ipB, ipC, tid )
    DO ibatch = 1, batchCount

!--DEBUG
#if EBUG > 2 && defined(_OPENMP)
       tid = OMP_GET_THREAD_NUM()
       WRITE(*,*) 'zgemm_batched_strided_omp: Thread ', tid, &
                  ' executes batch ', ibatch
#endif /* DEBUG */
!--DEBUG

       ! Set up pointer location
       ipA = 1 + ( ibatch - 1 ) * strideA
       ipB = 1 + ( ibatch - 1 ) * strideB
       ipC = 1 + ( ibatch - 1 ) * strideC

       ! Call ZGEMM (let BLAS check the arguments)
       CALL ZGEMM( transA, transB, m, n, k, &
                   alpha, A_array(ipA), ld1, &
                          B_array(ipB), ld2, &
                   beta,  C_array(ipC), ld3 )

    END DO ! ibatch

    RETURN
  END SUBROUTINE zgemm_strided_batched_omp

!==============================================================================

END MODULE mod_gpu
