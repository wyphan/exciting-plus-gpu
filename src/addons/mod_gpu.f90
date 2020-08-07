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
  INTEGER(C_INT), PARAMETER :: ngpus  = 0
  INTEGER(C_INT), PARAMETER :: devnum = -1
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
    INTEGER(KIND=C_INT) :: h_m, h_n, h_k, h_ldda, h_lddb, h_lddc, h_batchCount
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: dptr_a, dptr_b, dptr_c

    ! TODO: test thread safety
    !$OMP MASTER
    
    ! Map transA and transB to enum
    op_a = magma_trans_const( transA )
    op_b = magma_trans_const( transB )

    ALLOCATE( dptr_a( batchCount ))
    ALLOCATE( dptr_b( batchCount ))
    ALLOCATE( dptr_c( batchCount ))

    ! Check arguments
    !$ACC DATA PRESENT( dA_r, dB_r, dC_r ) CREATE( dptr_a, dptr_b, dptr_c )

    ! Convert integer arguments
    h_m = m
    h_n = n
    h_k = k
    h_ldda = ldda
    h_lddb = lddb
    h_lddc = lddc
    h_batchCount = batchCount

    ! Expose device pointers
    !$ACC HOST_DATA USE_DEVICE( dA_r, dB_r, dC_r )

    ! Extract device pointers
    !$ACC KERNELS LOOP PRIVATE(ibatch)
    DO ibatch = 1, batchCount
       dptr_a(ibatch) = C_LOC( dA_r( LBOUND(dA_r,1), LBOUND(dA_r,2), ibatch) )
       dptr_b(ibatch) = C_LOC( dB_r( LBOUND(dB_r,1), LBOUND(dB_r,2), ibatch) )
       dptr_c(ibatch) = C_LOC( dC_r( LBOUND(dC_r,1), LBOUND(dC_r,2), ibatch) )
    END DO
    !$ACC END KERNELS

    ! Expose both the device pointers and the array of pointers on device
    !$ACC END HOST_DATA
    !$ACC HOST_DATA USE_DEVICE( dA_r, dB_r, dC_r, dptr_a, dptr_b, dptr_c )

    ! Call MAGMA with extracted device pointers
    CALL magmablas_zgemm_batched( op_a, op_b, h_m, h_n, h_k, &
                                  alpha, C_LOC(dptr_a), h_ldda, &
                                         C_LOC(dptr_b), h_lddb, &
                                  beta,  C_LOC(dptr_c), h_lddc, &
                                  h_batchCount, queue )

    ! dA_r, dB_r, dC_r, dptr_a, dptr_b, dptr_c
    !$ACC END HOST_DATA

    ! everything else
    !$ACC END DATA
    DEALLOCATE( dptr_a )
    DEALLOCATE( dptr_b )
    DEALLOCATE( dptr_c )

    !$OMP END MASTER

#endif /* _MAGMA_ */

    RETURN
  END SUBROUTINE zgemm_batched_gpu_acc_magma

!==============================================================================
! Fallback mechanism: Batched ZGEMM on CPU using OpenMP parallel do
! Each thread operates on a different batch

  SUBROUTINE zgemm_batched_omp( transA, transB, m, n, k, &
                                alpha, A_array, lda, &
                                       B_array, ldb, &
                                beta,  C_array, ldc, &
                                batchCount )
    USE modmain, only: zzero, zone
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
    INTEGER :: lb11, lb12, lb21, lb22, lb31, lb32
    INTEGER :: ub11, ub12, ub21, ub22, ub31, ub32
    INTEGER :: ibatch, tid

    !WRITE(*,*) 'zgemm_batched_omp: batchCount=', batchCount

    ld1 = SIZE(A_array,1)
    lb11 = LBOUND(A_array,1);  ub11 = UBOUND(A_array,1)
    lb12 = LBOUND(A_array,2);  ub12 = UBOUND(A_array,2)

    ld2 = SIZE(B_array,1)
    lb21 = LBOUND(B_array,1);  ub21 = UBOUND(B_array,1)
    lb22 = LBOUND(B_array,2);  ub22 = UBOUND(B_array,2)

    ld3 = SIZE(C_array,1)
    lb31 = LBOUND(C_array,1);  ub31 = UBOUND(C_array,1)
    lb32 = LBOUND(C_array,2);  ub32 = UBOUND(C_array,2)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP   PRIVATE( ibatch, tid )
    DO ibatch = 1, batchCount

!--DEBUG
#if EBUG > 2 && defined(_OPENMP)
       tid = OMP_GET_THREAD_NUM()
       WRITE(*,*) 'zgemm_batched_omp: Thread ', tid, ' executes batch ', ibatch
#endif /* DEBUG */
!--DEBUG

       ! Call ZGEMM (let BLAS check the arguments)
       CALL zgemm( transA, transB, m, n, k, &
            alpha, A_array( lb11:ub11, lb12:ub12, ibatch), ld1, &
                   B_array( lb21:ub21, lb22:ub22, ibatch), ld2, &
            beta,  C_array( lb31:ub31, lb32:ub32, ibatch), ld3 )

    END DO ! ibatch

    RETURN
  END SUBROUTINE zgemm_batched_omp

!==============================================================================

END MODULE mod_gpu
