module cublas_f
  use ISO_C_BINDING

    enum, BIND(C)
        enumerator :: CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C
    end enum
   
    enum, BIND(C)
       enumerator :: cudaMemcpyHostToHost, cudaMemcpyHostToDevice, &
                     cudaMemcpyDeviceToHost, cudaMemcpyDeviceToDevice, &
                     cudaMemcpyDefault
    end enum

!------------------------------------------------------------------------------

    ENUM, BIND(C)
       ENUMERATOR :: CUBLAS_SIDE_LEFT, CUBLAS_SIDE_RIGHT
    END ENUM

    ENUM, BIND(C)
       ENUMERATOR :: CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_UPPER
    END ENUM

    ENUM, BIND(C)
       ENUMERATOR :: CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT
    END ENUM

    ! CUDA stream and handlers
    TYPE(C_PTR), POINTER :: stream, handleblas, handlesolv

    ! Data sizes in bytes
    INTEGER, PARAMETER :: sizeof_ptr = 8
    INTEGER, PARAMETER :: sizeof_int = 4
    INTEGER, PARAMETER :: sizeof_int8 = 8
    INTEGER, PARAMETER :: sizeof_complex = 16
    REAL(8), PARAMETER :: toMB = 1.d0/1024.d0/1024.d0

!------------------------------------------------------------------------------

  INTERFACE
    integer(C_INT) function cudaMalloc(ptr, bytes) BIND(C, NAME='cudaMalloc')
        use ISO_C_BINDING
        type(C_PTR) :: ptr
        integer(C_SIZE_T), value :: bytes
    end function

    integer(C_INT) function cudaMemcpy(dst, src, count, kind) BIND(C, NAME='cudaMemcpy')
        use ISO_C_BINDING
        type(C_PTR), value :: dst
        type(C_PTR), value :: src
        integer(C_SIZE_T), value :: count
        integer(C_INT), value :: kind
    end function

    integer(C_INT) function cublasCreate(handle_ptr) BIND(C, NAME='f_cublasCreate')
        use ISO_C_BINDING
        type(C_PTR) :: handle_ptr
    end function

    subroutine  cublasDestroy(handle_ptr) BIND(C, NAME='f_cublasDestroy')
        use ISO_C_BINDING
        type(C_PTR), value :: handle_ptr
    end subroutine

    subroutine cudaFree(ptr) BIND(C, NAME='cudaFree')
        use ISO_C_BINDING
        type(C_PTR), value :: ptr
    end subroutine

    integer(C_INT) function cublasSetMatrixAsync(rows, cols, elemSize, a_ptr, &
                                           lda, b_ptr, ldb, stream) &
                                           BIND(C, NAME='cublasSetMatrix')
        use ISO_C_BINDING
        integer(C_INT), value :: rows
        integer(C_INT), value :: cols
        integer(C_INT), value :: elemSize
        type(C_PTR),    value :: a_ptr
        integer(C_INT), value :: lda
        type(C_PTR),    value :: b_ptr
        integer(C_INT), value :: ldb
        type(C_PTR), value :: stream
    end function

    integer(C_INT) function cublasGetMatrixAsync(rows, cols, elemSize, a_ptr, &
                                            lda, b_ptr, ldb, stream) &
                                            BIND(C, NAME='cublasGetMatrix')
        use ISO_C_BINDING
        integer(C_INT), value :: rows
        integer(C_INT), value :: cols
        integer(C_INT), value :: elemSize
        type(C_PTR),    value :: a_ptr
        integer(C_INT), value :: lda
        type(C_PTR),    value :: b_ptr
        integer(C_INT), value :: ldb
        type(C_PTR), value :: stream

    end function

    integer(C_INT) function cublasZgemm(handle, transa, transb, m, n, k, alpha, &
                                        A, lda, B, ldb, beta, C, ldc) &
                                        BIND(C, NAME='f_cublasZgemm')
        use ISO_C_BINDING
        type(C_PTR), value    :: handle
        integer(C_INT), value :: transa
        integer(C_INT), value :: transb
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        complex(C_DOUBLE_COMPLEX)        :: alpha
        type(C_PTR), value    :: A
        integer(C_INT), value :: lda
        type(C_PTR), value    :: B
        integer(C_INT), value :: ldb
        complex(C_DOUBLE_COMPLEX)        :: beta
        type(C_PTR), value    :: C
        integer(C_INT), value :: ldc
    end function

    integer(C_INT) function cublasZgemmBatched(handle, transa, transb, m, n, k, alpha, &
                                        A, lda, B, ldb, beta, C, ldc, batch_count) &
                                        BIND(C, NAME='f_cublasZgemmBatched')
        use ISO_C_BINDING
        type(C_PTR), value    :: handle
        integer(C_INT), value :: transa
        integer(C_INT), value :: transb
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        complex(C_DOUBLE_COMPLEX)        :: alpha
        type(C_PTR), value    :: A
        integer(C_INT), value :: lda
        type(C_PTR), value    :: B
        integer(C_INT), value :: ldb
        complex(C_DOUBLE_COMPLEX)        :: beta
        type(C_PTR), value    :: C
        integer(C_INT), value :: ldc
        integer(C_INT), value :: batch_count
    end function

    integer(C_INT) function cudaStreamCreate(stream_ptr) BIND(C, NAME='f_cudaStreamCreate')
        use ISO_C_BINDING
        type(C_PTR) :: stream_ptr
    end function

    integer(C_INT) function cublasSetStream(handle, stream) BIND(C, NAME='f_cublasSetStream')
        use ISO_C_BINDING
        type(C_PTR), value :: handle
        type(C_PTR), value :: stream
    end function

    integer(C_INT) function cudaStreamDestroy(stream) BIND(C, NAME='f_cudaStreamDestroy')
        use ISO_C_BINDING
        type(C_PTR), value :: stream
    end function

    subroutine cudaDeviceSynchronize() BIND(C, NAME='cudaDeviceSynchronize')
    end subroutine

    subroutine cudaStreamSynchronize(stream) BIND(C, NAME='cudaStreamSynchronize')
        use ISO_C_BINDING
        type(C_PTR), value :: stream
    end subroutine

    integer(C_INT) function cudaMemset(ptr, val, bytes) BIND(C, NAME='cudaMemset')
        use ISO_C_BINDING
        type(C_PTR) :: ptr
        integer(C_SIZE_T), value :: val
        integer(C_SIZE_T), value :: bytes
    end function

    integer(C_INT) function addOffsetToPtr(array, offset) BIND(C, NAME='f_addOffsetToPtr')
        use ISO_C_BINDING
        type(C_PTR) :: array
        integer(C_SIZE_T), value :: offset
    end function

    integer(C_INT) function printValue(address) BIND(C, NAME='f_printValue')
        use ISO_C_BINDING
        type(C_PTR) :: address
    end function

!------------------------------------------------------------------------------

    INTEGER(C_INT) FUNCTION cusolverDnCreate(handle_ptr) &
                   BIND(C, NAME='f_cusolverDnCreate')
      USE ISO_C_BINDING
      TYPE(C_PTR) :: handle_ptr
    END FUNCTION cusolverDnCreate

    INTEGER(C_INT) FUNCTION cusolverDnDestroy(handle_ptr) &
                   BIND(C, NAME='f_cusolverDnDestroy')
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle_ptr
    END FUNCTION cusolverDnDestroy

    INTEGER(C_INT) FUNCTION cusolverDnSetStream(handle_ptr, stream) &
                   BIND(C, NAME='f_cusolverDnSetStream')
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle_ptr
      TYPE(C_PTR), VALUE :: stream
    END FUNCTION cusolverDnSetStream

    INTEGER(C_INT) FUNCTION cusolverDnZgetrf_bufferSize(handle_ptr, m, n, &
                                                        A, lda, lwork) &
                   BIND(C, NAME='f_cusolverDnZgetrf_bufferSize')
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle_ptr
      INTEGER(C_INT), VALUE :: m, n
      TYPE(C_PTR), VALUE :: A
      INTEGER(C_INT), VALUE :: lda
      INTEGER(C_INT) :: lwork
    END FUNCTION cusolverDnZgetrf_bufferSize

    INTEGER(C_INT) FUNCTION cusolverDnZgetrf(handle_ptr, m, n, A, lda, &
                                             Workspace, devIpiv, devInfo) &
                   BIND(C, NAME='f_cusolverDnZgetrf')
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle_ptr
      INTEGER(C_INT), VALUE :: m, n
      TYPE(C_PTR), VALUE :: A
      INTEGER(C_INT), VALUE :: lda
      TYPE(C_PTR), VALUE :: Workspace
      TYPE(C_PTR), VALUE :: devIpiv
      INTEGER(C_INT) :: devInfo
    END FUNCTION cusolverDnZgetrf

    INTEGER(C_INT) FUNCTION cublasZtrsm(handle_ptr, side, uplo, transa, diag, &
                                        m, n, alpha, A, lda, B, ldb) &
                                        BIND(C, NAME='f_cublasZtrsm')
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle_ptr
      INTEGER(C_INT), VALUE :: side, uplo, transa, diag
      INTEGER(C_INT), VALUE :: m, n
      COMPLEX(C_DOUBLE_COMPLEX), VALUE :: alpha
      TYPE(C_PTR), VALUE :: A
      INTEGER(C_INT), VALUE :: lda
      TYPE(C_PTR), VALUE :: B
      INTEGER(C_INT), VALUE :: ldb
    END FUNCTION cublasZtrsm

    ! TODO: Add more interfaces to helper functions for multiple GPU support
    INTEGER(C_INT) FUNCTION cudaSetDevice(deviceidx) &
                   BIND(C, NAME='cudaSetDevice')
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE :: deviceidx
    END FUNCTION cudaSetDevice

!------------------------------------------------------------------------------

   END INTERFACE

end module cublas_f
