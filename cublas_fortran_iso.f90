module cublas_f
  use ISO_C_BINDING

    enum, BIND(C)
        enumerator :: CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C
    end enum
   
    enum, BIND(C)
        enumerator :: cudaMemcpyHostToHost, cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost, &
                      cudaMemcpyDeviceToDevice, cudaMemcpyDefault
    end enum

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

  END INTERFACE

end module cublas_f
