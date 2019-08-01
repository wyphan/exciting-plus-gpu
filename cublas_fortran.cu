#include <stdio.h>
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cusolverDn.h"

extern "C" int f_cublasCreate(cublasHandle_t **handle)
{
    *handle = (cublasHandle_t*)malloc(sizeof(cublasHandle_t));
    return cublasCreate(*handle);
}

extern "C" int f_cublasZgemm(cublasHandle_t *handle,
               cublasOperation_t transa, cublasOperation_t transb, 
              int m, int n, int k, 
              const cuDoubleComplex *alpha,
              const cuDoubleComplex *A, int lda, 
              const cuDoubleComplex *B, int ldb,
              const cuDoubleComplex *beta, 
              cuDoubleComplex *C, int ldc)
{
    return cublasZgemm(*handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}

extern "C" int f_cublasZgemmBatched(cublasHandle_t *handle,
               cublasOperation_t transa, cublasOperation_t transb,
              int m, int n, int k,
              const cuDoubleComplex *alpha,
              const cuDoubleComplex **A, int lda,
              const cuDoubleComplex **B, int ldb,
              const cuDoubleComplex *beta,
              cuDoubleComplex **C, int ldc,
              int batch_count)
{
    return cublasZgemmBatched(*handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,batch_count);
}

extern "C" void f_cublasDestroy(cublasHandle_t *handle)
{
    cublasDestroy(*handle);
    free(handle);
}

extern "C" int f_cudaStreamCreate(cudaStream_t **stream)
{
    *stream = (cudaStream_t *) malloc(sizeof(cudaStream_t));
    return cudaStreamCreate(*stream);
}

extern "C" int f_cublasSetStream(cublasHandle_t *handle, cudaStream_t *streamid)
{
    return cublasSetStream(*handle, *streamid);
}

extern "C" void f_cudaStreamDestroy(cudaStream_t *stream)
{
    cudaStreamDestroy(*stream);
}

extern "C" void f_addOffsetToPtr(size_t* array, size_t offset)
{
    array[0] = array[0] + offset;
}

extern "C" void f_printValue(void** address)
{
    //printf(" Write something out!!! PLEASE!\n\n\n");
    printf("     value is %zu\n", ( (size_t*)(*address) ) );
}

//-----------------------------------------------------------------------------

extern "C" int f_cusolverDnCreate(cusolverDnHandle_t **handle)
{
    *handle = (cusolverDnHandle_t*)malloc(sizeof(cusolverDnHandle_t));
    return cusolverDnCreate(*handle);
}

extern "C" void f_cusolverDnDestroy(cusolverDnHandle_t *handle)
{
    cusolverDnDestroy(*handle);
    free(handle);
}

extern "C" int f_cusolverDnSetStream(cusolverDnHandle_t *handle,
				     cudaStream_t *streamId)
{
    return cusolverDnSetStream(*handle, *streamId);
}

extern "C" int f_cusolverDnZgetrf_bufferSize(cusolverDnHandle_t *handle,
					     int m, int n,
					     cuDoubleComplex *A, int lda,
					     int *Lwork)
{
  return cusolverDnZgetrf_bufferSize(*handle, m, n, A, lda, Lwork);
}

extern "C" int f_cusolverDnZgetrf(cusolverDnHandle_t *handle,
				  int m, int n,
				  cuDoubleComplex *A, int lda,
				  cuDoubleComplex *Workspace,
				  int *devIpiv,
				  int *devInfo)
{
  return cusolverDnZgetrf(*handle, m, n, A, lda, Workspace, devIpiv, devInfo);
}

extern "C" int f_cublasZtrsm(cublasHandle_t *handle,
			     cublasSideMode_t side, cublasFillMode_t uplo,
			     cublasOperation_t transa, cublasDiagType_t diag,
			     int m, int n,
			     const cuDoubleComplex *alpha,
			     const cuDoubleComplex *A, int lda, 
			     cuDoubleComplex *B, int ldb)
{
  return cublasZtrsm(*handle, side, uplo, transa, diag, m, n, alpha, A, lda,
		     B, ldb);
}
