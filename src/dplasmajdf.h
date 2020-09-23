#ifndef _DPLASMAJDF_H_
#define _DPLASMAJDF_H_

#include "dplasma.h"
#include "dplasmaaux.h"
#include "cores/core_blas.h"
#include "cores/dplasma_cores.h"
#include "parsec/private_mempool.h"
#include "floputils.h"

/* sqrt function; these names needed for the precision generator stage */
#define dplasma_zsqrt csqrt
#define dplasma_csqrt csqrtf
#define dplasma_dsqrt sqrt
#define dplasma_ssqrt sqrtf

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)


#ifdef DPLASMA_TRACE_KERNELS
#   include <stdlib.h>
#   include <stdio.h>
#   define printlog(str, ...) fprintf(stderr, "thread %d VP %d " str "\n", \
                                      es->th_id, es->virtual_process->vp_id, __VA_ARGS__)
#   define printlogcuda(str, ...) fprintf(stderr, "cuda %d " str "\n", \
                                          gpu_device->cuda_index, __VA_ARGS__)
#else
#   define printlog(...) do {} while(0)
#   define printlogcuda(...) do {} while(0)
#endif

#ifndef PARSEC_HAVE_MPI
#define TEMP_TYPE MPITYPE
#undef MPITYPE
#define MPITYPE ((parsec_datatype_t)QUOTEME(TEMP_TYPE))
#undef TEMP_TYPE
#endif  /* PARSEC_HAVE_MPI */

#if defined(DPLASMA_HAVE_CUDA)
#include <cublas.h>

typedef cublasStatus_t (*cublas_zgemm_t) ( char TRANSA, char TRANSB, int m, int n, int k,
                                 cuDoubleComplex alpha, cuDoubleComplex *d_A, int lda,
                                 cuDoubleComplex *d_B, int ldb,
                                 cuDoubleComplex beta,  cuDoubleComplex *d_C, int ldc );
typedef cublasStatus_t (*cublas_zgemm_v2_t) (cublasHandle_t handle,
        cublasOperation_t transa, cublasOperation_t transb, int m, int n, int k,
        cuDoubleComplex *alpha, cuDoubleComplex *A, int lda,
        cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta, cuDoubleComplex *C, int ldc);
typedef cublasStatus_t (*cublas_cgemm_t) ( char TRANSA, char TRANSB, int m, int n, int k,
                                 cuFloatComplex alpha, cuFloatComplex *d_A, int lda,
                                 cuFloatComplex *d_B, int ldb,
                                 cuFloatComplex beta,  cuFloatComplex *d_C, int ldc );
typedef cublasStatus_t (*cublas_cgemm_v2_t) (cublasHandle_t handle,
        cublasOperation_t transa, cublasOperation_t transb, int m, int n, int k,
        cuComplex *alpha, cuComplex *A, int lda,
        cuComplex *B, int ldb,
        cuComplex *beta, cuComplex *C, int ldc);
typedef cublasStatus_t (*cublas_dgemm_t) ( char TRANSA, char TRANSB, int m, int n, int k,
                                 double alpha, double *d_A, int lda,
                                 double *d_B, int ldb,
                                 double beta,  double *d_C, int ldc );
typedef cublasStatus_t (*cublas_dgemm_v2_t) (cublasHandle_t handle,
        cublasOperation_t transa, cublasOperation_t transb, int m, int n, int k,
        double *alpha, double *A, int lda,
        double *B, int ldb,
        double *beta, double *C, int ldc);
typedef cublasStatus_t (*cublas_sgemm_t) ( char TRANSA, char TRANSB, int m, int n, int k,
                                 float alpha, float *d_A, int lda,
                                              float *d_B, int ldb,
                                 float beta,  float *d_C, int ldc );
typedef cublasStatus_t (*cublas_sgemm_v2_t) (cublasHandle_t handle,
        cublasOperation_t transa, cublasOperation_t transb, int m, int n, int k,
        float *alpha, float *A, int lda,
        float *B, int ldb,
        float *beta, float *C, int ldc);

typedef cublasStatus_t (*cublas_strsm_v2_t) (cublasHandle_t handle,
                                   cublasSideMode_t side,
                                   cublasFillMode_t uplo,
                                   cublasOperation_t trans,
                                   cublasDiagType_t diag,
                                   int m, int n, float *alpha,
                                   float *A, int lda,
                                   float *B, int ldb);
typedef cublasStatus_t (*cublas_dtrsm_v2_t) (cublasHandle_t handle,
                                   cublasSideMode_t side,
                                   cublasFillMode_t uplo,
                                   cublasOperation_t trans,
                                   cublasDiagType_t diag,
                                   int m, int n, double *alpha,
                                   double *A, int lda,
                                   double *B, int ldb);
typedef cublasStatus_t (*cublas_ctrsm_v2_t) (cublasHandle_t handle,
                                   cublasSideMode_t side,
                                   cublasFillMode_t uplo,
                                   cublasOperation_t trans,
                                   cublasDiagType_t diag,
                                   int m, int n, cuFloatComplex *alpha,
                                   cuFloatComplex *A, int lda,
                                   cuFloatComplex *B, int ldb);
typedef cublasStatus_t (*cublas_ztrsm_v2_t) (cublasHandle_t handle,
                                   cublasSideMode_t side,
                                   cublasFillMode_t uplo,
                                   cublasOperation_t trans,
                                   cublasDiagType_t diag,
                                   int m, int n, cuDoubleComplex *alpha,
                                   cuDoubleComplex *A, int lda,
                                   cuDoubleComplex *B, int ldb);

typedef cublasStatus_t (*cublas_cherk_v2_t) (cublasHandle_t handle,
                                   cublasFillMode_t uplo,
                                   cublasOperation_t trans, int n, int k,
                                   float *alpha, cuFloatComplex *A, int lda,
                                   float *beta, cuFloatComplex *C, int ldc);
typedef cublasStatus_t (*cublas_zherk_v2_t) (cublasHandle_t handle,
                                   cublasFillMode_t uplo,
                                   cublasOperation_t trans, int n, int k,
                                   double *alpha, cuDoubleComplex *A, int lda,
                                   double *beta, cuDoubleComplex *C, int ldc);
typedef cublasStatus_t (*cublas_ssyrk_v2_t) (cublasHandle_t handle,
                                   cublasFillMode_t uplo,
                                   cublasOperation_t trans, int n, int k,
                                   float *alpha, float *A, int lda,
                                   float *beta, float *C, int ldc);
typedef cublasStatus_t (*cublas_dsyrk_v2_t) (cublasHandle_t handle,
                                   cublasFillMode_t uplo,
                                   cublasOperation_t trans, int n, int k,
                                   double *alpha, double *A, int lda,
                                   double *beta, double *C, int ldc);
#endif  /* defined(DPLASMA_HAVE_CUDA) */

#endif /* _DPLASMAJDF_H_ */

