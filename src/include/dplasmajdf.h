#ifndef _DPLASMAJDF_H_
#define _DPLASMAJDF_H_

#include "dplasma.h"
#include <core_blas.h>
#include "parsec/private_mempool.h"
#include "floputils.h"

/* sqrt function */
#define dplasma_zsqrt csqrt
#define dplasma_csqrt csqrtf
#define dplasma_dsqrt sqrt
#define dplasma_ssqrt sqrtf

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

#define plasma_const( x )  plasma_lapack_constants[x]

#ifdef PARSEC_CALL_TRACE
#   include <stdlib.h>
#   include <stdio.h>
#   define printlog(str, ...) fprintf(stderr, "thread %d VP %d " str "\n", \
                                      es->th_id, es->virtual_process->vp_id, __VA_ARGS__)
#   define printlogcuda(str, ...) fprintf(stderr, "cuda %d " str "\n", \
                                          gpu_device->cuda_index, __VA_ARGS__)
#   define OUTPUT(ARG)  printf ARG
#else
#   define printlog(...) do {} while(0)
#   define printlogcuda(...) do {} while(0)
#   define OUTPUT(ARG)
#endif

#ifndef PARSEC_HAVE_MPI
#define TEMP_TYPE MPITYPE
#undef MPITYPE
#define MPITYPE ((parsec_datatype_t)QUOTEME(TEMP_TYPE))
#undef TEMP_TYPE
#endif  /* PARSEC_HAVE_MPI */


#if defined(PARSEC_HAVE_CUDA)
#include <cublas.h>

typedef cublasStatus_t (*cublas_zgemm_t) ( char TRANSA, char TRANSB, int m, int n, int k,
                                 cuDoubleComplex alpha, cuDoubleComplex *d_A, int lda,
                                 cuDoubleComplex *d_B, int ldb,
                                 cuDoubleComplex beta,  cuDoubleComplex *d_C, int ldc );
typedef cublasStatus_t (*cublas_cgemm_t) ( char TRANSA, char TRANSB, int m, int n, int k,
                                 cuComplex alpha, cuComplex *d_A, int lda,
                                 cuComplex *d_B, int ldb,
                                 cuComplex beta,  cuComplex *d_C, int ldc );
typedef cublasStatus_t (*cublas_dgemm_t) ( char TRANSA, char TRANSB, int m, int n, int k,
                                 double alpha, double *d_A, int lda,
                                 double *d_B, int ldb,
                                 double beta,  double *d_C, int ldc );
typedef cublasStatus_t (*cublas_sgemm_t) ( char TRANSA, char TRANSB, int m, int n, int k,
                                 float alpha, float *d_A, int lda,
                                              float *d_B, int ldb,
                                 float beta,  float *d_C, int ldc );
#endif  /* defined(PARSEC_HAVE_CUDA) */

#endif /* _DPLASMAJDF_H_ */

