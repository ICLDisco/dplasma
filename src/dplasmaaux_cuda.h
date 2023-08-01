/*
 * Copyright (c) 2023-     The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * $COPYRIGHT
 *
 */

#ifndef _DPLASMAAAUX_CUDA_H_
#define _DPLASMAAAUX_CUDA_H_

#include "parsec/mca/device/cuda/device_cuda.h"

extern parsec_info_id_t CuHI;
extern parsec_info_id_t WoSI;

typedef struct {
    cublasHandle_t cublas_handle;
    void * cusolverDn_handle;
} dplasma_cuda_handles_t;

void *dplasma_create_cuda_handles(void *obj, void *user);

#define DPLASMA_CUBLAS_CHECK_STATUS( STR, STATUS, CODE )                     \
    do {                                                                     \
        cublasStatus_t __cublas_status = (cublasStatus_t) (STATUS);          \
        if( CUBLAS_STATUS_SUCCESS != __cublas_status ) {                     \
            parsec_warning( "%s:%d %s%s", __FILE__, __LINE__,                \
                            (STR), cublasGetStatusString(__cublas_status) ); \
            CODE;                                                            \
        }                                                                    \
    } while(0)

#define DPLASMA_CUSOLVER_CHECK_STATUS( STR, STATUS, CODE )                                \
    do {                                                                                  \
        cusolverStatus_t __cusolver_status = (cusolverStatus_t) (STATUS);                 \
        if( CUSOLVER_STATUS_SUCCESS != __cusolver_status ) {                              \
            parsec_warning( "%s:%d %s%s", __FILE__, __LINE__,                             \
                            (STR), dplasma_cusolver_error_to_string(__cusolver_status) ); \
            CODE;                                                                         \
        }                                                                                 \
    } while(0)

#endif /* __DPLAMAAUX_CUDA_H__ */