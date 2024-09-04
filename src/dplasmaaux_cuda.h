/*
 * Copyright (c) 2023-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * $COPYRIGHT
 *
 */

#ifndef _DPLASMAAAUX_CUDA_H_
#define _DPLASMAAAUX_CUDA_H_


#if defined(DPLASMA_HAVE_CUDA)
#include "parsec/mca/device/cuda/device_cuda.h"

/**
 * DPLASMA currently supports a mix of cublas v1 and v2, but not in the same source file. Thus,
 * the simplest way to provide common headers is to require the developer to manually specify
 * when legacy cublas is needed by including the header before dplasmaaux.h. Otherwise, we will include
 * cublas_v2.h (v2) automatically if CUDA is enabled.
 */
#if !defined(CUBLAS_H_)
#include <cublas_v2.h>
#include "dplasma/constants.h"

static inline cublasSideMode_t dplasma_cublas_side(int side) {
    assert( (side == dplasmaRight) || (side == dplasmaLeft) );
    return (side == dplasmaRight) ? CUBLAS_SIDE_RIGHT : CUBLAS_SIDE_LEFT;
}

static inline cublasDiagType_t dplasma_cublas_diag(int diag) {
    assert( (diag == dplasmaNonUnit) || (diag == dplasmaUnit) );
    return (diag == dplasmaNonUnit) ? CUBLAS_DIAG_NON_UNIT : CUBLAS_DIAG_UNIT;
}

static inline cublasFillMode_t dplasma_cublas_fill(int fill) {
    assert( (fill == dplasmaLower) || (fill == dplasmaUpper) );
    return (fill == dplasmaLower) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
}

static inline cublasOperation_t dplasma_cublas_op(int trans) {
#if defined(PRECISION_d) || defined(PRECISION_s)
    assert( (trans == dplasmaNoTrans) || (trans == dplasmaTrans) );
#else
    assert( (trans == dplasmaNoTrans) || (trans == dplasmaTrans) || (trans == dplasmaConjTrans) );
#endif /* PRECISION_d || PRECISION_s */
    return (trans == dplasmaConjTrans) ? CUBLAS_OP_C: ((trans == dplasmaTrans) ? CUBLAS_OP_T : CUBLAS_OP_N);
}
#endif  /* !defined(CUBLAS_V2_H_) */

extern parsec_info_id_t dplasma_dtd_cuda_infoid;
extern parsec_info_id_t dplasma_dtd_cuda_workspace_infoid;

typedef struct {
    cublasHandle_t cublas_handle;
    void * cusolverDn_handle;
} dplasma_cuda_handles_t;

void *dplasma_create_cuda_handles(void *obj, void *user);
void dplasma_destroy_cuda_handles(void *_h, void *_n);

const char *dplasma_cublas_error_to_string(cublasStatus_t cublas_status);

#define DPLASMA_CUBLAS_CHECK_STATUS( STR, STATUS, CODE )                     \
    do {                                                                     \
        cublasStatus_t __cublas_status = (cublasStatus_t) (STATUS);          \
        if( CUBLAS_STATUS_SUCCESS != __cublas_status ) {                     \
            parsec_warning( "%s:%d %s%s", __FILE__, __LINE__,                \
                            (STR), dplasma_cublas_error_to_string(__cublas_status) ); \
            CODE;                                                            \
        }                                                                    \
    } while(0)

#if defined(CUBLAS_V2_H_)
/* Support for cusolve requires cublas_v2 */
#include <cusolverDn.h>

const char *dplasma_cusolver_error_to_string(cusolverStatus_t cusolver_status);

#define DPLASMA_CUSOLVER_CHECK_STATUS( STR, STATUS, CODE )                                \
    do {                                                                                  \
        cusolverStatus_t __cusolver_status = (cusolverStatus_t) (STATUS);                 \
        if( CUSOLVER_STATUS_SUCCESS != __cusolver_status ) {                              \
            parsec_warning( "%s:%d %s%s", __FILE__, __LINE__,                             \
                            (STR), dplasma_cusolver_error_to_string(__cusolver_status) ); \
            CODE;                                                                         \
        }                                                                                 \
    } while(0)
#endif  /* defined(CUBLAS_V2_H_) */

#else
#warning "DPLASMA_HAVE_CUDA not defined, this file should not be included then."
#endif  /* defined(DPLASMA_HAVE_CUDA) */
#endif /* __DPLAMAAUX_CUDA_H__ */
