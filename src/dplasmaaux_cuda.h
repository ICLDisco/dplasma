/*
 * Copyright (c) 2023-     The University of Tennessee and The University
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
 * when cublas_v2 is needed by including the header before dplasmaaux.h. Otherwise, we will include
 * cublas.h (v1) automatically if CUDA is enabled.
 */
#if !defined(CUBLAS_V2_H_)
#include <cublas.h>
#endif  /* !defined(CUBLAS_V2_H_) */

#define dplasma_cublas_side(side)                                         \
    assert( (side == dplasmaRight) || (side == dplasmaLeft) );            \
    side = (side == dplasmaRight) ? CUBLAS_SIDE_RIGHT : CUBLAS_SIDE_LEFT;


#define dplasma_cublas_diag(diag)                                              \
    assert( (diag == dplasmaNonUnit) || (diag == dplasmaUnit) );               \
    diag = (diag == dplasmaNonUnit) ? CUBLAS_DIAG_NON_UNIT : CUBLAS_DIAG_UNIT;

#define dplasma_cublas_fill(fill)                                                    \
    assert( (fill == dplasmaLower) || (fill == dplasmaUpper) );                      \
    fill = (fill == dplasmaLower) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;

#if defined(PRECISION_z) || defined(PRECISION_c)
#define dplasma_cublas_op(trans)                 \
    assert( (trans == dplasmaNoTrans) || (trans == dplasmaTrans) || (trans == dplasmaConjTrans) ); \
    switch(trans){                               \
        case dplasmaNoTrans:                     \
            trans = CUBLAS_OP_N;                 \
            break;                               \
        case dplasmaTrans:                       \
            trans = CUBLAS_OP_T;                 \
            break;                               \
        case dplasmaConjTrans:                   \
            trans = CUBLAS_OP_C;                 \
            break;                               \
        default:                                 \
            trans = CUBLAS_OP_N;                 \
            break;                               \
    }
#else
#define dplasma_cublas_op(trans)                                    \
    assert( (trans == dplasmaNoTrans) || (trans == dplasmaTrans) ); \
    trans = (trans == dplasmaNoTrans) ? CUBLAS_OP_N : CUBLAS_OP_T;
#endif /* PRECISION_z || PRECISION_c */

extern parsec_info_id_t CuHI;
extern parsec_info_id_t WoSI;

typedef struct {
    cublasHandle_t cublas_handle;
    void * cusolverDn_handle;
} dplasma_cuda_handles_t;

void *dplasma_create_cuda_handles(void *obj, void *user);

char *dplasma_cublas_error_to_string(cublasStatus_t cublas_status);

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

char *dplasma_cusolver_error_to_string(cusolverStatus_t cusolver_status);

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
#endif  /* defined(DPLASMA_HAVE_CUDA) */

#endif /* __DPLAMAAUX_CUDA_H__ */
