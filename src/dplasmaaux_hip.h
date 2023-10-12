/*
 * Copyright (c) 2023-     The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * $COPYRIGHT
 *
 */

#ifndef _DPLASMAAAUX_HIP_H_
#define _DPLASMAAAUX_HIP_H_

#if defined(DPLASMA_HAVE_HIP)
#include "parsec/mca/device/hip/device_hip.h"

/**
 * DPLASMA currently supports a mix of hipblas v1 and v2, but not in the same source file. Thus,
 * the simplest way to provide common headers is to require the developer to manually specify
 * when hipblas_v2 is needed by including the header before dplasmaaux.h. Otherwise, we will include
 * hipblas.h (v1) automatically if HIP is enabled.
 */
#if !defined(HIPBLAS_V2_H_)
#include <hipblas/hipblas.h>
#endif  /* !defined(HIPBLAS_V2_H_) */

#define dplasma_hipblas_side(side)                                         \
    assert( (side == dplasmaRight) || (side == dplasmaLeft) );            \
    side = (side == dplasmaRight) ? HIPBLAS_SIDE_RIGHT : HIPBLAS_SIDE_LEFT;


#define dplasma_hipblas_diag(diag)                                              \
    assert( (diag == dplasmaNonUnit) || (diag == dplasmaUnit) );               \
    diag = (diag == dplasmaNonUnit) ? HIPBLAS_DIAG_NON_UNIT : HIPBLAS_DIAG_UNIT;

#define dplasma_hipblas_fill(fill)                                                    \
    assert( (fill == dplasmaLower) || (fill == dplasmaUpper) );                      \
    fill = (fill == dplasmaLower) ? HIPBLAS_FILL_MODE_LOWER : HIPBLAS_FILL_MODE_UPPER;

#if defined(PRECISION_z) || defined(PRECISION_c)
#define dplasma_hipblas_op(trans)                 \
    assert( (trans == dplasmaNoTrans) || (trans == dplasmaTrans) || (trans == dplasmaConjTrans) ); \
    switch(trans){                               \
        case dplasmaNoTrans:                     \
            trans = HIPBLAS_OP_N;                 \
            break;                               \
        case dplasmaTrans:                       \
            trans = HIPBLAS_OP_T;                 \
            break;                               \
        case dplasmaConjTrans:                   \
            trans = HIPBLAS_OP_C;                 \
            break;                               \
        default:                                 \
            trans = HIPBLAS_OP_N;                 \
            break;                               \
    }
#else
#define dplasma_hipblas_op(trans)                                    \
    assert( (trans == dplasmaNoTrans) || (trans == dplasmaTrans) ); \
    trans = (trans == dplasmaNoTrans) ? HIPBLAS_OP_N : HIPBLAS_OP_T;
#endif /* PRECISION_z || PRECISION_c */

extern parsec_info_id_t CuHI;
extern parsec_info_id_t WoSI;

typedef struct {
    hipblasHandle_t hipblas_handle;
} dplasma_hip_handles_t;

void *dplasma_create_hip_handles(void *obj, void *user);

#define DPLASMA_ROCBLAS_CHECK_ERROR(STR, ERROR, CODE) \
    do { \
        rocblas_status __error = (rocblas_status) (ERROR); \
        if(rocblas_status_success != __error) { \
            parsec_warning( "%s:%d %s%s", __FILE__, __LINE__, \
                            (STR), rocblas_status_to_string(__error)); \
            CODE; \
        } \
    } while(0)

/* For some reason the error values are not the same... */
#define DPLASMA_HIPBLAS_CHECK_ERROR(STR, ERROR, CODE) \
    do { \
        hipblasStatus_t __error = (hipblasStatus_t) (ERROR); \
        if(HIPBLAS_STATUS_SUCCESS != __error) { \
            parsec_warning( "%s:%d %s%s", __FILE__, __LINE__, \
                            (STR), hipblasStatusToString(__error)); \
            CODE; \
        } \
    } while(0)

#endif /* defined(DPLASMA_HAVE_HIP */
#endif /* __DPLAMAAUX_HIP_H__ */