/*
 * Copyright (c) 2023-     The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#ifndef __DTD_WRAPPERS_Z_H__
#define __DTD_WRAPPERS_Z_H__

#include "dplasma/types.h"
#include "parsec/interfaces/dtd/insert_function.h"
#include "cores/core_blas.h"

#if defined(DPLASMA_HAVE_CUDA)
#include "parsec/execution_stream.h"
#include "parsec/parsec_internal.h"
#include "parsec/mca/device/cuda/device_cuda.h"
#include "parsec/utils/zone_malloc.h"
#include "dplasmaaux.h"
#include "potrf_cublas_utils.h"
#include <cublas_v2.h>

/* probably need to add this to substitions */
#if defined(PRECISION_s)
#define CUSOLVER_COMPUTE_TYPE CUDA_R_32F
#elif defined(PRECISION_d)
#define CUSOLVER_COMPUTE_TYPE CUDA_R_64F
#elif defined(PRECISION_c)
#define CUSOLVER_COMPUTE_TYPE CUDA_C_32F
#elif defined(PRECISION_z)
#define CUSOLVER_COMPUTE_TYPE CUDA_C_64F
#endif

typedef struct zpotrf_dtd_workspace_info_s {
    int mb;
    int nb;
    dplasma_enum_t uplo;
} zpotrf_dtd_workspace_info_t;

void* zpotrf_dtd_create_workspace(void *obj, void *user);
void zpotrf_dtd_destroy_workspace(void *_ws, void *_n);

void* zpotrf_dtd_create_params(void *obj, void *user);
void zpotrf_dtd_destroy_params(void *params, void *_n);

#endif

parsec_task_class_t * parsec_dtd_create_zpotrf_task_class(parsec_taskpool_t * dtd_tp, int tile_full, int devices);
parsec_task_class_t * parsec_dtd_create_ztrsm_task_class(parsec_taskpool_t * dtd_tp, int tile_full, int devices);
parsec_task_class_t * parsec_dtd_create_zherk_task_class(parsec_taskpool_t * dtd_tp, int tile_full, int devices);
parsec_task_class_t * parsec_dtd_create_zgemm_task_class(parsec_taskpool_t * dtd_tp, int tile_full, int devices);

#endif /* __DTD_WRAPPERS_Z_H__ */