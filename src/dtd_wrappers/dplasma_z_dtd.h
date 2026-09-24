/*
 * Copyright (c) 2023-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2026      NVIDIA Corporation.  All rights reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#ifndef __DTD_WRAPPERS_Z_H__
#define __DTD_WRAPPERS_Z_H__

#include "dplasma/types.h"
#include "parsec/interfaces/dtd/insert_function.h"
#include "cores/core_blas.h"

#if defined(DPLASMA_HAVE_CUDA) || defined(DPLASMA_HAVE_HIP)
#include "parsec/execution_stream.h"
#include "parsec/parsec_internal.h"
#include "parsec/utils/zone_malloc.h"
/* Pulls in dplasmaaux_cuda.h and/or dplasmaaux_hip.h as available. */
#include "dplasmaaux.h"
#include "potrf_gpu_workspaces.h"
#endif

#if defined(DPLASMA_HAVE_CUDA)
#include "parsec/mca/device/cuda/device_cuda.h"

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

void* zpotrf_dtd_create_workspace(void *obj, void *user);
void zpotrf_dtd_destroy_workspace(void *_ws, void *_n);

void* zpotrf_dtd_create_params(void *obj, void *user);
void zpotrf_dtd_destroy_params(void *params, void *_n);

#endif /* defined(DPLASMA_HAVE_CUDA) */

#if defined(DPLASMA_HAVE_HIP)
#include "parsec/mca/device/hip/device_hip.h"

/* rocSOLVER sizes its own scratch, so the workspace only carries the device
 * info word and none of the mb/nb/uplo the cuSOLVER query needs. */
void* zpotrf_dtd_create_hip_workspace(void *obj, void *user);
void zpotrf_dtd_destroy_hip_workspace(void *_ws, void *_n);

#endif /* defined(DPLASMA_HAVE_HIP) */

#if defined(DPLASMA_HAVE_CUDA) || defined(DPLASMA_HAVE_HIP)
typedef struct zpotrf_dtd_workspace_info_s {
    int mb;
    int nb;
    dplasma_enum_t uplo;
} zpotrf_dtd_workspace_info_t;
#endif

parsec_task_class_t * parsec_dtd_create_zpotrf_task_class(parsec_taskpool_t * dtd_tp, int tile_full, int devices);
parsec_task_class_t * parsec_dtd_create_ztrsm_task_class(parsec_taskpool_t * dtd_tp, int tile_full, int devices);
parsec_task_class_t * parsec_dtd_create_zherk_task_class(parsec_taskpool_t * dtd_tp, int tile_full, int devices);
parsec_task_class_t * parsec_dtd_create_zgemm_task_class(parsec_taskpool_t * dtd_tp, int tile_full, int devices);

#endif /* __DTD_WRAPPERS_Z_H__ */
