/*
 * Copyright (c) 2023-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2026      NVIDIA Corporation.  All rights reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "dplasma/config.h"
#include "dplasma_z_dtd.h"

int
parsec_core_zpotrf(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int uplo;
    int m, lda, *info;
    dplasma_complex64_t *A;

    parsec_dtd_unpack_args(this_task, &uplo, &m, &A, &lda, &info);

    CORE_zpotrf(uplo, m, A, lda, info);

    return PARSEC_HOOK_RETURN_DONE;
}

#if defined(DPLASMA_HAVE_CUDA)

void*
zpotrf_dtd_create_workspace(void *obj, void *user)
{
    parsec_device_gpu_module_t *gpu_device = (parsec_device_gpu_module_t *)obj;
    cusolverDnHandle_t cusolverDnHandle;
    cusolverStatus_t status;
    zpotrf_dtd_workspace_info_t *infos = (zpotrf_dtd_workspace_info_t*) user;
    dplasma_potrf_gpu_workspaces_t *wp = NULL;
    void *tmpmem;
    size_t workspace_size;
    size_t host_size;
    int mb = infos->mb;
    int nb = infos->nb;
    size_t elt_size = sizeof(cuDoubleComplex);
    dplasma_enum_t uplo = infos->uplo;
    cusolverDnParams_t* params;

    status = cusolverDnCreate(&cusolverDnHandle);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    (void)status;


    params = malloc(sizeof(cusolverDnParams_t));

    status = cusolverDnCreateParams(params);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    status = cusolverDnXpotrf_bufferSize(cusolverDnHandle, *params, dplasma_cublas_fill(uplo), nb, CUSOLVER_COMPUTE_TYPE, NULL, mb, CUSOLVER_COMPUTE_TYPE, &workspace_size, &host_size);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    cusolverDnDestroy(cusolverDnHandle);

    /* As in the PTG zpotrf constructor, the scratch bypasses the zone PaRSEC manages
     * for tiles so that a zone saturated by data copies cannot stall potrf. */
    if( PARSEC_SUCCESS != gpu_device->memory_allocate(gpu_device,
                                                      workspace_size * elt_size + sizeof(int),
                                                      &tmpmem) ) {
        cusolverDnDestroyParams(*params);
        free(params);
        return NULL;
    }

    wp = (dplasma_potrf_gpu_workspaces_t*)malloc(sizeof(dplasma_potrf_gpu_workspaces_t));
    wp->tmpmem = tmpmem;
    wp->lwork = workspace_size;
    wp->gpu_device = gpu_device;
    wp->params = params;
    wp->host_size = host_size;
    wp->host_buffer = malloc(host_size*sizeof(dplasma_complex64_t));

    return wp;
}

void
zpotrf_dtd_destroy_workspace(void *_ws, void *_n)
{
    dplasma_potrf_gpu_workspaces_t *ws = (dplasma_potrf_gpu_workspaces_t*)_ws;
    parsec_device_gpu_module_t *gpu_device = (parsec_device_gpu_module_t*)ws->gpu_device;
    cusolverDnParams_t* params = ws->params;
    cusolverStatus_t status = cusolverDnDestroyParams(*params);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    gpu_device->memory_free(gpu_device, ws->tmpmem);
    free(params);
    free(ws->host_buffer);
    free(ws);
    (void)_n;
    (void)params;
    (void)status;
}

int
parsec_core_zpotrf_cuda(parsec_device_gpu_module_t* gpu_device,
                        parsec_gpu_task_t*          gpu_task,
                        parsec_gpu_exec_stream_t*   gpu_stream)
{
    int uplo;
    int m, lda, *info;
    dplasma_complex64_t *A, *Ag;
    parsec_task_t* this_task = gpu_task->ec;
    cusolverStatus_t status;
    dplasma_cuda_handles_t* handles;
    dplasma_potrf_gpu_workspaces_t *wp;
    cuDoubleComplex *workspace;
    cusolverDnParams_t* params;
    int             *d_iinfo;

    parsec_dtd_unpack_args(this_task, &uplo, &m, &A, &lda, &info);

    Ag = parsec_dtd_get_dev_ptr(this_task, 0);

    handles = parsec_info_get(&gpu_stream->infos, dplasma_dtd_cuda_infoid);
    assert(NULL != handles);
    wp = parsec_info_get(&gpu_device->super.infos, dplasma_dtd_cuda_workspace_infoid);
    /* The scratch allocation failed; retry once device memory frees up.
     * parsec_info_get does not cache a NULL, so the constructor runs again. */
    if( NULL == wp ) return PARSEC_HOOK_RETURN_AGAIN;

    workspace = (cuDoubleComplex*)wp->tmpmem;
    d_iinfo   = (int*)(wp->tmpmem + wp->lwork * sizeof(cuDoubleComplex));
    params = (cusolverDnParams_t*)wp->params;

#if defined(PARSEC_DEBUG_NOISIER)
    {
        char tmp[MAX_TASK_STRLEN];
        PARSEC_DEBUG_VERBOSE(10, parsec_gpu_output_stream, "GPU[%1d]:\tEnqueue on device %s priority %d",
                             gpu_device->super.device_index, parsec_task_snprintf(tmp, MAX_TASK_STRLEN,
                             (parsec_task_t *) this_task), this_task->priority);
    }
#endif /* defined(PARSEC_DEBUG_NOISIER) */

    status = cusolverDnXpotrf(handles->cusolverDn_handle, *params,
                              dplasma_cublas_fill(uplo), m, CUSOLVER_COMPUTE_TYPE,
                              (cuDoubleComplex *)Ag, lda, CUSOLVER_COMPUTE_TYPE,
                              workspace, wp->lwork, wp->host_buffer, wp->host_size, d_iinfo);

    DPLASMA_CUSOLVER_CHECK_STATUS( "cusolverDnZpotrf ", status,
                                   {return PARSEC_HOOK_RETURN_ERROR;} );

    (void)gpu_device;
    return PARSEC_HOOK_RETURN_DONE;
}

#endif /* DPLASMA_HAVE_CUDA */

#if defined(DPLASMA_HAVE_HIP)

void*
zpotrf_dtd_create_hip_workspace(void *obj, void *user)
{
    parsec_device_gpu_module_t *gpu_device = (parsec_device_gpu_module_t *)obj;
    dplasma_potrf_gpu_workspaces_t *wp = NULL;
    void *tmpmem;
    (void)user;

    /* As in the CUDA constructor, the scratch bypasses the zone PaRSEC manages for
     * tiles. rocSOLVER sizes its own scratch, so this only holds the info word. */
    if( PARSEC_SUCCESS != gpu_device->memory_allocate(gpu_device, sizeof(int), &tmpmem) )
        return NULL;

    wp = (dplasma_potrf_gpu_workspaces_t*)malloc(sizeof(dplasma_potrf_gpu_workspaces_t));
    wp->tmpmem = tmpmem;
    wp->lwork = 0;
    wp->gpu_device = gpu_device;
    wp->params = NULL;
    wp->host_size = 0;
    wp->host_buffer = NULL;

    return wp;
}

void
zpotrf_dtd_destroy_hip_workspace(void *_ws, void *_n)
{
    dplasma_potrf_gpu_workspaces_t *ws = (dplasma_potrf_gpu_workspaces_t*)_ws;
    parsec_device_gpu_module_t *gpu_device = (parsec_device_gpu_module_t*)ws->gpu_device;
    gpu_device->memory_free(gpu_device, ws->tmpmem);
    free(ws);
    (void)_n;
}

int
parsec_core_zpotrf_hip(parsec_device_gpu_module_t* gpu_device,
                       parsec_gpu_task_t*          gpu_task,
                       parsec_gpu_exec_stream_t*   gpu_stream)
{
    int uplo;
    int m, lda, *info;
    dplasma_complex64_t *A, *Ag;
    parsec_task_t* this_task = gpu_task->ec;
    rocblas_status status;
    rocblas_fill rocblas_uplo;
    dplasma_hip_handles_t* handles;
    dplasma_potrf_gpu_workspaces_t *wp;
    int *d_iinfo;

    parsec_dtd_unpack_args(this_task, &uplo, &m, &A, &lda, &info);

    Ag = parsec_dtd_get_dev_ptr(this_task, 0);

    /* hipblasFillMode_t and rocblas_fill are distinct enums, so this cannot go
     * through dplasma_hipblas_fill. */
    if( dplasmaLower == uplo )
        rocblas_uplo = rocblas_fill_lower;
    else
        rocblas_uplo = rocblas_fill_upper;

    handles = parsec_info_get(&gpu_stream->infos, dplasma_dtd_hip_infoid);
    assert(NULL != handles);
    wp = parsec_info_get(&gpu_device->super.infos, dplasma_dtd_hip_workspace_infoid);
    /* The scratch allocation failed; retry once device memory frees up.
     * parsec_info_get does not cache a NULL, so the constructor runs again. */
    if( NULL == wp ) return PARSEC_HOOK_RETURN_AGAIN;

    d_iinfo = (int*)wp->tmpmem;

#if defined(PARSEC_DEBUG_NOISIER)
    {
        char tmp[MAX_TASK_STRLEN];
        PARSEC_DEBUG_VERBOSE(10, parsec_gpu_output_stream, "GPU[%1d]:\tEnqueue on device %s priority %d",
                             gpu_device->super.device_index, parsec_task_snprintf(tmp, MAX_TASK_STRLEN,
                             (parsec_task_t *) this_task), this_task->priority);
    }
#endif /* defined(PARSEC_DEBUG_NOISIER) */

    parsec_hip_exec_stream_t* hip_stream = (parsec_hip_exec_stream_t*)gpu_stream;
    hipblasSetStream( handles->hipblas_handle, hip_stream->hip_stream );
    status = rocsolver_zpotrf( handles->hipblas_handle, rocblas_uplo, m, Ag, lda, d_iinfo );

    DPLASMA_ROCBLAS_CHECK_ERROR( "rocsolver_zpotrf ", status,
                                 {return PARSEC_HOOK_RETURN_ERROR;} );

    (void)gpu_device;
    return PARSEC_HOOK_RETURN_DONE;
}

#endif /* DPLASMA_HAVE_HIP */

parsec_task_class_t*
parsec_dtd_create_zpotrf_task_class(parsec_taskpool_t* dtd_tp, int tile_full, int devices)
{
    parsec_task_class_t* zpotrf_tc =  parsec_dtd_create_task_class( dtd_tp, "zpotrf",
                                      /* uplo    */   sizeof(int),   PARSEC_VALUE,
                                      /* m       */   sizeof(int),   PARSEC_VALUE,
                                      /* A(k, k) */   PASSED_BY_REF, PARSEC_INOUT | tile_full | PARSEC_AFFINITY,
                                      /* ldak    */   sizeof(int),   PARSEC_VALUE,
                                      /* info    */   sizeof(int *), PARSEC_SCRATCH,
                                      PARSEC_DTD_ARG_END );

#if defined(DPLASMA_HAVE_CUDA)
    if( devices & PARSEC_DEV_CUDA )
        parsec_dtd_task_class_add_chore(dtd_tp, zpotrf_tc, PARSEC_DEV_CUDA, parsec_core_zpotrf_cuda);
#endif

#if defined(DPLASMA_HAVE_HIP)
    if( devices & PARSEC_DEV_HIP )
        parsec_dtd_task_class_add_chore(dtd_tp, zpotrf_tc, PARSEC_DEV_HIP, parsec_core_zpotrf_hip);
#endif

    if( devices & PARSEC_DEV_CPU )
        parsec_dtd_task_class_add_chore(dtd_tp, zpotrf_tc, PARSEC_DEV_CPU, parsec_core_zpotrf);

    return zpotrf_tc;
}
