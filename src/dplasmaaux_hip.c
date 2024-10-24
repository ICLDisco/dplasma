/*
 * Copyright (c) 2023-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 */
#include "dplasma/config.h"
#include "parsec/utils/zone_malloc.h"
#include "parsec/utils/show_help.h"
#include "dplasmaaux_hip.h"
#include "potrf_gpu_workspaces.h"

#include <hipblas/hipblas.h>
#include <hipsolver/hipsolver.h>

/*
 * Global info ID's for cublas handles and workspaces
 * Should be initialized in the tests
 * with the return of parsec_info_register
 * or parsec_info_lookup
 */
parsec_info_id_t dplasma_dtd_hip_infoid = -1;

/* Unfortunately, hipSolver does not provide a error to string function */
const char *dplasma_hipsolver_error_to_string(hipsolverStatus_t hipsolver_status)
{
    switch(hipsolver_status) {
        case HIPSOLVER_STATUS_SUCCESS: return "HIPSOLVER_STATUS_SUCCESS";
        case HIPSOLVER_STATUS_NOT_INITIALIZED: return "HIPSOLVER_STATUS_NOT_INITIALIZED";
        case HIPSOLVER_STATUS_ALLOC_FAILED: return "HIPSOLVER_STATUS_ALLOC_FAILED";
        case HIPSOLVER_STATUS_INVALID_VALUE: return "HIPSOLVER_STATUS_INVALID_VALUE";
        case HIPSOLVER_STATUS_ARCH_MISMATCH: return "HIPSOLVER_STATUS_ARCH_MISMATCH";
        case HIPSOLVER_STATUS_EXECUTION_FAILED: return "HIPSOLVER_STATUS_EXECUTION_FAILED";
        case HIPSOLVER_STATUS_INTERNAL_ERROR: return "HIPSOLVER_STATUS_INTERNAL_ERROR";
        case HIPSOLVER_STATUS_MAPPING_ERROR: return "HIPSOLVER_STATUS_MAPPING_ERROR";
        case HIPSOLVER_STATUS_NOT_SUPPORTED: return "HIPSOLVER_STATUS_NOT_SUPPORTED";
        case HIPSOLVER_STATUS_HANDLE_IS_NULLPTR: return "HIPSOLVER_STATUS_HANDLE_IS_NULLPTR";
        case HIPSOLVER_STATUS_INVALID_ENUM: return "HIPSOLVER_STATUS_INVALID_ENUM";
        default: return "unknown hipsolver error";
    }
}

void *dplasma_create_hip_handles(void *obj, void *_n)
{
    parsec_hip_exec_stream_t *stream = (parsec_hip_exec_stream_t *)obj;
    dplasma_hip_handles_t *new;
    hipblasHandle_t hipblas_handle;
    hipblasStatus_t hipblas_status;

    (void)_n;

    /* No need to call hipSetDevice, as this has been done by PaRSEC before calling the task body */
    hipblas_status = hipblasCreate(&hipblas_handle);
    if(HIPBLAS_STATUS_SUCCESS != hipblas_status) {
        if( HIPBLAS_STATUS_ALLOC_FAILED == hipblas_status) {
            parsec_show_help("help-dplasma.txt", "gpu_alloc_failed", 1, "HIPBLAS");
        }
        parsec_fatal("Unable to create HIPBLAS Handle: %s",
                     hipblasStatusToString(hipblas_status));
        return NULL;
    }
    hipblas_status = hipblasSetStream(hipblas_handle, stream->hip_stream);
    assert(HIPBLAS_STATUS_SUCCESS == hipblas_status);

    new = malloc(sizeof(dplasma_hip_handles_t));
    new->hipblas_handle = hipblas_handle;

    return new;
}

void dplasma_destroy_hip_handles(void *_h, void *_n)
{
    dplasma_hip_handles_t *handles = (dplasma_hip_handles_t*)_h;
    (void)_n;
    hipblasDestroy(handles->hipblas_handle);
    free(handles);
}
