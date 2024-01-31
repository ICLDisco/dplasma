/*
 * Copyright (c) 2023-     The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * $COPYRIGHT
 *
 */
#include "dplasma/config.h"
#include "parsec/utils/zone_malloc.h"
#include "parsec/utils/show_help.h"
#include <cublas_v2.h>
#include "dplasmaaux_cuda.h"
#include "potrf_cublas_utils.h"

/* 
 * Global info ID's for cublas handles and workspaces
 * Should be initialized in the tests
 * with the return of parsec_info_register
 * or parsec_info_lookup
 */
parsec_info_id_t CuHI = -1;
parsec_info_id_t WoSI = -1;

/* Unfortunately, CUBLAS does not provide a error to string function */
static char *dplasma_cublas_error_to_string(cublasStatus_t cublas_status)
{
    switch(cublas_status)
    {
        case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR";
        default: return "unknown CUBLAS error";
    }
}

/* Unfortunately, cuSolver does not provide a error to string function */
char *dplasma_cusolver_error_to_string(cusolverStatus_t cusolver_status)
{
    switch(cusolver_status) {
        case CUSOLVER_STATUS_SUCCESS: return "CUSOLVER_STATUS_SUCCESS";
        case CUSOLVER_STATUS_NOT_INITIALIZED: return "CUSOLVER_STATUS_NOT_INITIALIZED";
        case CUSOLVER_STATUS_ALLOC_FAILED: return "CUSOLVER_STATUS_ALLOC_FAILED";
        case CUSOLVER_STATUS_INVALID_VALUE: return "CUSOLVER_STATUS_INVALID_VALUE";
        case CUSOLVER_STATUS_ARCH_MISMATCH: return "CUSOLVER_STATUS_ARCH_MISMATCH";
        case CUSOLVER_STATUS_EXECUTION_FAILED: return "CUSOLVER_STATUS_EXECUTION_FAILED";
        case CUSOLVER_STATUS_INTERNAL_ERROR: return "CUSOLVER_STATUS_INTERNAL_ERROR";
        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED: return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
        default: return "unknown cusolver error";
    }
}

void *dplasma_create_cuda_handles(void *obj, void *_n)
{
    parsec_cuda_exec_stream_t *cuda_stream = (parsec_cuda_exec_stream_t *)obj;
    dplasma_cuda_handles_t *new;
    cublasHandle_t cublas_handle;
    cublasStatus_t cublas_status;

    (void)_n;

    /* No need to call cudaSetDevice, as this has been done by PaRSEC before calling the task body */
    cublas_status = cublasCreate(&cublas_handle);
    if(CUBLAS_STATUS_SUCCESS != cublas_status) {
        if( CUBLAS_STATUS_ALLOC_FAILED == cublas_status ) {
            parsec_show_help("help-dplasma.txt", "cu*_alloc_failed", 1, "CUBLAS");
        }
        parsec_fatal("Unable to create CUBLAS Handle: %s",
                     dplasma_cublas_error_to_string(cublas_status));
        return NULL;
    }
    cublas_status = cublasSetStream(cublas_handle, cuda_stream->cuda_stream);
    assert(CUBLAS_STATUS_SUCCESS == cublas_status);

    cusolverDnHandle_t cusolver_handle;
    cusolverStatus_t   cusolver_status;
    cusolver_status = cusolverDnCreate(&cusolver_handle);
    if(CUSOLVER_STATUS_SUCCESS != cusolver_status) {
        cublasDestroy(cublas_handle);
        if( CUSOLVER_STATUS_ALLOC_FAILED == cusolver_status ) {
            parsec_show_help("help-dplasma.txt", "cu*_alloc_failed", 1, "cusolver");
        }
        parsec_fatal("Unable to create a cuSolver handle: %s",
                     dplasma_cusolver_error_to_string(cusolver_status));
        return NULL;
    }
    cusolver_status = cusolverDnSetStream(cusolver_handle, cuda_stream->cuda_stream);
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);

    new = malloc(sizeof(dplasma_cuda_handles_t));
    new->cublas_handle = cublas_handle;
    new->cusolverDn_handle = cusolver_handle;

    return new;
}

void dplasma_destroy_cuda_handles(void *_h, void *_n)
{
    dplasma_cuda_handles_t *handles = (dplasma_cuda_handles_t*)_h;
    (void)_n;
    cublasDestroy_v2(handles->cublas_handle);
    cusolverDnDestroy(handles->cusolverDn_handle);
    free(handles);
}
