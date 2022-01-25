/*
 * Copyright (c) 2011-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
 *
 */

#include "dplasma.h"
#include "parsec/vpmap.h"
#include <math.h>
#include <alloca.h>
#include <string.h>
#include "dplasmaaux.h"
#include "parsec/utils/show_help.h"

#if defined(PARSEC_HAVE_MPI)
/*
 * dplasma falls back to MPI_COMM_WORLD by default.
 * This is sub-optimal, as this does not provide an opportunity
 * to insulate dplasma communications from potential others, but
 * this allows to maintain the behavior that dplasma does not
 * need initialization / finalization.
 *
 * The dplasma API provides two functions to provide and free
 * a dplasma-specific communicator if needed (these should be called
 * before any other dplasma call and after any dplasma call, respectively)
 */

static MPI_Comm dplasma_comm = MPI_COMM_WORLD;
void *dplasma_pcomm = &dplasma_comm;

int dplasma_aux_dup_comm(void *_psrc)
{
    MPI_Comm *src = (MPI_Comm*)_psrc;
    return MPI_Comm_dup(*src, &dplasma_comm);
}

int dplasma_aux_free_comm(void)
{
    return MPI_Comm_free(&dplasma_comm);
}
#else
void *dplasma_pcomm = NULL;
int dplasma_aux_dup_comm(void *comm)
{
    return -1;
}

int dplasma_aux_free_comm(void)
{
    return -1;
}
#endif


int
dplasma_aux_get_priority_limit( char* function, const parsec_tiled_matrix_t* dc )
{
    char *v;
    char *keyword;

    if( NULL == function || NULL == dc )
        return 0;

    keyword = alloca( strlen(function)+2 );

    switch( dc->mtype ) {
    case PARSEC_MATRIX_FLOAT:
        sprintf(keyword, "S%s", function);
        break;
    case PARSEC_MATRIX_DOUBLE:
        sprintf(keyword, "D%s", function);
        break;
    case PARSEC_MATRIX_COMPLEX_FLOAT:
        sprintf(keyword, "C%s", function);
        break;
    case PARSEC_MATRIX_COMPLEX_DOUBLE:
        sprintf(keyword, "Z%s", function);
        break;
    default:
        return 0;
    }

    if( (v = getenv(keyword)) != NULL ) {
        return atoi(v);
    }
    return 0;
}

int
dplasma_aux_getGEMMLookahead( parsec_tiled_matrix_t *A )
{
    /**
     * Assume that the number of threads per node is constant, and compute the
     * look ahead based on the global information to get the same one on all
     * nodes.
     */
    int nbunits = vpmap_get_nb_total_threads() * A->super.nodes;
    double alpha =  3. * (double)nbunits / ( A->mt * A->nt );

    if ( A->super.nodes == 1 ) {
        /* No look ahaead */
        return dplasma_imax( A->mt, A->nt );
    }
    else {
        /* Look ahaed of at least 2, and that provides 3 tiles per computational units */
        return dplasma_imax( ceil( alpha ), 2 );
    }
}

#if defined(DPLASMA_HAVE_CUDA)
#include <cublas_v2.h>
#include <cusolverDn.h>
#include "potrf_cublas_utils.h"
#include "parsec/utils/zone_malloc.h"

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
static char *dplasma_cusolver_error_to_string(cusolverStatus_t cusolver_status)
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

#endif
