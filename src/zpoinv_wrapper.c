/*
 * Copyright (c) 2010-2022 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2026      NVIDIA Corporation.  All rights reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasma/types_lapack.h"
#include "dplasmaaux.h"
#include "potrf_gpu_workspaces.h"
#include "parsec/utils/zone_malloc.h"

#include "zpoinv_U.h"
#include "zpoinv_L.h"

#if defined(DPLASMA_HAVE_CUDA)
#include <cusolverDn.h>

static void *zpoinv_create_cuda_workspace(void *obj, void *user)
{
    parsec_device_module_t *mod = (parsec_device_module_t *)obj;
    zone_malloc_t *memory = ((parsec_device_gpu_module_t*)mod)->memory;
    cusolverDnHandle_t cusolverDnHandle;
    cusolverStatus_t status;
    parsec_zpoinv_U_taskpool_t *tp = (parsec_zpoinv_U_taskpool_t*)user;
    dplasma_potrf_gpu_workspaces_t *wp = NULL;
    int workspace_size;
    int mb = tp->_g_descA->mb;
    int nb = tp->_g_descA->nb;
    size_t elt_size = sizeof(cuDoubleComplex);
    dplasma_enum_t uplo = tp->_g_uplo;

    status = cusolverDnCreate(&cusolverDnHandle);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    (void)status;

    status = cusolverDnDpotrf_bufferSize(cusolverDnHandle, dplasma_cublas_fill(uplo), nb, NULL, mb, &workspace_size);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    cusolverDnDestroy(cusolverDnHandle);

    wp = (dplasma_potrf_gpu_workspaces_t*)malloc(sizeof(dplasma_potrf_gpu_workspaces_t));
    wp->tmpmem = zone_malloc(memory, workspace_size * elt_size + sizeof(int));
    assert(NULL != wp->tmpmem);
    wp->lwork = workspace_size;
    wp->memory = memory;

    return wp;
}

static void zpoinv_destroy_cuda_workspace(void *_ws, void *_n)
{
    dplasma_potrf_gpu_workspaces_t *ws = (dplasma_potrf_gpu_workspaces_t*)_ws;
    zone_free((zone_malloc_t*)ws->memory, ws->tmpmem);
    free(ws);
    (void)_n;
}
#endif

#if defined(DPLASMA_HAVE_HIP)
static void *zpoinv_create_hip_workspace(void *obj, void *user)
{
    parsec_device_module_t *mod = (parsec_device_module_t *)obj;
    zone_malloc_t *memory = ((parsec_device_gpu_module_t*)mod)->memory;
    dplasma_potrf_gpu_workspaces_t *wp = NULL;
    (void)user;

    wp = (dplasma_potrf_gpu_workspaces_t*)malloc(sizeof(dplasma_potrf_gpu_workspaces_t));
    wp->tmpmem = zone_malloc(memory, sizeof(int));
    assert(NULL != wp->tmpmem);
    wp->lwork = 0;
    wp->memory = memory;

    return wp;
}

static void zpoinv_destroy_hip_workspace(void *_ws, void *_n)
{
    dplasma_potrf_gpu_workspaces_t *ws = (dplasma_potrf_gpu_workspaces_t*)_ws;
    zone_free((zone_malloc_t*)ws->memory, ws->tmpmem);
    free(ws);
    (void)_n;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zpoinv_New - Generates the taskpool that computes the inverse of an
 * hermitian matrix through Cholesky factorization and inversion.
 *
 * WARNING: The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = dplasmaUpper: Upper triangle of A is referenced;
 *          = dplasmaLower: Lower triangle of A is referenced.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten by the inverse of A,
 *          A^(-1)
 *
 * @param[out] info
 *          Address where to store the output information of the factorization,
 *          this is not synchronized between the nodes, and might not be set
 *          when function exists.
 *          On DAG completion:
 *              - info = 0 on all nodes if successful.
 *              - info > 0 if the leading minor of order i of A is not positive
 *                definite, so the factorization could not be completed, and the
 *                solution has not been computed. Info will be equal to i on the
 *                node that owns the diagonal element (i,i), and 0 on all other
 *                nodes.
 *         TODO: support this correctly as it is made in Cholesky
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_zpoinv_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_zpoinv
 * @sa dplasma_zpoinv_Destruct
 * @sa dplasma_cpoinv_New
 * @sa dplasma_dpoinv_New
 * @sa dplasma_spoinv_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_zpoinv_New( dplasma_enum_t uplo,
                    parsec_tiled_matrix_t *A,
                    int *info )
{
    parsec_zpoinv_L_taskpool_t *parsec_zpoinv = NULL;
#if defined(DPLASMA_HAVE_CUDA) || defined(DPLASMA_HAVE_HIP)
    char workspace_info_name[64];
    static int uid = 0;
#endif
    parsec_taskpool_t *tp = NULL;

    /* Check input arguments */
    if ((uplo != dplasmaUpper) && (uplo != dplasmaLower)) {
        dplasma_error("dplasma_zpoinv_New", "illegal value of uplo");
        return NULL /*-1*/;
    }

    *info = 0;
    if ( uplo == dplasmaUpper ) {
        tp = (parsec_taskpool_t*)parsec_zpoinv_U_new( uplo, A, info );

        /* Upper part of A with diagonal part */
        /* dplasma_add2arena_upper( &((parsec_zpoinv_U_taskpool_t*)parsec_poinv)->arenas_datatypes[PARSEC_zpoinv_U_UPPER_TILE_ADT_IDX], */
        /*                          A->mb*A->nb*sizeof(dplasma_complex64_t), */
        /*                          PARSEC_ARENA_ALIGNMENT_SSE, */
        /*                          parsec_datatype_double_complex_t, A->mb, 1 ); */
    } else {
        tp = (parsec_taskpool_t*)parsec_zpoinv_L_new( uplo, A, info );

        /* Lower part of A with diagonal part */
        /* dplasma_add2arena_lower( &((parsec_zpoinv_L_taskpool_t*)parsec_poinv)->arenas_datatypes[PARSEC_zpoinv_L_LOWER_TILE_ADT_IDX], */
        /*                          A->mb*A->nb*sizeof(dplasma_complex64_t), */
        /*                          PARSEC_ARENA_ALIGNMENT_SSE, */
        /*                          parsec_datatype_double_complex_t, A->mb, 1 ); */
    }

    parsec_zpoinv = (parsec_zpoinv_L_taskpool_t*)tp;

#if defined(DPLASMA_HAVE_CUDA)
    /* It doesn't cost anything to define these infos if we have CUDA but
     * don't have GPUs on the current machine, so we do it non-conditionally. */
    parsec_zpoinv->_g_cuda_handles_infokey = parsec_info_lookup(&parsec_per_stream_infos, "DPLASMA::CUDA::HANDLES", NULL);
    snprintf(workspace_info_name, 64, "DPLASMA::ZPOINV(%d)::CUDA::WS", uid++);
    parsec_zpoinv->_g_cuda_workspaces_infokey = parsec_info_register(&parsec_per_device_infos, workspace_info_name,
                                                           zpoinv_destroy_cuda_workspace, NULL,
                                                           zpoinv_create_cuda_workspace, parsec_zpoinv,
                                                           NULL);
#else
    parsec_zpoinv->_g_cuda_handles_infokey = PARSEC_INFO_ID_UNDEFINED;
    parsec_zpoinv->_g_cuda_workspaces_infokey = PARSEC_INFO_ID_UNDEFINED;
#endif
#if defined(DPLASMA_HAVE_HIP)
    /* It doesn't cost anything to define these infos if we have HIP but
     * don't have GPUs on the current machine, so we do it non-conditionally. */
    parsec_zpoinv->_g_hip_handles_infokey = parsec_info_lookup(&parsec_per_stream_infos, "DPLASMA::HIP::HANDLES", NULL);
    snprintf(workspace_info_name, 64, "DPLASMA::ZPOINV(%d)::HIP::WS", uid++);
    parsec_zpoinv->_g_hip_workspaces_infokey = parsec_info_register(&parsec_per_device_infos, workspace_info_name,
                                                          zpoinv_destroy_hip_workspace, NULL,
                                                          zpoinv_create_hip_workspace, parsec_zpoinv,
                                                          NULL);
#else
    parsec_zpoinv->_g_hip_handles_infokey = PARSEC_INFO_ID_UNDEFINED;
    parsec_zpoinv->_g_hip_workspaces_infokey = PARSEC_INFO_ID_UNDEFINED;
#endif

    dplasma_add2arena_tile( &parsec_zpoinv->arenas_datatypes[PARSEC_zpoinv_L_DEFAULT_ADT_IDX],
                            A->mb*A->nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, A->mb );


    return tp;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zpoinv_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_zpoinv_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_zpoinv_New
 * @sa dplasma_zpoinv
 *
 ******************************************************************************/
void
dplasma_zpoinv_Destruct( parsec_taskpool_t *tp )
{
    parsec_zpoinv_L_taskpool_t *parsec_zpoinv = (parsec_zpoinv_L_taskpool_t *)tp;

    dplasma_matrix_del2arena( &parsec_zpoinv->arenas_datatypes[PARSEC_zpoinv_L_DEFAULT_ADT_IDX   ] );
    /* dplasma_matrix_del2arena( parsec_zpoinv->arenas_datatypes[PARSEC_zpoinv_L_LOWER_TILE_ADT_IDX] ); */
#if defined(DPLASMA_HAVE_CUDA)
    parsec_info_unregister(&parsec_per_device_infos, parsec_zpoinv->_g_cuda_workspaces_infokey, NULL);
#endif
#if defined(DPLASMA_HAVE_HIP)
    parsec_info_unregister(&parsec_per_device_infos, parsec_zpoinv->_g_hip_workspaces_infokey, NULL);
#endif
    parsec_taskpool_free(tp);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zpoinv - Computes the matrix inverse of an hermitian matrix through
 * Cholesky factorization and inversion.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] uplo
 *          = dplasmaUpper: Upper triangle of A is referenced;
 *          = dplasmaLower: Lower triangle of A is referenced.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten with inverse of A.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *          \retval > 0 if the leading minor of order i of A is not positive
 *          definite, so the factorization could not be completed, and the
 *          solution has not been computed. Info will be equal to i on the node
 *          that owns the diagonal element (i,i), and 0 on all other nodes.
 *         TODO: support this correctly as it is made in Cholesky
 *
 *******************************************************************************
 *
 * @sa dplasma_zpoinv_New
 * @sa dplasma_zpoinv_Destruct
 * @sa dplasma_cpoinv
 * @sa dplasma_dpoinv
 * @sa dplasma_spoinv
 *
 ******************************************************************************/
int
dplasma_zpoinv( parsec_context_t *parsec,
                dplasma_enum_t uplo,
                parsec_tiled_matrix_t *A )
{
    parsec_taskpool_t *parsec_zpoinv = NULL;
    int info = 0, ginfo = 0 ;

    parsec_zpoinv = dplasma_zpoinv_New( uplo, A, &info );

    if ( parsec_zpoinv != NULL )
    {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zpoinv);
        dplasma_wait_until_completion(parsec);
        dplasma_zpoinv_Destruct( parsec_zpoinv );
    }

    /* This covers both cases when we have not compiled with MPI, or we don't need to do the reduce */
    ginfo = info;
#if defined(PARSEC_HAVE_MPI)
    /* If we don't need to reduce, don't do it, this way we don't require MPI to be initialized */
    if( A->super.nodes > 1 )
        MPI_Allreduce( &info, &ginfo, 1, MPI_INT, MPI_MAX, *(MPI_Comm*)dplasma_pcomm);
#endif
    return ginfo;
}


/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zpoinv_sync - Computes the matrix inverse of an hermitian matrix
 * through Cholesky factorization and inversion as in dplasma_zpoinv. The
 * difference is in the fact that it calls successively three different DAGs
 * with intermediate synchronizations.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] uplo
 *          = dplasmaUpper: Upper triangle of A is referenced;
 *          = dplasmaLower: Lower triangle of A is referenced.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten with inverse of A.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *          \retval > 0 if the leading minor of order i of A is not positive
 *          definite, so the factorization could not be completed, and the
 *          solution has not been computed. Info will be equal to i on the node
 *          that owns the diagonal element (i,i), and 0 on all other nodes.
 *         TODO: support this correctly as it is made in Cholesky
 *
 *******************************************************************************
 *
 * @sa dplasma_zpoinv_New
 * @sa dplasma_zpoinv_Destruct
 * @sa dplasma_cpoinv
 * @sa dplasma_dpoinv
 * @sa dplasma_spoinv
 *
 ******************************************************************************/
int
dplasma_zpoinv_sync( parsec_context_t *parsec,
                     dplasma_enum_t uplo,
                     parsec_tiled_matrix_t* A )
{
    int info = 0;
    /* Check input arguments */
    if (uplo != dplasmaUpper && uplo != dplasmaLower) {
        dplasma_error("dplasma_zpoinv_sync", "illegal value of uplo");
        return -1;
    }

    info = dplasma_zpotrf( parsec, uplo, A );
    info = dplasma_ztrtri( parsec, uplo, dplasmaNonUnit, A );
    dplasma_zlauum( parsec, uplo, A );

    return info;
}
