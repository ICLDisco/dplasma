/*
 * Copyright (c) 2010-2023 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
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

#include "zpotrf_U.h"
#include "zpotrf_L.h"
#include "cores/dplasma_plasmatypes.h"

#define MAX_SHAPES 1

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zpotrf_setrecursive - Set the recursive size parameter to enable
 *  recursive DAGs.
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to modify.
 *          On exit, the modified taskpool.
 *
 * @param[in] hmb
 *          The tile size to use for the smaller recursive call.
 *          hmb must be > 0, otherwise nothing is changed.
 *
 *******************************************************************************
 *
 * @sa dplasma_zpotrf_New
 * @sa dplasma_zpotrf
 *
 ******************************************************************************/
void
dplasma_zpotrf_setrecursive( parsec_taskpool_t *tp, int hmb )
{
    parsec_zpotrf_L_taskpool_t *parsec_zpotrf = (parsec_zpotrf_L_taskpool_t*)tp;
    if (hmb > 0) {
        parsec_zpotrf->_g_smallnb = hmb;
    }
}

#if defined(DPLASMA_HAVE_CUDA)
static void *zpotrf_create_cuda_workspace(void *obj, void *user)
{
    parsec_device_module_t *mod = (parsec_device_module_t *)obj;
    zone_malloc_t *memory = ((parsec_device_gpu_module_t*)mod)->memory;
    cusolverDnHandle_t cusolverDnHandle;
    cusolverStatus_t status;
    parsec_zpotrf_U_taskpool_t *tp = (parsec_zpotrf_U_taskpool_t*)user;
    dplasma_potrf_workspace_t *wp = NULL;
    int workspace_size;
    int mb = tp->_g_descA->mb;
    int nb = tp->_g_descA->nb;
    size_t elt_size = sizeof(cuDoubleComplex);
    cublasFillMode_t cublas_uplo;
    dplasma_enum_t uplo = tp->_g_uplo;

    if( PlasmaLower == uplo )
        cublas_uplo = CUBLAS_FILL_MODE_LOWER;
    if( PlasmaUpper == uplo )
        cublas_uplo = CUBLAS_FILL_MODE_UPPER;

    status = cusolverDnCreate(&cusolverDnHandle);
    assert(CUSOLVER_STATUS_SUCCESS == status);
    (void)status;

    status = cusolverDnDpotrf_bufferSize(cusolverDnHandle, cublas_uplo, nb, NULL, mb, &workspace_size);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    cusolverDnDestroy(cusolverDnHandle);

    wp = (dplasma_potrf_workspace_t*)malloc(sizeof(dplasma_potrf_workspace_t));
    wp->tmpmem = zone_malloc(memory, workspace_size * elt_size + sizeof(int));
    assert(NULL != wp->tmpmem);
    wp->lwork = workspace_size;
    wp->memory = memory;

    return wp;
}

static void zpotrf_destroy_cuda_workspace(void *_ws, void *_n)
{
    dplasma_potrf_workspace_t *ws = (dplasma_potrf_workspace_t*)_ws;
    zone_free((zone_malloc_t*)ws->memory, ws->tmpmem);
    free(ws);
    (void)_n;
}
#endif

#if defined(DPLASMA_HAVE_HIP)
static void *zpotrf_create_hip_workspace(void *obj, void *user)
{
    parsec_device_module_t *mod = (parsec_device_module_t *)obj;
    zone_malloc_t *memory = ((parsec_device_gpu_module_t*)mod)->memory;
    dplasma_potrf_workspace_t *wp = NULL;
    (void)user;

    wp = (dplasma_potrf_workspace_t*)malloc(sizeof(dplasma_potrf_workspace_t));
    wp->tmpmem = zone_malloc(memory, sizeof(int));
    assert(NULL != wp->tmpmem);
    wp->lwork = 0;
    wp->memory = memory;

    return wp;
}

static void zpotrf_destroy_hip_workspace(void *_ws, void *_n)
{
    dplasma_potrf_workspace_t *ws = (dplasma_potrf_workspace_t*)_ws;
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
 * dplasma_zpotrf_New - Generates the taskpool that Computes the Cholesky
 * factorization of a symmetric positive definite (or Hermitian positive
 * definite in the complex case) matrix A, with or without recursive calls.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = dplasmaLower}^{U^H\times U, if uplo = dplasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 * WARNING: The computations are not done by this call.
 *
 * If you want to enable the recursive DAGs, don't forget to set the recursive
 * tile size and to synchonize the taskpool ids after the computations since those
 * are for now local. You can follow the code of dplasma_zpotrf_rec() as an
 * example to do this.
 *
 * Hierarchical DAG Scheduling for Hybrid Distributed Systems; Wu, Wei and
 * Bouteiller, Aurelien and Bosilca, George and Faverge, Mathieu and Dongarra,
 * Jack. 29th IEEE International Parallel & Distributed Processing Symposium,
 * May 2015. (https://hal.inria.fr/hal-0107835)
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = dplasmaUpper: Upper triangle of A is referenced;
 *          = dplasmaLower: Lower triangle of A is referenced.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten with the factorized
 *          matrix.
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
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_zpotrf_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_zpotrf
 * @sa dplasma_zpotrf_Destruct
 * @sa dplasma_cpotrf_New
 * @sa dplasma_dpotrf_New
 * @sa dplasma_spotrf_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_zpotrf_New( dplasma_enum_t uplo,
                    parsec_tiled_matrix_t *A,
                    int *info )
{
    parsec_zpotrf_L_taskpool_t *parsec_zpotrf = NULL;
#if defined(DPLASMA_HAVE_CUDA) || defined(DPLASMA_HAVE_HIP)
    char workspace_info_name[64];
    static int uid = 0;
#endif
    parsec_taskpool_t *tp = NULL;
    dplasma_data_collection_t * ddc_A = dplasma_wrap_data_collection(A);

    /* Check input arguments */
    if ((uplo != dplasmaUpper) && (uplo != dplasmaLower)) {
        dplasma_error("dplasma_zpotrf_New", "illegal value of uplo");
        return NULL /*-1*/;
    }

    *info = 0;

    if ( uplo == dplasmaUpper ) {
        tp = (parsec_taskpool_t*)parsec_zpotrf_U_new( uplo, ddc_A, info);
    } else {
        tp = (parsec_taskpool_t*)parsec_zpotrf_L_new( uplo, ddc_A, info);
    }
    parsec_zpotrf =  (parsec_zpotrf_L_taskpool_t*)tp;

    parsec_zpotrf->_g_PRI_CHANGE = dplasma_aux_get_priority_limit( "POTRF", A );
    if(0 == parsec_zpotrf->_g_PRI_CHANGE)
      parsec_zpotrf->_g_PRI_CHANGE = A->nt;
#if defined(DPLASMA_HAVE_CUDA)
    /* It doesn't cost anything to define these infos if we have CUDA but
     * don't have GPUs on the current machine, so we do it non-conditionally */
    parsec_zpotrf->_g_cuda_handles_infokey = parsec_info_lookup(&parsec_per_stream_infos, "DPLASMA::CUDA::HANDLES", NULL);
    snprintf(workspace_info_name, 64, "DPLASMA::ZPOTRF(%d)::WS", uid++);
    parsec_zpotrf->_g_cuda_workspaces_infokey = parsec_info_register(&parsec_per_device_infos, workspace_info_name,
                                                           zpotrf_destroy_cuda_workspace, NULL,
                                                           zpotrf_create_cuda_workspace, parsec_zpotrf,
                                                           NULL);
#else
    parsec_zpotrf->_g_cuda_handles_infokey = PARSEC_INFO_ID_UNDEFINED;
    parsec_zpotrf->_g_cuda_workspaces_infokey = PARSEC_INFO_ID_UNDEFINED;
    (void)uid; (void)workspace_info_name;
#endif
#if defined(PARSEC_HAVE_HIP)
    /* It doesn't cost anything to define these infos if we have HIP but
     * don't have GPUs on the current machine, so we do it non-conditionally */
    parsec_zpotrf->_g_hip_handles_infokey = parsec_info_lookup(&parsec_per_stream_infos, "DPLASMA::HIP::HANDLES", NULL);
    snprintf(workspace_info_name, 64, "DPLASMA::ZPOTRF(%d)::WS", uid++);
    parsec_zpotrf->_g_hip_workspaces_infokey = parsec_info_register(&parsec_per_device_infos, workspace_info_name,
                                                           zpotrf_destroy_hip_workspace, NULL,
                                                           zpotrf_create_hip_workspace, parsec_zpotrf,
                                                           NULL);
#else
    parsec_zpotrf->_g_hip_handles_infokey = PARSEC_INFO_ID_UNDEFINED;
    parsec_zpotrf->_g_hip_workspaces_infokey = PARSEC_INFO_ID_UNDEFINED;
#endif
    int shape = 0;
    dplasma_setup_adtt_all_loc( ddc_A,
                                parsec_datatype_double_complex_t,
                                PARSEC_MATRIX_FULL/*uplo*/, 1/*diag:for PARSEC_MATRIX_UPPER or PARSEC_MATRIX_LOWER types*/,
                                &shape);

    assert(shape == MAX_SHAPES);
    return tp;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zpotrf_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_zpotrf_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_zpotrf_New
 * @sa dplasma_zpotrf
 *
 ******************************************************************************/
void
dplasma_zpotrf_Destruct( parsec_taskpool_t *tp )
{
    parsec_zpotrf_L_taskpool_t *parsec_zpotrf = (parsec_zpotrf_L_taskpool_t *)tp;
    dplasma_clean_adtt_all_loc(parsec_zpotrf->_g_ddescA, MAX_SHAPES);
    dplasma_data_collection_t * ddc_A = parsec_zpotrf->_g_ddescA;

#if defined(DPLASMA_HAVE_CUDA)
    parsec_info_unregister(&parsec_per_device_infos, parsec_zpotrf->_g_cuda_workspaces_infokey, NULL);
#endif
#if defined(DPLASMA_HAVE_HIP)
    parsec_info_unregister(&parsec_per_device_infos, parsec_zpotrf->_g_hip_workspaces_infokey, NULL);
#endif

    parsec_taskpool_free(tp);
    /* free the dplasma_data_collection_t */
    dplasma_unwrap_data_collection(ddc_A);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zpotrf - Computes the Cholesky factorization of a symmetric positive
 * definite (or Hermitian positive definite in the complex case) matrix A.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = dplasmaLower}^{U^H\times U, if uplo = dplasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
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
 *          On exit, the uplo part of A is overwritten with the factorized
 *          matrix.
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
 *
 *******************************************************************************
 *
 * @sa dplasma_zpotrf_New
 * @sa dplasma_zpotrf_Destruct
 * @sa dplasma_cpotrf
 * @sa dplasma_dpotrf
 * @sa dplasma_spotrf
 *
 ******************************************************************************/
int
dplasma_zpotrf( parsec_context_t *parsec,
                dplasma_enum_t uplo,
                parsec_tiled_matrix_t *A )
{
    parsec_taskpool_t *parsec_zpotrf = NULL;
    int info = 0, ginfo = 0 ;

    parsec_zpotrf = dplasma_zpotrf_New( uplo, A, &info );

    if ( parsec_zpotrf != NULL )
    {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zpotrf);
        dplasma_wait_until_completion(parsec);
        dplasma_zpotrf_Destruct( parsec_zpotrf );
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
 * dplasma_zpotrf_rec - Computes the Cholesky factorization of a symmetric
 * positive definite (or Hermitian positive definite in the complex case) matrix
 * An, using the recursive DAGs feature if hmb is smaller than A.mb or A.nb.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = dplasmaLower}^{U^H\times U, if uplo = dplasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
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
 *          On exit, the uplo part of A is overwritten with the factorized
 *          matrix.
 *
 * @param[in] hmb
 *          The tile size to use for the smaller recursive call.
 *          If hmb <= 0 or hmb > A.mb, the classic algorithm without recursive
 *          calls is applied.
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
 *
 *******************************************************************************
 *
 * @sa dplasma_zpotrf_New
 * @sa dplasma_zpotrf_Destruct
 * @sa dplasma_cpotrf
 * @sa dplasma_dpotrf
 * @sa dplasma_spotrf
 *
 ******************************************************************************/
int
dplasma_zpotrf_rec( parsec_context_t *parsec,
                    dplasma_enum_t uplo,
                    parsec_tiled_matrix_t *A, int hmb )
{
    parsec_taskpool_t *parsec_zpotrf = NULL;
    int info = 0, ginfo = 0 ;

    parsec_zpotrf = dplasma_zpotrf_New( uplo, A, &info );
    if ( parsec_zpotrf != NULL )
    {
        dplasma_zpotrf_setrecursive( (parsec_taskpool_t*)parsec_zpotrf, hmb );
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zpotrf);
        dplasma_wait_until_completion(parsec);
        dplasma_zpotrf_Destruct( parsec_zpotrf );
        parsec_taskpool_sync_ids(); /* recursive DAGs are not synchronous on ids */
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
