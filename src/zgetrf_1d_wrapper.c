/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
 *
 * @precisions normal z -> s d c
 *
 */

#include "dplasma.h"
#include "parsec/vpmap.h"
#include "dplasmajdf.h"
#include "dplasma/types.h"
#include "dplasma/types_lapack.h"

#include "zgetrf_1d.h"

#define MAX_SHAPES 2

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zgetrf_1d_New - Generates the taskpool that computes the LU factorization
 * of a M-by-N matrix A: A = P * L * U by partial pivoting algorithm.
 *
 * This algorithm exploits the multi-threaded recursive kernels of the PLASMA
 * library and by consequence require a column-cyclic data distribution if used
 * in distributed memory.
 * This is not an optimal solution for distributed memory system, and should be
 * used only if no other possibiliies is available. Absolute priority scheduler
 * is known to improve the performance of this algorithm and should be prefered.
 *
 * Other variants of LU decomposition are available in the library wioth the
 * following function:
 *     - dplasma_zgetrf_incpiv_New() that performs tile incremental pivoting
 *       algorithm.
 *     - dplasma_zgetrf_nopiv_New() that performs LU decomposition with no pivoting
 *       if the matrix is known as beeing diagonal dominant.
 *     - dplasma_zgetrf_qrf_New() that performs an hybrid LU-QR decomposition.
 *
 * WARNING: The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A to be factorized.
 *          On entry, describes the M-by-N matrix A.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 * @param[out] IPIV
 *          Descriptor of the IPIV matrix. Should be of size 1-by-min(M,N).
 *          On exit, contains the pivot indices; for 1 <= i <= min(M,N), row i
 *          of the matrix was interchanged with row IPIV(i).
 *
 * @param[out] INFO
 *          On algorithm completion: equal to 0 on success, i if the ith
 *          diagonal value is equal to 0. That implies incoherent result.
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_zgetrf_1d_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_zgetrf_1d
 * @sa dplasma_zgetrf_1d_Destruct
 * @sa dplasma_cgetrf_1d_New
 * @sa dplasma_dgetrf_1d_New
 * @sa dplasma_sgetrf_1d_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_zgetrf_1d_New( parsec_tiled_matrix_t *A,
                       parsec_tiled_matrix_t *IPIV,
                       int *INFO )
{
    parsec_zgetrf_1d_taskpool_t *parsec_getrf_1d;
    int nbthreads = dplasma_imax( 1, parsec_vpmap_get_vp_threads(0) - 1 );
    dplasma_data_collection_t * ddc_A = dplasma_wrap_data_collection((parsec_tiled_matrix_t*)A);
    dplasma_data_collection_t * ddc_IPIV = dplasma_wrap_data_collection((parsec_tiled_matrix_t*)IPIV);

    if ( (IPIV->mt != 1) || (dplasma_imin(A->nt, A->mt) > IPIV->nt)) {
        dplasma_error("dplasma_zgetrf_1d_New", "IPIV doesn't have the correct number of tiles (1-by-min(A->mt,A->nt)");
        return NULL;
    }

    parsec_getrf_1d = parsec_zgetrf_1d_new( ddc_A,
                                            ddc_IPIV,
                                            INFO );

    if ( A->storage == PARSEC_MATRIX_TILE ) {
        parsec_getrf_1d->_g_getrfdata = CORE_zgetrf_rectil_init(nbthreads);
    } else {
        parsec_getrf_1d->_g_getrfdata = CORE_zgetrf_reclap_init(nbthreads);
    }
    parsec_getrf_1d->_g_nbmaxthrd = nbthreads;

    int shape = 0;
    dplasma_setup_adtt_all_loc( ddc_A,
                                parsec_datatype_double_complex_t,
                                PARSEC_MATRIX_FULL/*uplo*/, 1/*diag:for PARSEC_MATRIX_UPPER or PARSEC_MATRIX_LOWER types*/,
                                &shape);

    dplasma_setup_adtt_all_loc( ddc_IPIV,
                                parsec_datatype_int_t,
                                PARSEC_MATRIX_FULL/*uplo*/, 1/*diag:for PARSEC_MATRIX_UPPER or PARSEC_MATRIX_LOWER types*/,
                                &shape);
    assert(shape == MAX_SHAPES);

    return (parsec_taskpool_t*)parsec_getrf_1d;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zgetrf_1d_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_zgetrf_1d_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the tasdkpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_zgetrf_1d_New
 * @sa dplasma_zgetrf_1d
 *
 ******************************************************************************/
void
dplasma_zgetrf_1d_Destruct( parsec_taskpool_t *tp )
{
    parsec_zgetrf_1d_taskpool_t *parsec_zgetrf_1d = (parsec_zgetrf_1d_taskpool_t *)tp;
    dplasma_clean_adtt_all_loc(parsec_zgetrf_1d->_g_ddescA, MAX_SHAPES);
    dplasma_clean_adtt_all_loc(parsec_zgetrf_1d->_g_ddescIPIV, MAX_SHAPES);

    dplasma_data_collection_t * ddc_A = parsec_zgetrf_1d->_g_ddescA;
    dplasma_data_collection_t * ddc_IPIV = parsec_zgetrf_1d->_g_ddescIPIV;

    if ( parsec_zgetrf_1d->_g_getrfdata != NULL )
        free( parsec_zgetrf_1d->_g_getrfdata );

    parsec_taskpool_free(tp);

    /* free the dplasma_data_collection_t */
    dplasma_unwrap_data_collection(ddc_A);
    dplasma_unwrap_data_collection(ddc_IPIV);

}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zgetrf_1d - Computes the LU factorization of a M-by-N matrix A: A = P *
 * L * U by partial pivoting algorithm.
 *
 * This algorithm exploits the multi-threaded recursive kernels of the PLASMA
 * library and by consequence require a column-cyclic data distribution if used
 * in distributed memory.
 * This is not an optimal solution for distributed memory system, and should be
 * used only if no other possibiliies is available. Absolute priority scheduler
 * is known to improve the performance of this algorithm and should be prefered.
 *
 * Other variants of LU decomposition are available in the library wioth the
 * following function:
 *     - dplasma_zgetrf_incpiv() that performs tile incremental pivoting
 *       algorithm.
 *     - dplasma_zgetrf_nopiv() that performs LU decomposition with no pivoting
 *       if the matrix is known as beeing diagonal dominant.
 *     - dplasma_zgetrf_qrf() that performs an hybrid LU-QR decomposition.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A to be factorized.
 *          On entry, describes the M-by-N matrix A.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 * @param[out] IPIV
 *          Descriptor of the IPIV matrix. Should be of size 1-by-min(M,N).
 *          On exit, contains the pivot indices; for 1 <= i <= min(M,N), row i
 *          of the matrix was interchanged with row IPIV(i).
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *          \retval i if ith value is singular. Result is incoherent.
 *
 *******************************************************************************
 *
 * @sa dplasma_zgetrf_1d
 * @sa dplasma_zgetrf_1d_Destruct
 * @sa dplasma_cgetrf_1d_New
 * @sa dplasma_dgetrf_1d_New
 * @sa dplasma_sgetrf_1d_New
 *
 ******************************************************************************/
int
dplasma_zgetrf_1d( parsec_context_t *parsec,
                parsec_tiled_matrix_t *A,
                parsec_tiled_matrix_t *IPIV )
{
    parsec_taskpool_t *parsec_zgetrf_1d = NULL;

    int info = 0;

    if ( (IPIV->mt != 1) || (dplasma_imin(A->nt, A->mt) > IPIV->nt)) {
        dplasma_error("dplasma_zgetrf_1d", "IPIV doesn't have the correct number of tiles (1-by-min(A->mt,A->nt)");
        return -3;
    }

    parsec_zgetrf_1d = dplasma_zgetrf_1d_New(A, IPIV, &info);

    if ( parsec_zgetrf_1d != NULL ) {
        parsec_context_add_taskpool( parsec, parsec_zgetrf_1d );
        dplasma_wait_until_completion(parsec);
        dplasma_zgetrf_1d_Destruct( parsec_zgetrf_1d );
        return info;
    }
    else
        return -101;
}
