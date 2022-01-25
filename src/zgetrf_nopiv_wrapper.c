/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
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

#include "zgetrf_nopiv.h"

#define MAX_SHAPES 1

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zgetrf_nopiv_New - Generates the taskpool that computes the LU
 * factorization of a M-by-N matrix A: A = L * U by with no pivoting
 * strategy. The matrix has to be diaagonal dominant to use this
 * routine. Otherwise, the numerical stability of the result is not guaranted.
 *
 * Other variants of LU decomposition with pivoting stragies are available in
 * the library with the following function:
 *     - dplasma_zgetrf_New() that performs LU decomposition with partial pivoting.
 *       This is limited to matrices with column-cyclic distribution.
 *     - dplasma_zgetrf_incpiv_New() that performs tile incremental pivoting
 *       algorithm.
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
 *          A = L*U; the unit diagonal elements of L are not stored.
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
 *          destroy with dplasma_zgetrf_nopiv_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_zgetrf_nopiv
 * @sa dplasma_zgetrf_nopiv_Destruct
 * @sa dplasma_cgetrf_nopiv_New
 * @sa dplasma_dgetrf_nopiv_New
 * @sa dplasma_sgetrf_nopiv_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_zgetrf_nopiv_New( parsec_tiled_matrix_t *A,
                          int *INFO )
{

    dplasma_data_collection_t * ddc_A = dplasma_wrap_data_collection(A);

    parsec_zgetrf_nopiv_taskpool_t *parsec_getrf_nopiv;
    parsec_getrf_nopiv = parsec_zgetrf_nopiv_new( ddc_A, INFO );

    int shape = 0;
    dplasma_setup_adtt_all_loc( ddc_A,
                                parsec_datatype_double_complex_t,
                                PARSEC_MATRIX_FULL/*uplo*/, 1/*diag:for PARSEC_MATRIX_UPPER or PARSEC_MATRIX_LOWER types*/,
                                &shape);

    assert(shape == MAX_SHAPES);
    return (parsec_taskpool_t*)parsec_getrf_nopiv;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zgetrf_nopiv_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_zgetrf_nopiv_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_zgetrf_nopiv_New
 * @sa dplasma_zgetrf_nopiv
 *
 ******************************************************************************/
void
dplasma_zgetrf_nopiv_Destruct( parsec_taskpool_t *tp )
{
    parsec_zgetrf_nopiv_taskpool_t *parsec_zgetrf_nopiv = (parsec_zgetrf_nopiv_taskpool_t *)tp;
    dplasma_clean_adtt_all_loc(parsec_zgetrf_nopiv->_g_ddescA, MAX_SHAPES);
    dplasma_data_collection_t * ddc_A = parsec_zgetrf_nopiv->_g_ddescA;

    parsec_taskpool_free(tp);

    /* free the dplasma_data_collection_t */
    dplasma_unwrap_data_collection(ddc_A);

}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zgetrf_nopiv - Computes the LU factorization of a M-by-N matrix A: A
 * = L * U by with no pivoting strategy. The matrix has to be diaagonal dominant
 * to use this routine. Otherwise, the numerical stability of the result is not
 * guaranted.
 *
 * Other variants of LU decomposition with pivoting stragies are available in
 * the library with the following function:
 *     - dplasma_zgetrf_New() that performs LU decomposition with partial pivoting.
 *       This is limited to matrices with column-cyclic distribution.
 *     - dplasma_zgetrf_incpiv_New() that performs tile incremental pivoting
 *       algorithm.
 *     - dplasma_zgetrf_qrf_New() that performs an hybrid LU-QR decomposition.
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
 *          A = L*U; the unit diagonal elements of L are not stored.
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
 * @sa dplasma_zgetrf_nopiv
 * @sa dplasma_zgetrf_nopiv_Destruct
 * @sa dplasma_cgetrf_nopiv_New
 * @sa dplasma_dgetrf_nopiv_New
 * @sa dplasma_sgetrf_nopiv_New
 *
 ******************************************************************************/
int
dplasma_zgetrf_nopiv( parsec_context_t *parsec,
                      parsec_tiled_matrix_t *A )
{
    parsec_taskpool_t *parsec_zgetrf_nopiv = NULL;

    int info = 0;
    parsec_zgetrf_nopiv = dplasma_zgetrf_nopiv_New(A, &info);

    if ( parsec_zgetrf_nopiv != NULL ) {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zgetrf_nopiv);
        dplasma_wait_until_completion(parsec);
        dplasma_zgetrf_nopiv_Destruct( parsec_zgetrf_nopiv );
        return info;
    }
    else
        return -101;
}
