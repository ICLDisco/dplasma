/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 *
 * @precisions normal z -> c d s
 *
 */

#include <lapacke.h>
#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasmaaux.h"
#include "parsec/data_dist/matrix/apply.h"


static int
dplasma_zlaset_operator( parsec_execution_stream_t *es,
                         const parsec_tiled_matrix_t *descA,
                         void *_A,
                         dplasma_enum_t uplo, int m, int n,
                         void *args )
{
    int tempmm, tempnn, ldam;
    dplasma_complex64_t *alpha = (dplasma_complex64_t*)args;
    dplasma_complex64_t *A = (dplasma_complex64_t*)_A;
    (void)es;

    tempmm = ((m)==((descA->mt)-1)) ? ((descA->m)-(m*(descA->mb))) : (descA->mb);
    tempnn = ((n)==((descA->nt)-1)) ? ((descA->n)-(n*(descA->nb))) : (descA->nb);
    ldam = BLKLDD( descA, m );

    if (m == n) {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR, dplasma_lapack_const( uplo ), tempmm, tempnn,
            alpha[0], alpha[1], A, ldam);
    } else {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR, 'A', tempmm, tempnn,
            alpha[0], alpha[0], A, ldam);
    }
    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zlaset_New - Generates the taskpool that set the elements of the matrix
 * A on the diagonal to beta and the off-diagonals eklements to alpha.
 *
 * See parsec_apply_New() for further information.
 *
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which part of matrix A is set:
 *          = dplasmaUpperLower: All matrix is referenced.
 *          = dplasmaUpper:      Only upper part is referenced.
 *          = dplasmaLower:      Only lower part is referenced.
 *
 * @param[in] alpha
 *         The constant to which the off-diagonal elements are to be set.
 *
 * @param[in] beta
 *         The constant to which the diagonal elements are to be set.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A. Any tiled matrix
 *          descriptor can be used.
 *          On exit, A has been set accordingly.
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_zlaset_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_zlaset
 * @sa dplasma_zlaset_Destruct
 * @sa dplasma_claset_New
 * @sa dplasma_dlaset_New
 * @sa dplasma_slaset_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_zlaset_New( dplasma_enum_t uplo,
                    dplasma_complex64_t alpha,
                    dplasma_complex64_t beta,
                    parsec_tiled_matrix_t *A )
{
    dplasma_complex64_t *params = (dplasma_complex64_t*)malloc(2 * sizeof(dplasma_complex64_t));

    params[0] = alpha;
    params[1] = beta;

    return parsec_apply_New( uplo, A, dplasma_zlaset_operator, params );
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zlaset_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_zlaset_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_zlaset_New
 * @sa dplasma_zlaset
 *
 ******************************************************************************/
void
dplasma_zlaset_Destruct( parsec_taskpool_t *tp )
{
    if( ((parsec_apply_taskpool_t *)tp)->_g_op_args ) {
        free( ((parsec_apply_taskpool_t *)tp)->_g_op_args );
    }
    parsec_apply_Destruct(tp);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zlaset - Set the elements of the matrix
 * A on the diagonal to beta and the off-diagonals eklements to alpha.
 *
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] uplo
 *          Specifies which part of matrix A is set:
 *          = dplasmaUpperLower: All matrix is referenced.
 *          = dplasmaUpper:      Only upper part is refrenced.
 *          = dplasmaLower:      Only lower part is referenced.
 *
 * @param[in] alpha
 *         The constant to which the off-diagonal elements are to be set.
 *
 * @param[in] beta
 *         The constant to which the diagonal elements are to be set.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A. Any tiled matrix
 *          descriptor can be used.
 *          On exit, A has been set accordingly.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa dplasma_zlaset_New
 * @sa dplasma_zlaset_Destruct
 * @sa dplasma_claset
 * @sa dplasma_dlaset
 * @sa dplasma_slaset
 *
 ******************************************************************************/
int
dplasma_zlaset( parsec_context_t *parsec,
                dplasma_enum_t uplo,
                dplasma_complex64_t alpha,
                dplasma_complex64_t beta,
                parsec_tiled_matrix_t *A )
{
    parsec_taskpool_t *parsec_zlaset = NULL;

    /* Check input arguments */
    if ((uplo != dplasmaLower) &&
        (uplo != dplasmaUpper) &&
        (uplo != dplasmaUpperLower))
    {
        dplasma_error("dplasma_zlaset", "illegal value of type");
        return -1;
    }

    parsec_zlaset = dplasma_zlaset_New(uplo, alpha, beta, A);

    if ( parsec_zlaset != NULL ) {
        parsec_context_add_taskpool(parsec, (parsec_taskpool_t*)parsec_zlaset);
        dplasma_wait_until_completion(parsec);
        dplasma_zlaset_Destruct( parsec_zlaset );
    }
    return 0;
}
