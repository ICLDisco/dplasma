/*
 * Copyright (c) 2011-2022 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "parsec.h"
#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasmaaux.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

#include "ztrsmpl_ptgpanel.h"

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_ztrsmpl_ptgpanel_New - Generates the taskpool that solves U*x = b,
 * when U has been generated through LU factorization with partial pivoting
 * strategy implemented with ptgpanel strategy
 * See dplasma_zgetrf_ptgpanel_New().
 *
 * WARNING: The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A to be factorized.
 *          On entry, The factorized matrix through dplasma_zgetrf_ptgpanel_New()
 *          routine.  Elements on and above the diagonal are the elements of
 *          U. Elements below the diagonal are NOT the classic L, but the L
 *          factors obtaines by succesive pivoting.
 *
 * @param[in] IPIV
 *          Descriptor of the IPIV matrix. Should be distributed exactly as the
 *          A matrix. This matrix must be of size A.m - by - A.nt with IPIV.mb =
 *          A.mb and IPIV.nb = 1.
 *          On entry, contains the pivot indices of the successive row
 *          interchanged performed during the factorization.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, B is overwritten by the solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_ztrsmpl_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_ztrsmpl_ptgpanel
 * @sa dplasma_ztrsmpl_ptgpanel_Destruct
 * @sa dplasma_ctrsmpl_ptgpanel_New
 * @sa dplasma_dtrsmpl_ptgpanel_New
 * @sa dplasma_strsmpl_ptgpanel_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_ztrsmpl_ptgpanel_New( const parsec_tiled_matrix_t *A,
                            const parsec_tiled_matrix_t *IPIV,
                            parsec_tiled_matrix_t *B )
{
    parsec_ztrsmpl_ptgpanel_taskpool_t *parsec_ztrsmpl_ptgpanel = NULL;
    int nb = A->nb;
    int P = ((parsec_matrix_block_cyclic_t*)A)->grid.rows;

    parsec_ztrsmpl_ptgpanel = parsec_ztrsmpl_ptgpanel_new(A, IPIV, B, P);

    /* A */
    dplasma_add2arena_tile( &parsec_ztrsmpl_ptgpanel->arenas_datatypes[PARSEC_ztrsmpl_ptgpanel_DEFAULT_ADT_IDX],
                            A->mb*A->nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, A->mb );

    /* PERMUT */
    dplasma_add2arena_rectangle( &parsec_ztrsmpl_ptgpanel->arenas_datatypes[PARSEC_ztrsmpl_ptgpanel_PERMUT_ADT_IDX],
                                 2 * nb * sizeof(int),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_int_t, 2, nb, -1 );


    return (parsec_taskpool_t*)parsec_ztrsmpl_ptgpanel;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_ztrsmpl_ptgpanel_Destruct - Free the data structure associated to 
 *  a taskpool created with dplasma_ztrsmpl_ptgpanel_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_ztrsmpl_ptgpanel_New
 * @sa dplasma_ztrsmpl_ptgpanel
 *
 ******************************************************************************/
void
dplasma_ztrsmpl_ptgpanel_Destruct( parsec_taskpool_t *tp )
{
    parsec_ztrsmpl_ptgpanel_taskpool_t *parsec_ztrsmpl_ptgpanel = (parsec_ztrsmpl_ptgpanel_taskpool_t *)tp;

    dplasma_matrix_del2arena( &parsec_ztrsmpl_ptgpanel->arenas_datatypes[PARSEC_ztrsmpl_ptgpanel_DEFAULT_ADT_IDX] );
    dplasma_matrix_del2arena( &parsec_ztrsmpl_ptgpanel->arenas_datatypes[PARSEC_ztrsmpl_ptgpanel_PERMUT_ADT_IDX ] );

    parsec_taskpool_free(tp);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_ztrsmpl_ptgpanel - Solves U*x = b, when U has been generated
 * through LU factorization with partial pivoting implemented with the
 * ptgpanel strategy
 * See dplasma_zgetrf_ptgpanel_New().
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A to be factorized.
 *          On entry, The factorized matrix through dplasma_zgetrf_ptgpanel_New()
 *          routine.  Elements on and above the diagonal are the elements of
 *          U. Elements below the diagonal are NOT the classic L, but the L
 *          factors obtaines by succesive pivoting.
 *
 * @param[in] IPIV
 *          Descriptor of the IPIV matrix. Should be distributed exactly as the
 *          A matrix. This matrix must be of size A.m - by - A.nt with IPIV.mb =
 *          A.mb and IPIV.nb = 1.
 *          On entry, contains the pivot indices of the successive row
 *          interchanged performed during the factorization.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, B is overwritten by the solution matrix X.
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
 * @sa dplasma_ztrsmpl_ptgpanel
 * @sa dplasma_ztrsmpl_ptgpanel_Destruct
 * @sa dplasma_ctrsmpl_ptgpanel_New
 * @sa dplasma_dtrsmpl_ptgpanel_New
 * @sa dplasma_strsmpl_ptgpanel_New
 *
 ******************************************************************************/
int
dplasma_ztrsmpl_ptgpanel( parsec_context_t *parsec,
                        const parsec_tiled_matrix_t *A,
                        const parsec_tiled_matrix_t *IPIV,
                        parsec_tiled_matrix_t *B )
{
    parsec_taskpool_t *parsec_ztrsmpl_ptgpanel = NULL;

    parsec_ztrsmpl_ptgpanel = dplasma_ztrsmpl_ptgpanel_New(A, IPIV, B );

    if ( parsec_ztrsmpl_ptgpanel != NULL )
    {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_ztrsmpl_ptgpanel);
        dplasma_wait_until_completion(parsec);
        dplasma_ztrsmpl_ptgpanel_Destruct( parsec_ztrsmpl_ptgpanel );
    }

    return 0;
}

