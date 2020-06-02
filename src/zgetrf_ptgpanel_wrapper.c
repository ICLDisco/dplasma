/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "parsec.h"
#include "dplasma/types.h"
#include "dplasmaaux.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

#include "zgetrf_ptgpanel.h"

#define IB  32

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zgetrf_ptgpanel_New - Generates the taskpool that computes the LU
 * factorization of a M-by-N matrix A: A = P * L * U by partial pivoting
 * algorithm.
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
 *          destroy with dplasma_zgetrf_ptgpanel_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_zgetrf_ptgpanel
 * @sa dplasma_zgetrf_ptg_panel_Destruct
 * @sa dplasma_cgetrf_ptgpanel_New
 * @sa dplasma_dgetrf_ptgpanel_New
 * @sa dplasma_sgetrf_ptgpanel_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_zgetrf_ptgpanel_New( parsec_tiled_matrix_dc_t *A,
                           parsec_tiled_matrix_dc_t *IPIV,
                           int *info )
{
    parsec_zgetrf_ptgpanel_taskpool_t *parsec_zgetrf_ptgpanel = NULL;
    int nb = A->nb;
    int P = ((two_dim_block_cyclic_t*)A)->grid.rows;
    int Q = ((two_dim_block_cyclic_t*)A)->grid.cols;

    /* The code has to be fixed for N >> M */
    assert( A->m >= A->n );

    *info = 0;
    parsec_zgetrf_ptgpanel = parsec_zgetrf_ptgpanel_new( A, IPIV, IB,
                                                     P, Q, info);

    /* A */
    dplasma_add2arena_tile( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_DEFAULT_ARENA],
                            A->mb*A->nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, A->mb );

    /* SWAP */
    dplasma_add2arena_rectangle( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_SWAP_ARENA],
                                 (2*nb+1)*sizeof(dplasma_complex64_t),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_double_complex_t, 2*nb+1, 1, -1 );

    /* MAXL */
    dplasma_add2arena_rectangle( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_MAXL_ARENA],
                                 (nb+1)*sizeof(dplasma_complex64_t),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_double_complex_t, 1, nb+1, -1 );

    /* UMES */
    dplasma_add2arena_rectangle( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_UMES_ARENA],
                                 IB*nb*sizeof(dplasma_complex64_t),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_double_complex_t, IB, nb, -1 );

    /* PIVOT */
    dplasma_add2arena_rectangle( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_PIVOT_ARENA],
                                 A->mb*sizeof(int),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_int_t, IPIV->mb, IPIV->nb, -1 );

    /* PERMUT */
    dplasma_add2arena_rectangle( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_PERMUT_ARENA],
                                 2 * nb * sizeof(int),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_int_t, 2, nb, -1 );

    return (parsec_taskpool_t*)parsec_zgetrf_ptgpanel;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zgetrf_ptgpanel_Destruct - Free the data structure associated to a
 *  taskpool created with dplasma_zgetrf_ptgpanel_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the tasdkpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_zgetrf_ptgpanel_New
 * @sa dplasma_zgetrf_ptgpanel
 *
 ******************************************************************************/
void
dplasma_zgetrf_ptgpanel_Destruct( parsec_taskpool_t *tp )
{
    parsec_zgetrf_ptgpanel_taskpool_t *parsec_zgetrf_ptgpanel = (parsec_zgetrf_ptgpanel_taskpool_t *)tp;

    parsec_matrix_del2arena( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_DEFAULT_ARENA] );
    parsec_matrix_del2arena( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_SWAP_ARENA   ] );
    parsec_matrix_del2arena( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_MAXL_ARENA   ] );
    parsec_matrix_del2arena( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_UMES_ARENA   ] );
    parsec_matrix_del2arena( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_PIVOT_ARENA  ] );
    parsec_matrix_del2arena( parsec_zgetrf_ptgpanel->arenas[PARSEC_zgetrf_ptgpanel_PERMUT_ARENA ] );

    parsec_taskpool_free(tp);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zgetrf_ptgpanel - Computes the LU factorization of a M-by-N 
 * matrix A: A = P * L * U by partial pivoting algorithm.
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
 * @sa dplasma_zgetrf_ptgpanel
 * @sa dplasma_zgetrf_ptgpanel_Destruct
 * @sa dplasma_cgetrf_ptgpanel_New
 * @sa dplasma_dgetrf_ptgpanel_New
 * @sa dplasma_sgetrf_ptgpanel_New
 *
 ******************************************************************************/
int
dplasma_zgetrf_ptgpanel( parsec_context_t *parsec,
                       parsec_tiled_matrix_dc_t *A,
                       parsec_tiled_matrix_dc_t *IPIV )
{
    int info = 0, ginfo = 0 ;
    parsec_taskpool_t *parsec_zgetrf_ptgpanel = NULL;

    parsec_zgetrf_ptgpanel = dplasma_zgetrf_ptgpanel_New(A, IPIV, &info);

    if ( parsec_zgetrf_ptgpanel != NULL )
    {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zgetrf_ptgpanel);
        dplasma_wait_until_completion(parsec);
        dplasma_zgetrf_ptgpanel_Destruct( parsec_zgetrf_ptgpanel );
    }

#if defined(HAVE_MPI)
    MPI_Allreduce( &info, &ginfo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    ginfo = info;
#endif
    return ginfo;
}
