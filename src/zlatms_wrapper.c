/*
 * Copyright (c) 2011-2022 The University of Tennessee and The University
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
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

/**
 *******************************************************************************
 *
 *  Generic case
 *
 *******************************************************************************
 */
static int
dplasma_zlatms_operator( parsec_execution_stream_t *es,
                         const parsec_tiled_matrix_t *descA,
                         void *_A,
                         dplasma_enum_t uplo, int m, int n,
                         void *args )
{
    int tempmm, tempnn, ldam, i;
    double            *cond = (double*)args;
    dplasma_complex64_t *A    = (dplasma_complex64_t*)_A;
    (void)es;

    tempmm = ((m)==((descA->mt)-1)) ? ((descA->m)-(m*(descA->mb))) : (descA->mb);
    tempnn = ((n)==((descA->nt)-1)) ? ((descA->n)-(n*(descA->nb))) : (descA->nb);
    ldam = BLKLDD( descA, m );

    if (m == n) {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR, dplasma_lapack_const( uplo ), tempmm, tempnn,
            0., 0., A, ldam);

        /* D(i) = 1 - (i-1)/(N-1)*(1 - 1/COND) */
        if ( m == 0 ) {
            A[0] = 1.;
            i = 1;
        }
        else {
            i = 0;
        }
        if ( descA->n > 1 ) {
            double tmp = 1. / (*cond);
            double alp = ( 1. - tmp ) / ((double)( descA->n - 1 ));
            int minmn = dplasma_imin( tempmm, tempnn );
            for(; i < minmn; i++){
                A[i+i*ldam] = (dplasma_complex64_t)( (double)(descA->n-(descA->nb*n+i+1)) * alp + tmp );
            }
        }
    } else {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR, 'A', tempmm, tempnn,
            0., 0., A, ldam);
    }
    return 0;
}



/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zpltmg - Generates a special test matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] mtxtype
 *           - dplasmaGeneral:   Generate a general matrix
 *           - dplasmaSymmetric: Generate a symmetric matrix
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A to generate. Any tiled matrix
 *          descriptor can be used.
 *          On exit, the matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa dplasma_cpltmg
 * @sa dplasma_dpltmg
 * @sa dplasma_spltmg
 *
 ******************************************************************************/
int
dplasma_zlatms( parsec_context_t *parsec,
                dplasma_enum_t mtxtype, double cond,
                parsec_tiled_matrix_t *A,
                unsigned long long int seed)
{
    parsec_matrix_block_cyclic_t Q, T;
    int nodes, rank, mb, nb, m, n, mt, nt, P, IP, JQ;

    /* Init the diagonal of A */
    {
        parsec_taskpool_t *tp;
        double *condptr = malloc(sizeof( double ));
        *condptr = cond;
        tp = parsec_apply_New( dplasmaUpperLower, A, dplasma_zlatms_operator, condptr );
        if ( tp != NULL ) {
            parsec_context_add_taskpool(parsec, tp);
            parsec_context_start( parsec );
            parsec_context_wait( parsec );
            parsec_apply_Destruct(tp);
            free(condptr);
        }
        else {
            free(condptr);
            return -1;
        }
    }

    nodes = A->super.nodes;
    rank  = A->super.myrank;
    mb    = A->mb;
    nb    = A->nb;
    m     = A->m;
    n     = A->n;
    mt    = A->mt;
    nt    = A->nt;

    if ( A->dtype & parsec_matrix_block_cyclic_type ) {
        P = ((parsec_matrix_block_cyclic_t*)A)->grid.rows;
        IP = ((parsec_matrix_block_cyclic_t*)A)->grid.ip;
        JQ = ((parsec_matrix_block_cyclic_t*)A)->grid.jq;
    }
    else if ( A->dtype & parsec_matrix_sym_block_cyclic_type ) {
        P = ((parsec_matrix_sym_block_cyclic_t*)A)->grid.rows;
        IP = ((parsec_matrix_sym_block_cyclic_t*)A)->grid.ip;
        JQ = ((parsec_matrix_sym_block_cyclic_t*)A)->grid.jq;
    }
    else {
        P = 1;
        IP = 0;
        JQ = 0;
    }

    /* Init the random matrix R */
    parsec_matrix_block_cyclic_init( &Q, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, mb, nb, m, n, 0, 0, m, n, P, nodes/P, 1, 1, IP, JQ );
    Q.mat = parsec_data_allocate((size_t)Q.super.nb_local_tiles *
                                (size_t)Q.super.bsiz *
                                (size_t)parsec_datadist_getsizeoftype(Q.super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t*)&Q, "Q");

    /* Init the T matrix */
    parsec_matrix_block_cyclic_init( &T, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, 32, nb, mt*32, n, 0, 0, mt*32, n, P, nodes/P, 1, 1, IP, JQ );
    T.mat = parsec_data_allocate((size_t)T.super.nb_local_tiles *
                                (size_t)T.super.bsiz *
                                (size_t)parsec_datadist_getsizeoftype(T.super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t*)&T, "T");

    if ( mtxtype == dplasmaGeneral ) {
        if ( m >= n ) {
            parsec_tiled_matrix_t *subA = parsec_tiled_matrix_submatrix( A, 0, 0, n, n );
            parsec_tiled_matrix_t *subQ = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&Q,
                                                                0, 0, n, n );
            parsec_tiled_matrix_t *subT = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&T,
                                                                0, 0, nt*32, n );


            /* Multiply on the right by an unitary matrix */
            dplasma_zplrnt( parsec, 0, subQ, seed + 1 );
            dplasma_zgeqrf( parsec, subQ, subT );
            dplasma_zunmqr( parsec, dplasmaRight, dplasmaNoTrans,
                            subQ, subT, subA );

            /* Multiply on the left by an unitary matrix */
            dplasma_zplrnt( parsec, 0,
                            (parsec_tiled_matrix_t *)&Q, seed );
            dplasma_zgeqrf( parsec,
                            (parsec_tiled_matrix_t*)&Q,
                            (parsec_tiled_matrix_t*)&T );
            dplasma_zunmqr( parsec, dplasmaLeft, dplasmaNoTrans,
                            (parsec_tiled_matrix_t*)&Q,
                            (parsec_tiled_matrix_t*)&T, A );

            free(subA); free(subQ); free(subT);
        }
        else {
            parsec_tiled_matrix_t *subA = parsec_tiled_matrix_submatrix( A, 0, 0, m, m );
            parsec_tiled_matrix_t *subQ = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&Q,
                                                                0, 0, m, m );
            parsec_tiled_matrix_t *subT = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&T,
                                                                0, 0, mt*32, m );


            /* Multiply on the left by an unitary matrix */
            dplasma_zplrnt( parsec, 0, subQ, seed );
            dplasma_zgeqrf( parsec, subQ, subT );
            dplasma_zunmqr( parsec, dplasmaLeft, dplasmaNoTrans,
                            subQ, subT, subA );

            /* Multiply on the right by an unitary matrix */
            dplasma_zplrnt( parsec, 0,
                            (parsec_tiled_matrix_t *)&Q, seed );
            dplasma_zgeqrf( parsec,
                            (parsec_tiled_matrix_t*)&Q,
                            (parsec_tiled_matrix_t*)&T );
            dplasma_zunmqr( parsec, dplasmaRight, dplasmaNoTrans,
                            (parsec_tiled_matrix_t*)&Q,
                            (parsec_tiled_matrix_t*)&T, A );
            free(subA); free(subQ); free(subT);
        }
    }
    else {
        assert( mtxtype == dplasmaHermitian );

        /* Init the unitary matrix */
        dplasma_zplrnt( parsec, 0,
                        (parsec_tiled_matrix_t *)&Q, seed );
        dplasma_zgeqrf( parsec,
                        (parsec_tiled_matrix_t*)&Q,
                        (parsec_tiled_matrix_t*)&T );

        dplasma_zunmqr( parsec, dplasmaLeft, dplasmaNoTrans,
                        (parsec_tiled_matrix_t*)&Q,
                        (parsec_tiled_matrix_t*)&T, A );
        dplasma_zunmqr( parsec, dplasmaRight, dplasmaConjTrans,
                        (parsec_tiled_matrix_t*)&Q,
                        (parsec_tiled_matrix_t*)&T, A );
    }

    free(Q.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&Q );
    free(T.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&T );

    return 0;
}
