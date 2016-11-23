/*
 * Copyright (c) 2011-2013 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 *
 * @precisions normal z -> c d s
 *
 */

#include <lapacke.h>
#include "dplasma.h"
#include "dplasma/lib/dplasmatypes.h"
#include "data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "data_dist/matrix/two_dim_rectangle_cyclic.h"

#include "map.h"

/**
 *******************************************************************************
 *
 *  Generic case
 *
 *******************************************************************************
 */
static int
dplasma_zlatms_operator( dague_execution_unit_t *eu,
                         const tiled_matrix_desc_t *descA,
                         void *_A,
                         PLASMA_enum uplo, int m, int n,
                         void *args )
{
    int tempmm, tempnn, ldam, i;
    double            *cond = (double*)args;
    dague_complex64_t *A    = (dague_complex64_t*)_A;
    (void)eu;

    tempmm = ((m)==((descA->mt)-1)) ? ((descA->m)-(m*(descA->mb))) : (descA->mb);
    tempnn = ((n)==((descA->nt)-1)) ? ((descA->n)-(n*(descA->nb))) : (descA->nb);
    ldam = BLKLDD( descA, m );

    if (m == n) {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR, lapack_const( uplo ), tempmm, tempnn,
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
                A[i+i*ldam] = (dague_complex64_t)( (double)(descA->n-(descA->nb*n+i+1)) * alp + tmp );
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
 * @param[in,out] dague
 *          The dague context of the application that will run the operation.
 *
 * @param[in] mtxtype
 *           - PlasmaGeneral:   Generate a general matrix
 *           - PlasmaSymmetric: Generate a symmetric matrix
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
dplasma_zlatms( dague_context_t *dague,
                PLASMA_enum mtxtype, double cond,
                tiled_matrix_desc_t *A,
                unsigned long long int seed)
{
    two_dim_block_cyclic_t Q, T;
    int nodes, rank, mb, nb, m, n, mt, nt, P;

    /* Init the diagonal of A */
    {
        dague_handle_t *handle;
        double *condptr = malloc(sizeof( double ));
        *condptr = cond;
        handle = dplasma_map_New( PlasmaUpperLower, A, dplasma_zlatms_operator, condptr );
        if ( handle != NULL ) {
            dague_enqueue(dague, handle);
            dague_context_wait( dague );
            dplasma_map_Destruct( handle );
        }
        else {
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

    if ( A->dtype & two_dim_block_cyclic_type ) {
        P = ((two_dim_block_cyclic_t*)A)->grid.rows;
    }
    else if ( A->dtype & sym_two_dim_block_cyclic_type ) {
        P = ((sym_two_dim_block_cyclic_t*)A)->grid.rows;
    }
    else {
        P = 1;
    }

    /* Init the random matrix R */
    two_dim_block_cyclic_init( &Q, matrix_ComplexDouble, matrix_Tile,
                               nodes, rank, mb, nb, m, n, 0, 0, m, n, 1, 1, P );
    Q.mat = dague_data_allocate((size_t)Q.super.nb_local_tiles *
                                (size_t)Q.super.bsiz *
                                (size_t)dague_datadist_getsizeoftype(Q.super.mtype));
    dague_ddesc_set_key((dague_ddesc_t*)&Q, "Q");

    /* Init the T matrix */
    two_dim_block_cyclic_init( &T, matrix_ComplexDouble, matrix_Tile,
                               nodes, rank, 32, nb, mt*32, n, 0, 0, mt*32, n, 1, 1, P );
    T.mat = dague_data_allocate((size_t)T.super.nb_local_tiles *
                                (size_t)T.super.bsiz *
                                (size_t)dague_datadist_getsizeoftype(T.super.mtype));
    dague_ddesc_set_key((dague_ddesc_t*)&T, "T");

    if ( mtxtype == PlasmaGeneral ) {
        if ( m >= n ) {
            tiled_matrix_desc_t *subA = tiled_matrix_submatrix( A, 0, 0, n, n );
            tiled_matrix_desc_t *subQ = tiled_matrix_submatrix( (tiled_matrix_desc_t *)&Q,
                                                                0, 0, n, n );
            tiled_matrix_desc_t *subT = tiled_matrix_submatrix( (tiled_matrix_desc_t *)&T,
                                                                0, 0, nt*32, n );


            /* Multiply on the right by an unitary matrix */
            dplasma_zplrnt( dague, 0, subQ, seed + 1 );
            dplasma_zgeqrf( dague, subQ, subT );
            dplasma_zunmqr( dague, PlasmaRight, PlasmaNoTrans,
                            subQ, subT, subA );

            /* Multiply on the left by an unitary matrix */
            dplasma_zplrnt( dague, 0,
                            (tiled_matrix_desc_t *)&Q, seed );
            dplasma_zgeqrf( dague,
                            (tiled_matrix_desc_t*)&Q,
                            (tiled_matrix_desc_t*)&T );
            dplasma_zunmqr( dague, PlasmaLeft, PlasmaNoTrans,
                            (tiled_matrix_desc_t*)&Q,
                            (tiled_matrix_desc_t*)&T, A );
        }
        else {
            tiled_matrix_desc_t *subA = tiled_matrix_submatrix( A, 0, 0, m, m );
            tiled_matrix_desc_t *subQ = tiled_matrix_submatrix( (tiled_matrix_desc_t *)&Q,
                                                                0, 0, m, m );
            tiled_matrix_desc_t *subT = tiled_matrix_submatrix( (tiled_matrix_desc_t *)&T,
                                                                0, 0, mt*32, m );


            /* Multiply on the left by an unitary matrix */
            dplasma_zplrnt( dague, 0, subQ, seed );
            dplasma_zgeqrf( dague, subQ, subT );
            dplasma_zunmqr( dague, PlasmaLeft, PlasmaNoTrans,
                            subQ, subT, subA );

            /* Multiply on the right by an unitary matrix */
            dplasma_zplrnt( dague, 0,
                            (tiled_matrix_desc_t *)&Q, seed );
            dplasma_zgeqrf( dague,
                            (tiled_matrix_desc_t*)&Q,
                            (tiled_matrix_desc_t*)&T );
            dplasma_zunmqr( dague, PlasmaRight, PlasmaNoTrans,
                            (tiled_matrix_desc_t*)&Q,
                            (tiled_matrix_desc_t*)&T, A );
        }
    }
    else {
        assert( mtxtype == PlasmaHermitian );

        /* Init the unitary matrix */
        dplasma_zplrnt( dague, 0,
                        (tiled_matrix_desc_t *)&Q, seed );
        dplasma_zgeqrf( dague,
                        (tiled_matrix_desc_t*)&Q,
                        (tiled_matrix_desc_t*)&T );

        dplasma_zunmqr( dague, PlasmaLeft, PlasmaNoTrans,
                        (tiled_matrix_desc_t*)&Q,
                        (tiled_matrix_desc_t*)&T, A );
        dplasma_zunmqr( dague, PlasmaRight, PlasmaConjTrans,
                        (tiled_matrix_desc_t*)&Q,
                        (tiled_matrix_desc_t*)&T, A );
    }
    return 0;
}