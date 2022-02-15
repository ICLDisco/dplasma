/*
 * Copyright (c) 2009-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "flops.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    dplasma_enum_t uplo = dplasmaUpper;
    int info = 0;
    int ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 0, 180, 180);

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZPOTRF, ((DagDouble_t)N));

    /* initializing matrix structure */
    LDA = dplasma_imax( LDA, N );
    LDB = dplasma_imax( LDB, N );
    KP = 1;
    KQ = 1;

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   N, N, P, nodes/P, uplo));
    int t;
    for(t = 0; t < nruns+1; t++) {
        /* matrix (re)generation */
        if(loud > 3) printf("+++ Generate matrices ... ");
        dplasma_zplghe( parsec, (double)(N), uplo,
                        (parsec_tiled_matrix_t *)&dcA, random_seed);
        if(loud > 3) printf("Done\n");

        parsec_devices_release_memory();

        if((iparam[IPARAM_HNB] != iparam[IPARAM_NB]) || (iparam[IPARAM_HMB] != iparam[IPARAM_MB]))
        {

            SYNC_TIME_START();
            parsec_taskpool_t* PARSEC_zpotrf = dplasma_zpotrf_New( uplo, (parsec_tiled_matrix_t*)&dcA, &info );
            /* Set the recursive size */
            dplasma_zpotrf_setrecursive( PARSEC_zpotrf, iparam[IPARAM_HMB] );
            parsec_context_add_taskpool(parsec, PARSEC_zpotrf);
            if( loud > 2 && t > 0) {
                SYNC_TIME_PRINT(rank, ( "zpotrf\tDAG created\n"));
            }

            PASTE_CODE_PROGRESS_KERNEL(parsec, zpotrf, t);
            dplasma_zpotrf_Destruct( PARSEC_zpotrf );

            parsec_taskpool_sync_ids(); /* recursive DAGs are not synchronous on ids */

        }
        else
        {
            PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec, zpotrf, 
                                      ( uplo, (parsec_tiled_matrix_t*)&dcA, &info),
                                      dplasma_zpotrf_Destruct( PARSEC_zpotrf ), t);
        }
        parsec_devices_reset_load(parsec);
    }
    PASTE_CODE_PERF_LOOP_DONE();

    if( 0 == rank && info != 0 ) {
        printf("-- Factorization is suspicious (info = %d) ! \n", info);
        ret |= 1;
    }
    if( !info && check ) {
        /* Check the factorization */
        PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
            parsec_matrix_sym_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                       rank, MB, NB, LDA, N, 0, 0,
                                       N, N, P, nodes/P, uplo));
        dplasma_zplghe( parsec, (double)(N), uplo,
                        (parsec_tiled_matrix_t *)&dcA0, random_seed);

        ret |= check_zpotrf( parsec, (rank == 0) ? loud : 0, uplo,
                             (parsec_tiled_matrix_t *)&dcA,
                             (parsec_tiled_matrix_t *)&dcA0);

        /* Check the solution */
        PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
            parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, NRHS, 0, 0,
                                   N, NRHS, P, nodes/P, KP, KQ, IP, JQ));
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, random_seed+1);

        PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
            parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, NRHS, 0, 0,
                                   N, NRHS, P, nodes/P, KP, KQ, IP, JQ));
        dplasma_zlacpy( parsec, dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcB, (parsec_tiled_matrix_t *)&dcX );

        dplasma_zpotrs(parsec, uplo,
                       (parsec_tiled_matrix_t *)&dcA,
                       (parsec_tiled_matrix_t *)&dcX );

        ret |= check_zaxmb( parsec, (rank == 0) ? loud : 0, uplo,
                            (parsec_tiled_matrix_t *)&dcA0,
                            (parsec_tiled_matrix_t *)&dcB,
                            (parsec_tiled_matrix_t *)&dcX);

        /* Cleanup */
        parsec_data_free(dcA0.mat); dcA0.mat = NULL;
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0 );
        parsec_data_free(dcB.mat); dcB.mat = NULL;
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB );
        parsec_data_free(dcX.mat); dcX.mat = NULL;
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcX );
    }

    parsec_data_free(dcA.mat); dcA.mat = NULL;
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);

    cleanup_parsec(parsec, iparam);
    return ret;
}
