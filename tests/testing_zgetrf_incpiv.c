/*
 * Copyright (c) 2009-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX );

static int check_inverse( parsec_context_t *parsec, int loud,
                          parsec_tiled_matrix_t *dcA,
                          parsec_tiled_matrix_t *dcInvA,
                          parsec_tiled_matrix_t *dcI );

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int info = 0;
    int ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 40, 200, 200);
    iparam[IPARAM_KP] = 4;
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'m';

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);

    dplasma_warmup(parsec);

    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZGETRF, ((DagDouble_t)M,(DagDouble_t)N));

    LDA = max(M, LDA);

    if ( M != N && check ) {
        fprintf(stderr, "Check is impossible if M != N\n");
        check = 0;
    }

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcL, 1,
        parsec_matrix_block_cyclic, (&dcL, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcIPIV, 1,
        parsec_matrix_block_cyclic, (&dcIPIV, PARSEC_MATRIX_INTEGER, PARSEC_MATRIX_TILE,
                               rank, MB, 1, M, NT, 0, 0,
                               M, NT, P, nodes/P, KP, KQ, IP, JQ));

    int t;
    for(t = 0; t < nruns; t++) {
        /* matrix (re)generation */
        if(loud > 2) printf("+++ Generate matrices ... ");
        dplasma_zpltmg( parsec, matrix_init, (parsec_tiled_matrix_t *)&dcA, random_seed );
        if(loud > 2) printf("Done\n");

        /* Create PaRSEC */
        if(loud > 2) printf("+++ Computing getrf_incpiv ... ");
        PASTE_CODE_ENQUEUE_KERNEL(parsec,
            zgetrf_incpiv, ((parsec_tiled_matrix_t*)&dcA,
                            (parsec_tiled_matrix_t*)&dcL,
                            (parsec_tiled_matrix_t*)&dcIPIV,
                            &info));
        /* lets rock! */
        PASTE_CODE_PROGRESS_KERNEL(parsec, zgetrf_incpiv);
        dplasma_zgetrf_incpiv_Destruct( PARSEC_zgetrf_incpiv );
        if(loud > 2) printf("Done.\n");
    }
    PASTE_CODE_PERF_LOOP_DONE();

    if ( info != 0 ) {
        if( rank == 0 && loud ) printf("-- Factorization is suspicious (info = %d) ! \n", info );
        ret |= 1;
    }
    else if ( check ) {
        /* regenerate A0 from seed */
        PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
            parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   M, N, P, nodes/P, KP, KQ, IP, JQ));
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA0, random_seed);

        /* First check: Ax=B */
        PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
            parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, NRHS, 0, 0,
                                   M, NRHS, P, nodes/P, KP, KQ, IP, JQ));
        PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
            parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, NRHS, 0, 0,
                                   M, NRHS, P, nodes/P, KP, KQ, IP, JQ));
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, random_seed+1);
        dplasma_zlacpy( parsec, dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcB,
                        (parsec_tiled_matrix_t *)&dcX );

        dplasma_zgetrs_incpiv(parsec, dplasmaNoTrans,
                              (parsec_tiled_matrix_t *)&dcA,
                              (parsec_tiled_matrix_t *)&dcL,
                              (parsec_tiled_matrix_t *)&dcIPIV,
                              (parsec_tiled_matrix_t *)&dcX );
        ret |= check_solution( parsec, (rank == 0) ? loud : 0,
                               (parsec_tiled_matrix_t *)&dcA0,
                               (parsec_tiled_matrix_t *)&dcB,
                               (parsec_tiled_matrix_t *)&dcX);

        parsec_data_free(dcB.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
        parsec_data_free(dcX.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcX);

        /* Second check: Inverse check */
        if ( check_inv ) {
            PASTE_CODE_ALLOCATE_MATRIX(dcI, check_inv,
                parsec_matrix_block_cyclic, (&dcI, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                       rank, MB, NB, LDA, N, 0, 0,
                                       M, N, P, nodes/P, KP, KQ, IP, JQ));
            PASTE_CODE_ALLOCATE_MATRIX(dcInvA, check_inv,
                parsec_matrix_block_cyclic, (&dcInvA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                       rank, MB, NB, LDA, N, 0, 0,
                                       M, N, P, nodes/P, KP, KQ, IP, JQ));
            dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&dcI);
            dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&dcInvA);

            dplasma_zgetrs_incpiv(parsec, dplasmaNoTrans,
                           (parsec_tiled_matrix_t *)&dcA,
                           (parsec_tiled_matrix_t *)&dcL,
                           (parsec_tiled_matrix_t *)&dcIPIV,
                           (parsec_tiled_matrix_t *)&dcInvA );
            ret |= check_inverse(parsec, (rank == 0) ? loud : 0,
                                 (parsec_tiled_matrix_t *)&dcA0,
                                 (parsec_tiled_matrix_t *)&dcInvA,
                                 (parsec_tiled_matrix_t *)&dcI);

            parsec_data_free(dcInvA.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcInvA);
            parsec_data_free(dcI.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcI);
        }

        parsec_data_free(dcA0.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0);
    }

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    parsec_data_free(dcL.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcL);
    parsec_data_free(dcIPIV.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcIPIV);

    cleanup_parsec(parsec, iparam);

    return ret;
}

static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX )
{
    int info_solution;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double Bnorm = 0.0;
    double Xnorm, result;
    int m = dcB->m;
    double eps = LAPACKE_dlamch_work('e');

    Anorm = dplasma_zlange(parsec, dplasmaInfNorm, dcA);
    Bnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);
    Xnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcX);

    /* Compute b - A*x */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans, -1.0, dcA, dcX, 1.0, dcB);

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);

    result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * m * eps ) ;

    if ( loud > 2 ) {
        printf("============\n");
        printf("Checking the Residual of the solution \n");
        if ( loud > 3 )
            printf( "-- ||A||_oo = %e, ||X||_oo = %e, ||B||_oo= %e, ||A X - B||_oo = %e\n",
                    Anorm, Xnorm, Bnorm, Rnorm );
        printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);
    }

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        if( loud ) printf("-- Solution is suspicious ! \n");
        info_solution = 1;
    }
    else{
        if( loud ) printf("-- Solution is CORRECT ! \n");
        info_solution = 0;
    }

    return info_solution;
}

static int check_inverse( parsec_context_t *parsec, int loud,
                          parsec_tiled_matrix_t *dcA,
                          parsec_tiled_matrix_t *dcInvA,
                          parsec_tiled_matrix_t *dcI )
{
    int info_solution;
    double Anorm    = 0.0;
    double InvAnorm = 0.0;
    double Rnorm, result;
    int m = dcA->m;
    double eps = LAPACKE_dlamch_work('e');

    Anorm    = dplasma_zlange(parsec, dplasmaInfNorm, dcA   );
    InvAnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcInvA);

    /* Compute I - A*A^{-1} */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans, -1.0, dcA, dcInvA, 1.0, dcI);

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcI);

    result = Rnorm / ( ( Anorm * InvAnorm ) * m * eps ) ;

    if ( loud > 2 ) {
        printf("============\n");
        printf("Checking the Residual of the solution \n");
        if ( loud > 3 )
            printf( "-- ||A||_oo = %e, ||A^{-1}||_oo = %e, ||A A^{-1} - I||_oo = %e\n",
                    Anorm, InvAnorm, Rnorm );
        printf("-- ||AA^{-1}-I||_oo/((||A||_oo||A^{-1}||_oo).N.eps) = %e \n", result);
    }

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        if( loud ) printf("-- Solution is suspicious ! \n");
        info_solution = 1;
    }
    else{
        if( loud ) printf("-- Solution is CORRECT ! \n");
        info_solution = 0;
    }

    return info_solution;
}
