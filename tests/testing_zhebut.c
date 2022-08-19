/*
 * Copyright (c) 2009-2014 The University of Tennessee and The University
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

#if defined(CHECK_B)
static int check_solution( parsec_context_t *parsec, int loud, dplasma_enum_t uplo,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX );
#endif

static int check_inverse( parsec_context_t *parsec, int loud, dplasma_enum_t uplo,
                          parsec_tiled_matrix_t *dcA,
                          parsec_tiled_matrix_t *dcInvA,
                          parsec_tiled_matrix_t *dcI );

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;
    dplasma_enum_t uplo = dplasmaLower;
    dplasma_complex64_t *U_but_vec;
    DagDouble_t time_butterfly, time_facto, time_total;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 0, 180, 180);

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZHETRF, ((DagDouble_t)N));

    /* initializing matrix structure */
    LDA = max( LDA, N );
    LDB = max( LDB, N );
    KP = 1;
    KQ = 1;

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   N, N, P, nodes/P, uplo));

    PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
        parsec_matrix_sym_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   N, N, P, nodes/P, uplo));

#if defined(CHECK_B)
    PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
        parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, NRHS, 0, 0,
                                   N, NRHS, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
        parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, NRHS, 0, 0,
                                   N, NRHS, P, nodes/P, KP, KQ, IP, JQ));
#endif

    PASTE_CODE_ALLOCATE_MATRIX(dcInvA, check_inv,
        parsec_matrix_block_cyclic, (&dcInvA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               N, N, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcI, check_inv,
        parsec_matrix_block_cyclic, (&dcI, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               N, N, P, nodes/P, KP, KQ, IP, JQ));

    for(int t = 0; t < nruns+1; t++) {
        if(loud > 2 && rank==0) printf("+++ Generate matrices ... ");
        dplasma_zplghe( parsec, (double)(0), uplo,
                        (parsec_tiled_matrix_t *)&dcA, 1358);
        if(0 ==t && check) {
            dplasma_zlacpy( parsec, uplo,
                            (parsec_tiled_matrix_t *)&dcA,
                            (parsec_tiled_matrix_t *)&dcA0);
        }
        if(loud > 2 && rank == 0) printf("Done\n");

        if(loud > 2 && rank == 0) printf("+++ Computing Butterfly ... \n");
        SYNC_TIME_START();
        if (loud) TIME_START();

        ret = dplasma_zhebut(parsec, (parsec_tiled_matrix_t *)&dcA, &U_but_vec, butterfly_level);
        if( ret < 0 )
            return ret;

        if(t > 0) {
            SYNC_TIME_PRINT(rank, ("zhebut computation N= %d NB= %d\n", N, NB));

            /* Backup butterfly time */
            time_butterfly = sync_time_elapsed;
        }
        if(loud > 2) printf("... Done\n");


        if(loud > 2) printf("+++ Computing Factorization ... \n");
        SYNC_TIME_START();
        dplasma_zhetrf(parsec, (parsec_tiled_matrix_t *)&dcA);

        if(t > 0) {
            SYNC_TIME_PRINT(rank, ("zhetrf computation N= %d NB= %d : %f gflops\n", N, NB,
                    (gflops = (flops/1e9)/(sync_time_elapsed))));
            gflops_avg += gflops/nruns;
        }
        if(loud > 2) printf(" ... Done.\n");

        time_facto = sync_time_elapsed;
        time_total = time_butterfly + time_facto;

        if(0 == rank && t>0) {
            printf( "zhebut+zhetrf computation : %f gflops (%02.2f%% / %02.2f%%)\n",
                    (gflops = (flops/1e9)/(time_total)),
                    time_butterfly / time_total * 100.,
                    time_facto     / time_total * 100.);
        }
        (void)gflops;
    }
    PASTE_CODE_PERF_LOOP_DONE();

    if(check) {
        /* matrix generation */
        if(loud > 2) printf("+++ Generate additional checking matrices ... ");
#if defined(CHECK_B)
        dplasma_zplrnt( parsec, 0,
                        (parsec_tiled_matrix_t *)&dcB, 3872);

        dplasma_zlacpy( parsec, dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcB,
                        (parsec_tiled_matrix_t *)&dcX);
#endif
        if (check_inv) {
            dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&dcInvA);
            dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&dcI);
        }
        if(loud > 2) printf("Done\n");

#if defined(CHECK_B)
        dplasma_zhetrs(parsec, uplo,
                       (parsec_tiled_matrix_t *)&dcA,
                       (parsec_tiled_matrix_t *)&dcX,
                       U_but_vec, butterfly_level);

        /* Check the solution */
        ret |= check_solution( parsec, (rank == 0) ? loud : 0, uplo,
                               (parsec_tiled_matrix_t *)&dcA0,
                               (parsec_tiled_matrix_t *)&dcB,
                               (parsec_tiled_matrix_t *)&dcX);

        parsec_data_free(dcB.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);

        parsec_data_free(dcX.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcX);
#endif

        if (check_inv) {

            fprintf(stderr, "Start hetrs\n");
            dplasma_zhetrs(parsec, uplo,
                           (parsec_tiled_matrix_t *)&dcA,
                           (parsec_tiled_matrix_t *)&dcInvA,
                           U_but_vec, butterfly_level);

            fprintf(stderr, "Start check_inv\n");
            /* Check the solution against the inverse */
            ret |= check_inverse(parsec, (rank == 0) ? loud : 0, uplo,
                                 (parsec_tiled_matrix_t *)&dcA0,
                                 (parsec_tiled_matrix_t *)&dcInvA,
                                 (parsec_tiled_matrix_t *)&dcI);

            parsec_data_free(dcInvA.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcInvA);

            parsec_data_free(dcI.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcI);
        }
    }

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);

    cleanup_parsec(parsec, iparam);

    return ret;
}

#if defined(CHECK_B)
/*
 * This function destroys B
 */
static int check_solution( parsec_context_t *parsec, int loud, dplasma_enum_t uplo,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX )
{
    int info_solution;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double Bnorm = 0.0;
    double Xnorm, result;
    int N = dcB->m;
    double eps = LAPACKE_dlamch_work('e');

    Anorm = dplasma_zlanhe(parsec, dplasmaInfNorm, uplo, dcA);
    Bnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);
    Xnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcX);

    /* Compute A*x */
    dplasma_zhemm( parsec, dplasmaLeft, uplo, -1.0, dcA, dcX, 1.0, dcB);

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);

    result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * N * eps ) ;

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
#endif

static int check_inverse( parsec_context_t *parsec, int loud, dplasma_enum_t uplo,
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

    Anorm    = dplasma_zlanhe(parsec, dplasmaInfNorm, uplo, dcA);
    InvAnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcInvA);

    /* Compute I - A*A^{-1} */
    dplasma_zhemm( parsec, dplasmaLeft, uplo, -1.0, dcA, dcInvA, 1.0, dcI);

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
