/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

static int check_orthogonality(parsec_context_t *parsec, int loud,
                               parsec_tiled_matrix_t *Q);
static int check_factorization(parsec_context_t *parsec, int loud,
                               parsec_tiled_matrix_t *Aorig,
                               parsec_tiled_matrix_t *A,
                               parsec_tiled_matrix_t *Q);
static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX );

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;
    dplasma_qrtree_t qrtree;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 32, 200, 200);
    iparam[IPARAM_KP] = 1;
    iparam[IPARAM_KQ] = 1;
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'n';
    iparam[IPARAM_QR_HLVL_SZE] = -'Q';

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);

    /* Make sure KP and KQ are set to 1, since it conflicts with HQR */
    iparam[IPARAM_KP] = 1;
    iparam[IPARAM_KQ] = 1;

    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZGELQF, ((DagDouble_t)M, (DagDouble_t)N));

    LDA = max(M, LDA);
    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcTS, 1,
        parsec_matrix_block_cyclic, (&dcTS, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcTT, 1,
        parsec_matrix_block_cyclic, (&dcTT, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
        parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcQ, check,
        parsec_matrix_block_cyclic, (&dcQ, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    /* Check the solution */
    PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
        parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDB, NRHS, 0, 0,
                               N, NRHS, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
        parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDB, NRHS, 0, 0,
                               N, NRHS, P, nodes/P, KP, KQ, IP, JQ));

    for(int t = 0; t < nruns+1; t++) {
        if(loud > 3) printf("+++ Generate matrices ... ");
        dplasma_zpltmg( parsec, matrix_init, (parsec_tiled_matrix_t *)&dcA, random_seed );
        dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0., (parsec_tiled_matrix_t *)&dcTS);
        dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0., (parsec_tiled_matrix_t *)&dcTT);
        if(t == 0 && check)
            dplasma_zlacpy( parsec, dplasmaUpperLower,
                            (parsec_tiled_matrix_t *)&dcA, (parsec_tiled_matrix_t *)&dcA0 );
        if(loud > 3) printf("Done\n");

        dplasma_hqr_init( &qrtree,
                          dplasmaConjTrans, (parsec_tiled_matrix_t *)&dcA,
                          iparam[IPARAM_LOWLVL_TREE], iparam[IPARAM_HIGHLVL_TREE],
                          iparam[IPARAM_QR_TS_SZE],   iparam[IPARAM_QR_HLVL_SZE],
                          iparam[IPARAM_QR_DOMINO],   iparam[IPARAM_QR_TSRR] );

        /* Create PaRSEC */
        PASTE_CODE_ENQUEUE_KERNEL(parsec, zgelqf_param,
                                  (&qrtree,
                                          (parsec_tiled_matrix_t*)&dcA,
                                          (parsec_tiled_matrix_t*)&dcTS,
                                          (parsec_tiled_matrix_t*)&dcTT));

        /* lets rock! This code should be copy the PASTE_CODE_PROGRESS_KERNEL macro */
        SYNC_TIME_START();
        parsec_context_start(parsec);
        TIME_START();
        parsec_context_wait(parsec);

        if(t > 0) {
            SYNC_TIME_PRINT(rank,
                            ("zgelqf HQR computation NP= %d NC= %d P= %d IB= %d MB= %d NB= %d qr_a= %d qr_p = %d treel= %d treeh= %d domino= %d RR= %d M= %d N= %d : %f gflops\n",
                                    iparam[IPARAM_NNODES],
                                    iparam[IPARAM_NCORES],
                                    iparam[IPARAM_P],
                                    iparam[IPARAM_IB],
                                    iparam[IPARAM_MB],
                                    iparam[IPARAM_NB],
                                    iparam[IPARAM_QR_TS_SZE],
                                    iparam[IPARAM_QR_HLVL_SZE],
                                    iparam[IPARAM_LOWLVL_TREE],
                                    iparam[IPARAM_HIGHLVL_TREE],
                                    iparam[IPARAM_QR_DOMINO],
                                    iparam[IPARAM_QR_TSRR],
                                    iparam[IPARAM_M],
                                    iparam[IPARAM_N],
                                    gflops = (flops/1e9)/(sync_time_elapsed)));
            gflops_avg += gflops/nruns;
        }

        dplasma_zgelqf_param_Destruct( PARSEC_zgelqf_param );
    }
    PASTE_CODE_PERF_LOOP_DONE();

#if defined(PARSEC_SIM)
    if ( rank == 0 ) {
        printf("zgelqf HQR simulation NP= %d NC= %d P= %d qr_a= %d qr_p = %d treel= %d treeh= %d domino= %d RR= %d MT= %d NT= %d : %d \n",
               iparam[IPARAM_NNODES],
               iparam[IPARAM_NCORES],
               iparam[IPARAM_P],
               iparam[IPARAM_QR_TS_SZE],
               iparam[IPARAM_QR_HLVL_SZE],
               iparam[IPARAM_LOWLVL_TREE],
               iparam[IPARAM_HIGHLVL_TREE],
               iparam[IPARAM_QR_DOMINO],
               iparam[IPARAM_QR_TSRR],
               MT, NT,
               parsec_getsimulationdate( parsec ));
    }
#endif

    if( check ) {
        if (N >= M) {
            if(loud > 2) printf("+++ Generate the Q ...");
            dplasma_zunglq_param( parsec, &qrtree,
                                  (parsec_tiled_matrix_t *)&dcA,
                                  (parsec_tiled_matrix_t *)&dcTS,
                                  (parsec_tiled_matrix_t *)&dcTT,
                                  (parsec_tiled_matrix_t *)&dcQ);
            if(loud > 2) printf("Done\n");

            if(loud > 2) printf("+++ Solve the system ...");
            dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcX, random_seed+1);
            dplasma_zlacpy( parsec, dplasmaUpperLower,
                            (parsec_tiled_matrix_t *)&dcX,
                            (parsec_tiled_matrix_t *)&dcB );
            dplasma_zgelqs_param( parsec, &qrtree,
                                  (parsec_tiled_matrix_t *)&dcA,
                                  (parsec_tiled_matrix_t *)&dcTS,
                                  (parsec_tiled_matrix_t *)&dcTT,
                                  (parsec_tiled_matrix_t *)&dcX );
            if(loud > 2) printf("Done\n");

            /* Check the orthogonality, factorization and the solution */
            ret |= check_orthogonality( parsec, (rank == 0) ? loud : 0,
                                        (parsec_tiled_matrix_t *)&dcQ);
            ret |= check_factorization( parsec, (rank == 0) ? loud : 0,
                                        (parsec_tiled_matrix_t *)&dcA0,
                                        (parsec_tiled_matrix_t *)&dcA,
                                        (parsec_tiled_matrix_t *)&dcQ );
            ret |= check_solution( parsec, (rank == 0) ? loud : 0,
                                   (parsec_tiled_matrix_t *)&dcA0,
                                   (parsec_tiled_matrix_t *)&dcB,
                                   (parsec_tiled_matrix_t *)&dcX );

        } else {
            printf("Check cannot be performed when M > N\n");
        }

        parsec_data_free(dcA0.mat);
        parsec_data_free(dcQ.mat);
        parsec_data_free(dcB.mat);
        parsec_data_free(dcX.mat);
        parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&dcA0);
        parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&dcQ);
        parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&dcB);
        parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&dcX);
    }

    dplasma_hqr_finalize( &qrtree );

    parsec_data_free(dcA.mat);
    parsec_data_free(dcTS.mat);
    parsec_data_free(dcTT.mat);
    parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&dcA);
    parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&dcTS);
    parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&dcTT);

    cleanup_parsec(parsec, iparam);

    return ret;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int check_orthogonality(parsec_context_t *parsec, int loud, parsec_tiled_matrix_t *Q)
{
    parsec_matrix_block_cyclic_t *twodQ = (parsec_matrix_block_cyclic_t *)Q;
    double normQ = 999999.0;
    double result;
    double eps = LAPACKE_dlamch_work('e');
    int info_ortho;
    int M = Q->m;
    int N = Q->n;
    int minMN = min(M, N);

    PASTE_CODE_ALLOCATE_MATRIX(Id, 1,
        parsec_matrix_block_cyclic, (&Id, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               twodQ->grid.rank,
                               Q->mb, Q->nb, minMN, minMN, 0, 0,
                               minMN, minMN, twodQ->grid.rows, twodQ->grid.cols, twodQ->grid.krows, twodQ->grid.kcols, twodQ->grid.ip, twodQ->grid.jq));

    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&Id);

    /* Perform Id - Q'Q */
    if ( M >= N ) {
        dplasma_zherk( parsec, dplasmaUpper, dplasmaConjTrans,
                       1.0, Q, -1.0, (parsec_tiled_matrix_t*)&Id );
    } else {
        dplasma_zherk( parsec, dplasmaUpper, dplasmaNoTrans,
                       1.0, Q, -1.0, (parsec_tiled_matrix_t*)&Id );
    }

    normQ = dplasma_zlanhe(parsec, dplasmaInfNorm, dplasmaUpper, (parsec_tiled_matrix_t*)&Id);

    result = normQ / (minMN * eps);
    if ( loud ) {
        printf("============\n");
        printf("Checking the orthogonality of Q \n");
        printf("||Id-Q'*Q||_oo / (M*eps) = %e \n", result);
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        if ( loud ) printf("-- Orthogonality is suspicious ! \n");
        info_ortho=1;
    }
    else {
        if ( loud ) printf("-- Orthogonality is CORRECT ! \n");
        info_ortho=0;
    }

    parsec_data_free(Id.mat);
    parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&Id);
    return info_ortho;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int
check_factorization(parsec_context_t *parsec, int loud,
                    parsec_tiled_matrix_t *Aorig,
                    parsec_tiled_matrix_t *A,
                    parsec_tiled_matrix_t *Q)
{
    parsec_tiled_matrix_t *subA;
    parsec_matrix_block_cyclic_t *twodA = (parsec_matrix_block_cyclic_t *)A;
    double Anorm, Rnorm;
    double result;
    double eps = LAPACKE_dlamch_work('e');
    int info_factorization;
    int M = A->m;
    int N = A->n;
    int minMN = min(M, N);

    PASTE_CODE_ALLOCATE_MATRIX(Residual, 1,
        parsec_matrix_block_cyclic, (&Residual, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               twodA->grid.rank,
                               A->mb, A->nb, M, N, 0, 0,
                               M, N, twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols, twodA->grid.ip, twodA->grid.jq));

    PASTE_CODE_ALLOCATE_MATRIX(L, 1,
        parsec_matrix_block_cyclic, (&L, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               twodA->grid.rank,
                               A->mb, A->nb, M, M, 0, 0,
                               M, M, twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols, twodA->grid.ip, twodA->grid.jq));

    /* Copy the original A in Residual */
    dplasma_zlacpy( parsec, dplasmaUpperLower, Aorig, (parsec_tiled_matrix_t *)&Residual );

    /* Extract the L */
    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0., (parsec_tiled_matrix_t *)&L);

    subA = parsec_tiled_matrix_submatrix( A, 0, 0, M, M );
    dplasma_zlacpy( parsec, dplasmaLower, subA, (parsec_tiled_matrix_t *)&L );
    free(subA);

    /* Perform Residual = Aorig - L*Q */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans,
                   -1.0, (parsec_tiled_matrix_t *)&L, Q,
                    1.0, (parsec_tiled_matrix_t *)&Residual);

    /* Free R */
    parsec_data_free(L.mat);
    parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&L);

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&Residual);
    Anorm = dplasma_zlange(parsec, dplasmaInfNorm, Aorig);

    result = Rnorm / ( Anorm * minMN * eps);

    if ( loud ) {
        printf("============\n");
        printf("Checking the LQ Factorization \n");
        printf("-- ||A-LQ||_oo/(||A||_oo.M.eps) = %e \n", result );
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        if ( loud ) printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else {
        if ( loud ) printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    parsec_data_free(Residual.mat);
    parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&Residual);
    return info_factorization;
}

static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX )
{
    parsec_tiled_matrix_t *subB;
    int info_solution;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double Bnorm = 0.0;
    double Xnorm, result;
    double eps = LAPACKE_dlamch_work('e');

    subB = parsec_tiled_matrix_submatrix( dcB, 0, 0, dcA->m, dcB->n );

    Anorm = dplasma_zlange(parsec, dplasmaInfNorm, dcA);
    Bnorm = dplasma_zlange(parsec, dplasmaInfNorm, subB);
    Xnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcX);

    /* Compute A*x-b */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans, 1.0, dcA, dcX, -1.0, subB);

    /* Compute A' * ( A*x - b ) */
    dplasma_zgemm( parsec, dplasmaConjTrans, dplasmaNoTrans,
                   1.0, dcA, subB, 0., dcX );

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcX );
    free(subB);

    result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * dcA->m * eps ) ;

    if ( loud > 2 ) {
        printf("============\n");
        printf("Checking the Residual of the solution \n");
        if ( loud > 3 )
            printf( "-- ||A||_oo = %e, ||X||_oo = %e, ||B||_oo= %e, ||A X - B||_oo = %e\n",
                    Anorm, Xnorm, Bnorm, Rnorm );

        printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).M.eps) = %e \n", result);
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
