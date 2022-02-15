/*
 * Copyright (c) 2009-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

static int check_solution( parsec_context_t *parsec, int loud,
                           int M, int N,
                           parsec_tiled_matrix_t *dcC,
                           parsec_tiled_matrix_t *dcCfinal );

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;
    int Aseed = 3872;
    int Cseed = 2873;
    dplasma_complex64_t alpha = 0.98;
    parsec_tiled_matrix_t *dcA;

#if defined(PRECISION_z) || defined(PRECISION_c)
    alpha -= I * 0.32;
#endif

    /* Set defaults for non argv iparams */
    iparam_default_gemm(iparam);
    iparam_default_ibnbmb(iparam, 0, 200, 200);

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    /* initializing matrix structure */
    int Am = dplasma_imax(M, N);
    LDA = dplasma_imax(LDA, Am);
    LDC = dplasma_imax(LDC, M);
    PASTE_CODE_ALLOCATE_MATRIX(dcA0, 1,
        parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, Am, 0, 0,
                               Am, Am, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcC0, check,
        parsec_matrix_block_cyclic, (&dcC0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    if(!check)
    {
        dplasma_enum_t side  = dplasmaLeft;
        dplasma_enum_t uplo  = dplasmaLower;
        dplasma_enum_t trans = dplasmaNoTrans;
        dplasma_enum_t diag  = dplasmaUnit;

        PASTE_CODE_FLOPS(FLOPS_ZTRSM, (side, (DagDouble_t)M, (DagDouble_t)N));

        for(int t = 0; t < nruns+1; t++) {
            /* matrix generation */
            if(loud > 2) printf("+++ Generate matrices ... ");
            /* Generate matrix A with diagonal dominance to keep stability during computation */
            dplasma_zplrnt( parsec, 1, (parsec_tiled_matrix_t *)&dcA0, Aseed);
            /* Scale down the full matrix to keep stability in diag = dplasmaUnit case */
            dplasma_zlascal( parsec, dplasmaUpperLower,
                             1. / (dplasma_complex64_t)Am,
                             (parsec_tiled_matrix_t *)&dcA0 );
            dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed);
            if (check)
                dplasma_zlacpy( parsec, dplasmaUpperLower,
                                (parsec_tiled_matrix_t *)&dcC, (parsec_tiled_matrix_t *)&dcC0 );
            if(loud > 2) printf("Done\n");

            /* Make A square */
            if (side == dplasmaLeft) {
                dcA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA0, 0, 0, M, M );
            } else {
                dcA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA0, 0, 0, N, N );
            }

            /* Compute b = 1/alpha * A * x */
            dplasma_ztrmm(parsec, side, uplo, trans, diag, 1. / alpha,
                          dcA, (parsec_tiled_matrix_t *)&dcC);

            /* Create PaRSEC */
            PASTE_CODE_ENQUEUE_KERNEL(parsec, ztrsm,
                                      (side, uplo, trans, diag, alpha,
                                              dcA, (parsec_tiled_matrix_t *)&dcC));

            /* lets rock! */
            PASTE_CODE_PROGRESS_KERNEL(parsec, ztrsm, t);

            dplasma_ztrsm_Destruct( PARSEC_ztrsm );
        }
        PASTE_CODE_PERF_LOOP_DONE();

        free(dcA);
    }
    else
    {
        int s, u, t, d;
        int info_solution;

        /* matrix generation */
        if(loud > 2) printf("+++ Generate matrices ... ");
        /* Generate matrix A with diagonal dominance to keep stability during computation */
        dplasma_zplrnt( parsec, 1, (parsec_tiled_matrix_t *)&dcA0, Aseed);
        /* Scale down the full matrix to keep stability in diag = dplasmaUnit case */
        dplasma_zlascal( parsec, dplasmaUpperLower,
                         1. / (dplasma_complex64_t)Am,
                         (parsec_tiled_matrix_t *)&dcA0 );
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed);
        dplasma_zlacpy( parsec, dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcC, (parsec_tiled_matrix_t *)&dcC0 );
        if(loud > 2) printf("Done\n");

        for (s=0; s<2; s++) {
            /* Make A square */
            if (side[s] == dplasmaLeft) {
                dcA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA0, 0, 0, M, M );
            } else {
                dcA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA0, 0, 0, N, N );
            }

            for (u=0; u<2; u++) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                for (t=0; t<3; t++) {
#else
                for (t=0; t<2; t++) {
#endif
                    for (d=0; d<2; d++) {

                        if ( rank == 0 ) {
                            printf("***************************************************\n");
                            printf(" ----- TESTING ZTRSM (%s, %s, %s, %s) -------- \n",
                                   sidestr[s], uplostr[u], transstr[t], diagstr[d]);
                        }

                        /* matrix generation */
                        printf("Generate matrices ... ");
                        dplasma_zlacpy( parsec, dplasmaUpperLower,
                                        (parsec_tiled_matrix_t *)&dcC0,
                                        (parsec_tiled_matrix_t *)&dcC );
                        dplasma_ztrmm(parsec, side[s], uplo[u], trans[t], diag[d], 1./alpha,
                                      dcA, (parsec_tiled_matrix_t *)&dcC);
                        printf("Done\n");

                        /* Compute */
                        printf("Compute ... ... ");
                        dplasma_ztrsm(parsec, side[s], uplo[u], trans[t], diag[d], alpha,
                                      dcA, (parsec_tiled_matrix_t *)&dcC);
                        printf("Done\n");

                        /* Check the solution */
                        info_solution = check_solution(parsec, rank == 0 ? loud : 0, M, N,
                                                       (parsec_tiled_matrix_t*)&dcC0,
                                                       (parsec_tiled_matrix_t*)&dcC);
                        if ( rank == 0 ) {
                            if (info_solution == 0) {
                                printf(" ---- TESTING ZTRSM (%s, %s, %s, %s) ...... PASSED !\n",
                                       sidestr[s], uplostr[u], transstr[t], diagstr[d]);
                            }
                            else {
                                printf(" ---- TESTING ZTRSM (%s, %s, %s, %s) ... FAILED !\n",
                                       sidestr[s], uplostr[u], transstr[t], diagstr[d]);
                                ret |= 1;
                            }
                            printf("***************************************************\n");
                        }
                    }
                }
            }
            free(dcA);
        }
        parsec_data_free(dcC0.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC0);
    }

    parsec_data_free(dcA0.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0);
    parsec_data_free(dcC.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC);

    cleanup_parsec(parsec, iparam);

    return ret;
}


/**********************************
 * static functions
 **********************************/

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
static int check_solution( parsec_context_t *parsec, int loud,
                           int M, int N,
                           parsec_tiled_matrix_t *dcC,
                           parsec_tiled_matrix_t *dcCfinal )
{
    int info_solution = 1;
    double Cinitnorm, Cdplasmanorm, Rnorm;
    double eps, result;

    eps = LAPACKE_dlamch_work('e');

    Cinitnorm    = dplasma_zlange( parsec, dplasmaInfNorm, dcC );
    Cdplasmanorm = dplasma_zlange( parsec, dplasmaInfNorm, dcCfinal );

    dplasma_zgeadd( parsec, dplasmaNoTrans,
                    -1.0, (parsec_tiled_matrix_t*)dcC,
                     1.0, (parsec_tiled_matrix_t*)dcCfinal );

    Rnorm = dplasma_zlange( parsec, dplasmaMaxNorm, dcCfinal );

    result = Rnorm / (Cinitnorm * eps * dplasma_imax(M, N));

    if ( loud > 2 ) {
        printf("  ||x||_inf = %e, ||dplasma(A^(-1) b||_inf = %e, ||R||_m = %e, res = %e\n",
               Cinitnorm, Cdplasmanorm, Rnorm, result);
    }

    if (  isinf(Cdplasmanorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    return info_solution;
}
