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
                           dplasma_enum_t side, dplasma_enum_t uplo, dplasma_enum_t trans, dplasma_enum_t diag,
                           dplasma_complex64_t alpha,
                           int Am, int An, int Aseed,
                           int M,  int N,  int Cseed,
                           parsec_matrix_block_cyclic_t *dcCfinal );

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;
    int Aseed = 3872;
    int Cseed = 2873;
    dplasma_complex64_t alpha = 3.5;
    parsec_tiled_matrix_t *dcA;

#if defined(PRECISION_z) || defined(PRECISION_c)
    alpha -= I * 4.2;
#endif

    /* Set defaults for non argv iparams */
    iparam_default_gemm(iparam);
    iparam_default_ibnbmb(iparam, 0, 200, 200);
    iparam[IPARAM_NGPUS] = DPLASMA_ERR_NOT_SUPPORTED;

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    /* initializing matrix structure */
    int Am = max(M, N);
    LDA = max(LDA, Am);
    LDC = max(LDC, M);
    PASTE_CODE_ALLOCATE_MATRIX(dcA0, 1,
        parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, Am, 0, 0,
                               Am, Am, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    if(!check)
    {
        dplasma_enum_t side  = dplasmaLeft;
        dplasma_enum_t uplo  = dplasmaLower;
        dplasma_enum_t trans = dplasmaNoTrans;
        dplasma_enum_t diag  = dplasmaUnit;

        PASTE_CODE_FLOPS(FLOPS_ZTRMM, (side, (DagDouble_t)M, (DagDouble_t)N));

        /* Make A square */
        if (side == dplasmaLeft) {
            dcA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA0, 0, 0, M, M );
        } else {
            dcA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA0, 0, 0, N, N );
        }

        for(int t = 0; t < nruns+1; t++) {
            /* matrix generation */
            if(loud > 2) printf("+++ Generate matrices ... ");
            dplasma_zplghe( parsec, 0., uplo, dcA, Aseed);
            dplasma_zplrnt( parsec, 0,        (parsec_tiled_matrix_t *)&dcC, Cseed);
            if(loud > 2) printf("Done\n");

            parsec_devices_release_memory();
            /* Create PaRSEC */
            PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec, ztrmm,
                                                        (side, uplo, trans, diag,
                                                                1.0, dcA,
                                                                (parsec_tiled_matrix_t *)&dcC),
                                                                dplasma_ztrmm_Destruct( PARSEC_ztrmm ), t);

            parsec_devices_reset_load(parsec);
        }
        PASTE_CODE_PERF_LOOP_DONE();
        free(dcA);
    }
    else
    {
        int s, u, t, d;
        int info_solution;

        PASTE_CODE_ALLOCATE_MATRIX(dcC2, 1,
            parsec_matrix_block_cyclic, (&dcC2, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDC, N, 0, 0,
                                   M, N, P, nodes/P, KP, KQ, IP, JQ));

        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC2, Cseed);

        for (s=0; s<2; s++) {
            /* Make A square */
            if (side[s] == dplasmaLeft) {
                Am = M;
                dcA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA0, 0, 0, M, M );
            } else {
                Am = N;
                dcA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA0, 0, 0, N, N );
            }
            dplasma_zplghe( parsec, 0., dplasmaUpperLower, dcA, Aseed);
            for (u=0; u<2; u++) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                for (t=0; t<3; t++) {
#else
                for (t=0; t<2; t++) {
#endif
                    for (d=0; d<2; d++) {

                        if ( rank == 0 ) {
                            printf("***************************************************\n");
                            printf(" ----- TESTING ZTRMM (%s, %s, %s, %s) -------- \n",
                                   sidestr[s], uplostr[u], transstr[t], diagstr[d]);
                        }

                        /* matrix generation */
                        printf("Generate matrices ... ");
                        dplasma_zlacpy( parsec, dplasmaUpperLower,
                                        (parsec_tiled_matrix_t *)&dcC2, (parsec_tiled_matrix_t *)&dcC );
                        printf("Done\n");

                        /* Compute */
                        printf("Compute ... ... ");
                        dplasma_ztrmm(parsec, side[s], uplo[u], trans[t], diag[d],
                                      alpha, dcA, (parsec_tiled_matrix_t *)&dcC);
                        printf("Done\n");

                        /* Check the solution */
                        info_solution = check_solution(parsec, rank == 0 ? loud : 0,
                                                       side[s], uplo[u], trans[t], diag[d],
                                                       alpha, Am, Am, Aseed,
                                                              M,  N,  Cseed,
                                                       &dcC);
                        if ( rank == 0 ) {
                            if (info_solution == 0) {
                                printf(" ---- TESTING ZTRMM (%s, %s, %s, %s) ...... PASSED !\n",
                                       sidestr[s], uplostr[u], transstr[t], diagstr[d]);
                            }
                            else {
                                printf(" ---- TESTING ZTRMM (%s, %s, %s, %s) ... FAILED !\n",
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
        parsec_data_free(dcC2.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC2);
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
                           dplasma_enum_t side, dplasma_enum_t uplo, dplasma_enum_t trans, dplasma_enum_t diag,
                           dplasma_complex64_t alpha,
                           int Am, int An, int Aseed,
                           int M,  int N,  int Cseed,
                           parsec_matrix_block_cyclic_t *dcCfinal )
{
    int info_solution = 1;
    double Anorm, Cinitnorm, Cdplasmanorm, Clapacknorm, Rnorm;
    double eps, result;
    int MB = dcCfinal->super.mb;
    int NB = dcCfinal->super.nb;
    int LDA = Am;
    int LDC = M;
    int rank  = dcCfinal->super.super.myrank;

    eps = LAPACKE_dlamch_work('e');

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDA, An, 0, 0,
                               Am, An, 1, 1, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, 1, 1, 1, 1, 0, 0));

    dplasma_zplghe( parsec, 0., dplasmaUpperLower, (parsec_tiled_matrix_t *)&dcA, Aseed);
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed );

    Anorm        = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcA );
    Cinitnorm    = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcC );
    Cdplasmanorm = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)dcCfinal );

    if ( rank == 0 ) {
        cblas_ztrmm(CblasColMajor,
                    (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                    (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                    M, N,
                    CBLAS_SADDR(alpha), dcA.mat, LDA,
                                        dcC.mat, LDC );
    }

    Clapacknorm = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcC );

    dplasma_zgeadd( parsec, dplasmaNoTrans,
                    -1.0, (parsec_tiled_matrix_t*)dcCfinal,
                     1.0, (parsec_tiled_matrix_t*)&dcC );

    Rnorm = dplasma_zlange( parsec, dplasmaMaxNorm, (parsec_tiled_matrix_t*)&dcC );

    result = Rnorm / (Clapacknorm * max(M,N) * eps);

    if ( rank == 0 ) {
        if ( loud > 2 ) {
            printf("  ||A||_inf = %e, ||C||_inf = %e\n"
                   "  ||lapack(a*A*C)||_inf = %e, ||dplasma(a*A*C)||_inf = %e, ||R||_m = %e, res = %e\n",
                   Anorm, Cinitnorm, Clapacknorm, Cdplasmanorm, Rnorm, result);
        }

        if (  isinf(Clapacknorm) || isinf(Cdplasmanorm) ||
              isnan(result) || isinf(result) || (result > 10.0) ) {
            info_solution = 1;
        }
        else {
            info_solution = 0;
        }
    }

#if defined(PARSEC_HAVE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    parsec_data_free(dcC.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC);

    return info_solution;
}
