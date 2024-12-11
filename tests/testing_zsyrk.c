/*
 * Copyright (c) 2009-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> z c
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"

static int check_solution( parsec_context_t *parsec, int loud,
                           dplasma_enum_t uplo, dplasma_enum_t trans,
                           dplasma_complex64_t alpha, int Am, int An, int Aseed,
                           dplasma_complex64_t beta,  int M,  int N,  int Cseed,
                           parsec_matrix_sym_block_cyclic_t *dcCfinal );

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;
    int Aseed = 3872;
    int Cseed = 2873;
    dplasma_complex64_t alpha =  3.5 - I * 4.2;
    dplasma_complex64_t beta  = -2.8 + I * 0.7;

    /* Set defaults for non argv iparams */
    iparam_default_gemm(iparam);
    iparam_default_ibnbmb(iparam, 0, 200, 200);

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);

    M = N;
    LDC = max(LDC, N);

    if(!check)
    {
        dplasma_enum_t trans = dplasmaNoTrans;
        int Am = ( trans == dplasmaNoTrans ? N : K );
        int An = ( trans == dplasmaNoTrans ? K : N );
        LDA = max(LDA, Am);

        PASTE_CODE_FLOPS(FLOPS_ZSYRK, ((DagDouble_t)K, (DagDouble_t)N));

        PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
            parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDA, An, 0, 0,
                                   Am, An, P, nodes/P, KP, KQ, IP, JQ));

        PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
            parsec_matrix_sym_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                       rank, MB, NB, LDC, N, 0, 0,
                                       N, N, P, nodes/P, uplo));

        /* matrix generation */
        if(loud > 2) printf("+++ Generate matrices ... ");
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA,  Aseed);
        dplasma_zplgsy( parsec, 0., uplo, (parsec_tiled_matrix_t *)&dcC, Cseed);
        if(loud > 2) printf("Done\n");

        /* Create PaRSEC */
        PASTE_CODE_ENQUEUE_KERNEL(parsec, zsyrk,
                                  (uplo, trans,
                                   alpha, (parsec_tiled_matrix_t *)&dcA,
                                   beta,  (parsec_tiled_matrix_t *)&dcC));

        /* lets rock! */
        PASTE_CODE_PROGRESS_KERNEL(parsec, zsyrk);

        dplasma_zsyrk_Destruct( PARSEC_zsyrk );

        parsec_data_free(dcA.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
        parsec_data_free(dcC.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC);
    }
    else
    {
        int u, t;
        int info_solution;

        PASTE_CODE_ALLOCATE_MATRIX(dcC2, check,
            parsec_matrix_block_cyclic, (&dcC2, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDC, N, 0, 0,
                                   N, N, P, nodes/P, KP, KQ, IP, JQ));
        if (loud > 2) printf("Generate matrices ... ");
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC2, Cseed);
        if (loud > 2) printf("Done\n");

        for (u=0; u<2; u++) {
            PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
                parsec_matrix_sym_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                           rank, MB, NB, LDC, N, 0, 0,
                                           N, N, P, nodes/P, uplos[u]));

            for (t=0; t<2; t++) {

                /* initializing matrix structure */
                int Am = ( trans[t] == dplasmaNoTrans ? N : K );
                int An = ( trans[t] == dplasmaNoTrans ? K : N );
                LDA = max(LDA, Am);

                PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
                    parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                           rank, MB, NB, LDA, An, 0, 0,
                                           Am, An, P, nodes/P, KP, KQ, IP, JQ));

                if (loud > 2) printf("Generate matrices ... ");
                dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, Aseed);
                dplasma_zlacpy( parsec, uplos[u],
                                (parsec_tiled_matrix_t *)&dcC2, (parsec_tiled_matrix_t *)&dcC );
                if (loud > 2) printf("Done\n");

                /* Compute */
                if (loud > 2) printf("Compute ... ... ");
                dplasma_zsyrk(parsec, uplos[u], trans[t],
                              alpha, (parsec_tiled_matrix_t *)&dcA,
                              beta,  (parsec_tiled_matrix_t *)&dcC);
                if (loud > 2) printf("Done\n");

                /* Check the solution */
                info_solution = check_solution(parsec, rank == 0 ? loud : 0,
                                               uplos[u], trans[t],
                                               alpha, Am, An, Aseed,
                                               beta,  N,  N,  Cseed,
                                               &dcC);

                if ( rank == 0 ) {
                    if (info_solution == 0) {
                        printf(" ---- TESTING ZSYRK (%s, %s) ...... PASSED !\n",
                               uplostr[u], transstr[t]);
                    }
                    else {
                        printf(" ---- TESTING ZSYRK (%s, %s) ... FAILED !\n",
                               uplostr[u], transstr[t]);
                        ret |= 1;
                    }
                    printf("***************************************************\n");
                }

                parsec_data_free(dcA.mat);
                parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
            }
            parsec_data_free(dcC.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC);
        }

        parsec_data_free(dcC2.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC2);
    }

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
                           dplasma_enum_t uplo, dplasma_enum_t trans,
                           dplasma_complex64_t alpha, int Am, int An, int Aseed,
                           dplasma_complex64_t beta,  int M,  int N,  int Cseed,
                           parsec_matrix_sym_block_cyclic_t *dcCfinal )
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

    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, Aseed);
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed );

    Anorm        = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcA );
    Cinitnorm    = dplasma_zlansy( parsec, dplasmaInfNorm, uplo, (parsec_tiled_matrix_t*)&dcC     );
    Cdplasmanorm = dplasma_zlansy( parsec, dplasmaInfNorm, uplo, (parsec_tiled_matrix_t*)dcCfinal );

    if ( rank == 0 ) {
        cblas_zsyrk(CblasColMajor,
                    (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                    N, (trans == dplasmaNoTrans) ? An : Am,
                    CBLAS_SADDR(alpha), dcA.mat, LDA,
                    CBLAS_SADDR(beta),  dcC.mat, LDC);
    }

    Clapacknorm = dplasma_zlansy( parsec, dplasmaInfNorm, uplo, (parsec_tiled_matrix_t*)&dcC );

    dplasma_ztradd( parsec, uplo, dplasmaNoTrans,
                    -1.0, (parsec_tiled_matrix_t*)dcCfinal,
                     1.0, (parsec_tiled_matrix_t*)&dcC );

    Rnorm = dplasma_zlansy( parsec, dplasmaMaxNorm, uplo, (parsec_tiled_matrix_t*)&dcC );

    result = Rnorm / (Clapacknorm * max(M,N) * eps);

    if ( rank == 0 ) {
        if ( loud > 2 ) {
            printf("  ||A||_inf = %e, ||C||_inf = %e\n"
                   "  ||lapack(a*A*At+b*C)||_inf = %e, ||dplasma(a*A*At+b*C)||_inf = %e, ||R||_m = %e, res = %e\n",
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
