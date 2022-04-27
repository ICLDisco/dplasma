/*
 * Copyright (c) 2009-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> z c d s
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"

static int check_solution( parsec_context_t *parsec, int loud,
                           dplasma_enum_t side, dplasma_enum_t uplo,
                           dplasma_complex64_t alpha, int Am, int An, int Aseed,
                           dplasma_complex64_t beta,  int M,  int N,  int Bseed, int Cseed,
                           parsec_matrix_block_cyclic_t *dcCfinal );

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;
    int Aseed = 3872;
    int Bseed = 4674;
    int Cseed = 2873;
    dplasma_complex64_t alpha = 3.5;
    dplasma_complex64_t beta  = -2.8;

#if defined(PRECISION_z) || defined(PRECISION_c)
    alpha -= I * 4.2;
    beta += I * 0.7;
#endif

    /* Set defaults for non argv iparams */
    iparam_default_gemm(iparam);
    iparam_default_ibnbmb(iparam, 0, 200, 200);
    iparam[IPARAM_NGPUS] = DPLASMA_ERR_NOT_SUPPORTED;

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);

    dplasma_warmup(parsec);

    LDB = max(LDB, M);
    LDC = max(LDC, M);

    PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
        parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDB, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcC2, check,
        parsec_matrix_block_cyclic, (&dcC2, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    if(!check)
    {
        dplasma_enum_t side  = dplasmaLeft;
        dplasma_enum_t uplo  = dplasmaLower;
        int Am = ( side == dplasmaLeft ? M : N );
        LDA = max(LDA, Am);

        PASTE_CODE_FLOPS(FLOPS_ZSYMM, (side, (DagDouble_t)M, (DagDouble_t)N));

        PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
            parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                       rank, MB, NB, LDA, Am, 0, 0,
                                       Am, Am, P, nodes/P, uplo));

        for(int t = 0; t < nruns; t++) {
            /* matrix generation */
            if (loud > 2) printf("Generate matrices ... ");
            dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB,  Bseed);
            dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC,  Cseed);
            dplasma_zplgsy( parsec, 0., uplo, (parsec_tiled_matrix_t *)&dcA, Aseed);
            if(loud > 2) printf("Done\n");

            /* Create PaRSEC */
            PASTE_CODE_ENQUEUE_KERNEL(parsec, zsymm,
                                      (side, uplo,
                                              alpha, (parsec_tiled_matrix_t *)&dcA,
                                              (parsec_tiled_matrix_t *)&dcB,
                                              beta,  (parsec_tiled_matrix_t *)&dcC));

            /* lets rock! */
            PASTE_CODE_PROGRESS_KERNEL(parsec, zsymm);

            dplasma_zsymm_Destruct( PARSEC_zsymm );
        }
        PASTE_CODE_PERF_LOOP_DONE();

        parsec_data_free(dcA.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    }
    else
    {
        int s, u;
        int info_solution;

        if (loud > 2) printf("Generate matrices ... ");
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB,  Bseed);
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC,  Cseed);
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC2, Cseed);

        for (s=0; s<2; s++) {
            /* initializing matrix structure */
            int Am = ( side[s] == dplasmaLeft ? M : N );
            LDA = max(LDA, Am);

            for (u=0; u<2; u++) {

                PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
                    parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                               rank, MB, NB, LDA, Am, 0, 0,
                                               Am, Am, P, nodes/P, uplo[u]));

                dplasma_zplgsy( parsec, 0., uplo[u], (parsec_tiled_matrix_t *)&dcA, Aseed);
                dplasma_zlacpy( parsec, dplasmaUpperLower,
                                (parsec_tiled_matrix_t *)&dcC2, (parsec_tiled_matrix_t *)&dcC );
                if (loud > 2) printf("Done\n");

                /* Compute */
                if (loud > 2) printf("Compute ... ... ");
                dplasma_zsymm(parsec, side[s], uplo[u],
                              (dplasma_complex64_t)alpha, (parsec_tiled_matrix_t *)&dcA,
                                                        (parsec_tiled_matrix_t *)&dcB,
                              (dplasma_complex64_t)beta,  (parsec_tiled_matrix_t *)&dcC);
                if (loud > 2) printf("Done\n");

                /* Check the solution */
                info_solution = check_solution(parsec, rank == 0 ? loud : 0,
                                               side[s], uplo[u],
                                               alpha, Am, Am, Aseed,
                                               beta,  M,  N,  Bseed, Cseed,
                                               &dcC);

                if ( rank == 0 ) {
                    if (info_solution == 0) {
                        printf(" ---- TESTING ZSYMM (%s, %s) ...... PASSED !\n",
                               sidestr[s], uplostr[u]);
                    }
                    else {
                        printf(" ---- TESTING ZSYMM (%s, %s) ... FAILED !\n",
                               sidestr[s], uplostr[u]);
                        ret |= 1;
                    }
                    printf("***************************************************\n");
                }

                parsec_data_free(dcA.mat);
                parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
            }
        }

        parsec_data_free(dcC2.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC2);
    }

    parsec_data_free(dcB.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
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
                           dplasma_enum_t side, dplasma_enum_t uplo,
                           dplasma_complex64_t alpha, int Am, int An, int Aseed,
                           dplasma_complex64_t beta,  int M,  int N,  int Bseed, int Cseed,
                           parsec_matrix_block_cyclic_t *dcCfinal )
{
    int info_solution = 1;
    double Anorm, Bnorm, Cinitnorm, Cdplasmanorm, Clapacknorm, Rnorm;
    double eps, result;
    int MB = dcCfinal->super.mb;
    int NB = dcCfinal->super.nb;
    int LDA = Am;
    int LDC = M;
    int LDB = LDC;
    int rank  = dcCfinal->super.super.myrank;

    eps = LAPACKE_dlamch_work('e');

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDA, An, 0, 0,
                               Am, An, 1, 1, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
        parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, 1, 1, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, 1, 1, 1, 1, 0, 0));

    dplasma_zplgsy( parsec, 0., dplasmaUpperLower, (parsec_tiled_matrix_t *)&dcA, Aseed);
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, Bseed );
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed );

    Anorm        = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcA );
    Bnorm        = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcB );
    Cinitnorm    = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcC );
    Cdplasmanorm = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)dcCfinal );

    if ( rank == 0 ) {
        cblas_zsymm(CblasColMajor,
                    (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                    M, N,
                    CBLAS_SADDR(alpha), dcA.mat, LDA,
                                        dcB.mat, LDB,
                    CBLAS_SADDR(beta),  dcC.mat, LDC);
    }

    Clapacknorm = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcC );

    dplasma_zgeadd( parsec, dplasmaNoTrans, -1.0, (parsec_tiled_matrix_t*)dcCfinal,
                                           1.0, (parsec_tiled_matrix_t*)&dcC );

    Rnorm = dplasma_zlange( parsec, dplasmaMaxNorm, (parsec_tiled_matrix_t*)&dcC );

    result = Rnorm / (Clapacknorm * max(M,N) * eps);

    if ( rank == 0 ) {
        if ( loud > 2 ) {
            printf("  ||A||_inf = %e, ||B||_inf = %e, ||C||_inf = %e\n"
                   "  ||lapack(a*A*B+b*C)||_inf = %e, ||dplasma(a*A*B+b*C)||_inf = %e, ||R||_m = %e, res = %e\n",
                   Anorm, Bnorm, Cinitnorm, Clapacknorm, Cdplasmanorm, Rnorm, result);
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
    parsec_data_free(dcB.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
    parsec_data_free(dcC.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC);

    return info_solution;
}
