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
                           dplasma_enum_t transA, dplasma_enum_t transB,
                           dplasma_complex64_t alpha, int Am, int An, int Aseed,
                                                      int Bm, int Bn, int Bseed,
                           dplasma_complex64_t beta,  int M,  int N,  int Cseed,
                           parsec_matrix_block_cyclic_t *dcCfinal );

static void warmup_zgemm(int rank, int nodes, int random_seed, parsec_context_t *parsec)
{
    int MB = 64;
    int NB = 64;
    int KB = 64;
    int MT = nodes;
    int NT = 1;
    int KT = 1;
    int M = MT*MB;
    int N = NT*NB;
    int K = KT*KB;
    int did;
    unsigned int rs = (unsigned int)random_seed;
    int Aseed = rand_r(&rs);
    int Bseed = rand_r(&rs);
    int Cseed = rand_r(&rs);
    int tA = dplasmaNoTrans;
    int tB = dplasmaNoTrans;
    dplasma_complex64_t alpha =  0.51;
    dplasma_complex64_t beta  = -0.42;

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
            parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, KB, M, K, 0, 0,
                                   M, K, nodes, 1, 1, 1, 0, 0));

    PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
            parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, KB, NB, K, N, 0, 0,
                                   K, N, 1, 1, 1, 1, 0, 0));

    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
            parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, M, N, 0, 0,
                                   M, N, nodes, 1, 1, 1, 0, 0));

    /* Do the CPU warmup first */
    dplasma_zplrnt( parsec, 0, &dcA.super, Aseed);
    dplasma_zplrnt( parsec, 0, &dcB.super, Bseed);
    dplasma_zplrnt( parsec, 0, &dcC.super, Cseed);
    parsec_taskpool_t *zgemm = dplasma_zgemm_New(tA, tB, alpha, &dcA.super, &dcB.super, beta, &dcC.super);
    zgemm->devices_index_mask = 1<<0; /* Only CPU ! */
    parsec_context_add_taskpool(parsec, zgemm);
    parsec_context_start(parsec);
    parsec_context_wait(parsec);
    dplasma_zgemm_Destruct(zgemm);

    /* Now do the other devices, skipping RECURSIVE */
    /* We know that there is a GPU-enabled version of this operation, so warm it up if some device is enabled */
    for(did = 2; did < (int)parsec_nb_devices; did++) {
        for(int i = 0; i < MT; i++) {
            for(int j = 0; j < NT; j++) {
                if( rank == (int)dcC.super.super.rank_of(&dcC.super.super, i, j) ) {
                    parsec_data_t *dta = dcC.super.super.data_of(&dcC.super.super, i, j);
                    parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                }
            }
        }
        dplasma_zplrnt( parsec, 0, &dcA.super, Aseed);
        dplasma_zplrnt( parsec, 0, &dcB.super, Bseed);
        dplasma_zplrnt( parsec, 0, &dcC.super, Cseed);
        dplasma_zgemm(parsec, tA, tB, alpha, &dcA.super, &dcB.super, beta, &dcC.super);
        parsec_devices_release_memory();
    }

    parsec_data_free(dcA.mat); dcA.mat = NULL;
    parsec_tiled_matrix_destroy( &dcA.super );   
    parsec_data_free(dcB.mat); dcB.mat = NULL;
    parsec_tiled_matrix_destroy( &dcB.super );   
    parsec_data_free(dcC.mat); dcC.mat = NULL;
    parsec_tiled_matrix_destroy( &dcC.super );   
}

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int info_solution = 0;
    int Aseed = 3872;
    int Bseed = 4674;
    int Cseed = 2873;
    int tA = dplasmaNoTrans;
    int tB = dplasmaNoTrans;
    dplasma_complex64_t alpha =  0.51;
    dplasma_complex64_t beta  = -0.42;

#if defined(PRECISION_z) || defined(PRECISION_c)
    alpha -= I * 0.32;
    beta  += I * 0.21;
#endif

    /* Set defaults for non argv iparams */
    iparam_default_gemm(iparam);
    iparam_default_ibnbmb(iparam, 0, 200, 200);

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);

    PASTE_CODE_FLOPS(FLOPS_ZGEMM, ((DagDouble_t)M,(DagDouble_t)N,(DagDouble_t)K));

    LDA = max(LDA, max(M, K));
    LDB = max(LDB, max(K, N));
    LDC = max(LDC, M);

    warmup_zgemm(rank, nodes, random_seed, parsec);

    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    /* initializing matrix structure */
    if(!check)
    {
        PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
            parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDA, K, 0, 0,
                                   M, K, P, nodes/P, KP, KQ, IP, JQ));
        PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
            parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, N, 0, 0,
                                   K, N, P, nodes/P, KP, KQ, IP, JQ));

        /* matrix generation */
        if(loud > 2) printf("+++ Generate matrices ... ");
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, Aseed);
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, Bseed);
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed);
        if(loud > 2) printf("Done\n");

        int t;
        for(t = 0; t < nruns; t++) {
            parsec_devices_release_memory();
            /* Create PaRSEC */
            PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec, zgemm,
                                      (tA, tB, alpha,
                                       (parsec_tiled_matrix_t *)&dcA,
                                       (parsec_tiled_matrix_t *)&dcB,
                                       beta,
                                       (parsec_tiled_matrix_t *)&dcC),
                                      dplasma_zgemm_Destruct( PARSEC_zgemm ));

            parsec_devices_reset_load(parsec);
        }

        parsec_data_free(dcA.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
        parsec_data_free(dcB.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
    }
    else {
        int Am, An, Bm, Bn;

#if defined(PRECISION_z) || defined(PRECISION_c)
        for(tA=0; tA<3; tA++) {
            for(tB=0; tB<3; tB++) {
#else
        for(tA=0; tA<2; tA++) {
            for(tB=0; tB<2; tB++) {
#endif
                if ( trans[tA] == dplasmaNoTrans ) {
                    Am = M; An = K;
                } else {
                    Am = K; An = M;
                }
                if ( trans[tB] == dplasmaNoTrans ) {
                    Bm = K; Bn = N;
                } else {
                    Bm = N; Bn = K;
                }
                PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
                    parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                           rank, MB, NB, LDA, LDA, 0, 0,
                                           Am, An, P, nodes/P, KP, KQ, IP, JQ));
                PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
                    parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                           rank, MB, NB, LDB, LDB, 0, 0,
                                           Bm, Bn, P, nodes/P, KP, KQ, IP, JQ));

                if ( rank == 0 ) {
                    printf("***************************************************\n");
                    printf(" ----- TESTING ZGEMM (%s, %s) -------- \n",
                           transstr[tA], transstr[tB]);
                }

                /* matrix generation */
                if(loud) printf("Generate matrices ... ");
                dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, Aseed);
                dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, Bseed);
                dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed);
                if(loud) printf("Done\n");

                /* reset device load incurred by prior rounds, if any */
                parsec_devices_release_memory();
                parsec_devices_reset_load(parsec);

                /* Create GEMM PaRSEC */
                if(loud) printf("Compute ... ... ");
                    PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec, zgemm,
                              (trans[tA], trans[tB],
                               (dplasma_complex64_t)alpha,
                               (parsec_tiled_matrix_t *)&dcA,
                               (parsec_tiled_matrix_t *)&dcB,
                               (dplasma_complex64_t)beta,
                               (parsec_tiled_matrix_t *)&dcC),
                              dplasma_zgemm_Destruct( PARSEC_zgemm ));

                if(loud) printf("Done\n");

                parsec_data_free(dcA.mat);
                parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
                parsec_data_free(dcB.mat);
                parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);

                /* Check the solution */
                info_solution = check_solution( parsec, (rank == 0) ? loud : 0,
                                                trans[tA], trans[tB],
                                                alpha, Am, An, Aseed,
                                                       Bm, Bn, Bseed,
                                                beta,  M,  N,  Cseed,
                                                &dcC);
                if ( rank == 0 ) {
                    if (info_solution == 0) {
                        printf(" ---- TESTING ZGEMM (%s, %s) ...... PASSED !\n",
                               transstr[tA], transstr[tB]);
                    }
                    else {
                        printf(" ---- TESTING ZGEMM (%s, %s) ... FAILED !\n",
                               transstr[tA], transstr[tB]);
                    }
                    printf("***************************************************\n");
                }
            }
        }
    }

    parsec_data_free(dcC.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC);

    cleanup_parsec(parsec, iparam);

    return info_solution;
}

/**********************************
 * static functions
 **********************************/

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
static int check_solution( parsec_context_t *parsec, int loud,
                           dplasma_enum_t transA, dplasma_enum_t transB,
                           dplasma_complex64_t alpha, int Am, int An, int Aseed,
                                                      int Bm, int Bn, int Bseed,
                           dplasma_complex64_t beta,  int M,  int N,  int Cseed,
                           parsec_matrix_block_cyclic_t *dcCfinal )
{
    int info_solution = 1;
    double Anorm, Bnorm, Cinitnorm, Cdplasmanorm, Clapacknorm, Rnorm;
    double eps, result;
    int K  = ( transA == dplasmaNoTrans ) ? An : Am ;
    int MB = dcCfinal->super.mb;
    int NB = dcCfinal->super.nb;
    int LDA = Am;
    int LDB = Bm;
    int LDC = M;
    int rank  = dcCfinal->super.super.myrank;

    eps = LAPACKE_dlamch_work('e');

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDA, An, 0, 0,
                               Am, An, 1, 1, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
        parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDB, Bn, 0, 0,
                               Bm, Bn, 1, 1, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, 1, 1, 1, 1, 0, 0));

    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, Aseed );
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, Bseed );
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC, Cseed );

    Anorm        = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcA );
    Bnorm        = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcB );
    Cinitnorm    = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcC );
    Cdplasmanorm = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)dcCfinal );

    if ( rank == 0 ) {
        cblas_zgemm(CblasColMajor,
                    (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                    M, N, K,
                    CBLAS_SADDR(alpha), dcA.mat, LDA,
                                        dcB.mat, LDB,
                    CBLAS_SADDR(beta),  dcC.mat, LDC );
    }

    Clapacknorm = dplasma_zlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&dcC );

    dplasma_zgeadd( parsec, dplasmaNoTrans, -1.0, (parsec_tiled_matrix_t*)dcCfinal,
                                           1.0, (parsec_tiled_matrix_t*)&dcC );

    Rnorm = dplasma_zlange( parsec, dplasmaMaxNorm, (parsec_tiled_matrix_t*)&dcC);

    if ( rank == 0 ) {
        if ( loud > 2 ) {
            printf("  ||A||_inf = %e, ||B||_inf = %e, ||C||_inf = %e\n"
                   "  ||lapack(a*A*B+b*C)||_inf = %e, ||dplasma(a*A*B+b*C)||_inf = %e, ||R||_m = %e\n",
                   Anorm, Bnorm, Cinitnorm, Clapacknorm, Cdplasmanorm, Rnorm);
        }

        result = Rnorm / ((Anorm + Bnorm + Cinitnorm) * max(M,N) * eps);
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
