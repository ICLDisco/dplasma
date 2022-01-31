/*
 * Copyright (c) 2009-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    double *work = NULL;
    double result;
    double normlap = 0.0;
    double normdag = 0.0;
    double eps = LAPACKE_dlamch_work('e');
    int iparam[IPARAM_SIZEOF];
    int An, i, u, ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 40, 200, 200);
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'m';
    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);

    An = dplasma_imax(M, N);
    LDA = max( LDA, M );

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA0, 1,
        parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB, NB, LDA, An, 0, 0,
                               M, An, 1, 1, KP, KQ, IP, JQ));

    if( rank == 0 ) {
        work = (double *)malloc( max(M,N) * sizeof(double));
    }

    /*
     * General cases LANGE
     */
    {
        PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
            parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   M, N, P, nodes/P, KP, KQ, IP, JQ));

        /* matrix generation */
        if(loud > 2) printf("+++ Generate matrices ... ");
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA0, 3872);
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA,  3872);
        if(loud > 2) printf("Done\n");

        for(i=0; i<4; i++) {
            if ( rank == 0 ) {
                printf("***************************************************\n");
            }
            if(loud > 2) printf("+++ Computing norm %s ... ", normsstr[i]);
            normdag = dplasma_zlange(parsec, norms[i],
                                     (parsec_tiled_matrix_t *)&dcA);

            if ( rank == 0 ) {
                normlap = LAPACKE_zlange_work(LAPACK_COL_MAJOR, normsstr[i][0], M, N,
                                              (dplasma_complex64_t*)(dcA0.mat), dcA0.super.lm, work);
            }
            if(loud > 2) printf("Done.\n");

            if ( loud > 3 ) {
                printf( "%d: The norm %s of A is %e\n",
                        rank, normsstr[i], normdag);
            }

            if ( rank == 0 ) {
                result = fabs(normdag - normlap) / (normlap * eps) ;

                if ( loud > 3 ) {
                    printf( "%d: The norm %s of A is %e (LAPACK)\n",
                            rank, normsstr[i], normlap);
                }

                switch(norms[i]) {
                case dplasmaMaxNorm:
                    /* result should be perfectly equal */
                    break;
                case dplasmaInfNorm:
                    /* Sum order on the line can differ */
                    result = result / (double)N;
                    break;
                case dplasmaOneNorm:
                    /* Sum order on the column can differ */
                    result = result / (double)M;
                    break;
                case dplasmaFrobeniusNorm:
                    /* Sum order on every element can differ */
                    result = result / ((double)M * (double)N);
                    break;
                }

                if ( result < 1. ) {
                    printf(" ----- TESTING ZLANGE (%s) ... SUCCESS !\n", normsstr[i]);
                } else {
                    printf("       Ndag = %e, Nlap = %e\n", normdag, normlap );
                    printf("       | Ndag - Nlap | / Nlap = %e\n", result);
                    printf(" ----- TESTING ZLANGE (%s) ... FAILED !\n", normsstr[i]);
                    ret |= 1;
                }
            }
        }

        parsec_data_free(dcA.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    }

    /*
     * Triangular cases LANTR
     */
    {
        int d;

        /* matrix generation */
        if(loud > 2) printf("+++ Generate matrices ... ");
        dplasma_zplrnt( parsec, 0., (parsec_tiled_matrix_t *)&dcA0, 3872);
        if(loud > 2) printf("Done\n");

        /* Computing the norm */
        PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
            parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   M, N, P, nodes/P, KP, KQ, IP, JQ));

        for(u=0; u<2; u++) {
            dplasma_zplrnt( parsec, 0., (parsec_tiled_matrix_t *)&dcA, 3872);

            for(d=0; d<2; d++) {
                for(i=0; i<4; i++) {
                    if ( rank == 0 ) {
                        printf("***************************************************\n");
                    }
                    if(loud > 2) printf("+++ Computing norm %s ... ", normsstr[i]);
                    normdag = dplasma_zlantr(parsec, norms[i], uplo[u], diag[d],
                                             (parsec_tiled_matrix_t *)&dcA);

                    if ( rank == 0 ) {
                        normlap = LAPACKE_zlantr_work(LAPACK_COL_MAJOR, normsstr[i][0], uplostr[u][0], diagstr[d][0], M, N,
                                                      (dplasma_complex64_t*)(dcA0.mat), dcA0.super.lm, work);
                    }
                    if(loud > 2) printf("Done.\n");

                    if ( loud > 3 ) {
                        printf( "%d: The norm %s of A is %e\n",
                                rank, normsstr[i], normdag);
                    }

                    if ( rank == 0 ) {
                        result = fabs(normdag - normlap) / (normlap * eps);

                        if ( loud > 3 ) {
                            printf( "%d: The norm %s of A is %e (LAPACK)\n",
                                    rank, normsstr[i], normlap);
                        }

                        switch(norms[i]) {
                        case dplasmaMaxNorm:
                            /* result should be perfectly equal */
                            break;
                        case dplasmaInfNorm:
                            /* Sum order on the line can differ */
                            result = result / (double)N;
                            break;
                        case dplasmaOneNorm:
                            /* Sum order on the column can differ */
                            result = result / (double)M;
                            break;
                        case dplasmaFrobeniusNorm:
                            /* Sum oreder on every element can differ */
                            result = result / ((double)M * (double)N);
                            break;
                        }

                        if ( result < 1. ) {
                            printf(" ----- TESTING ZLANTR (%s, %s, %s) ... SUCCESS !\n",
                                   normsstr[i], uplostr[u], diagstr[d]);
                        } else {
                            printf(" ----- TESTING ZLANTR (%s, %s, %s) ... FAILED !\n",
                                   normsstr[i], uplostr[u], diagstr[d]);
                            printf("       | Ndag - Nlap | / Nlap = %e\n", normdag);
                            ret |= 1;
                        }
                    }
                }
            }
        }
        parsec_data_free(dcA.mat);
        parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&dcA);
    }

    /* Let set N=M for the triangular cases */
    N = M;

    /*
     * Symmetric cases LANSY
     */
    {
        /* matrix generation */
        if(loud > 2) printf("+++ Generate matrices ... ");
        dplasma_zplgsy( parsec, 0., dplasmaUpperLower, (parsec_tiled_matrix_t *)&dcA0, 3872);
        if(loud > 2) printf("Done\n");

        for(u=0; u<2; u++) {

            /* Computing the norm */
            PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
                parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                           rank, MB, NB, LDA, N, 0, 0,
                                           M, N, P, nodes/P, uplo[u]));

            dplasma_zplgsy( parsec, 0., uplo[u], (parsec_tiled_matrix_t *)&dcA, 3872);

            for(i=0; i<4; i++) {
                if ( rank == 0 ) {
                    printf("***************************************************\n");
                }
                if(loud > 2) printf("+++ Computing norm %s ... ", normsstr[i]);
                normdag = dplasma_zlansy(parsec, norms[i], uplo[u],
                                         (parsec_tiled_matrix_t *)&dcA);

                if ( rank == 0 ) {
                    normlap = LAPACKE_zlansy_work(LAPACK_COL_MAJOR, normsstr[i][0], uplostr[u][0], M,
                                                  (dplasma_complex64_t*)(dcA0.mat), dcA0.super.lm, work);
                }
                if(loud > 2) printf("Done.\n");

                if ( loud > 3 ) {
                    printf( "%d: The norm %s of A is %e\n",
                            rank, normsstr[i], normdag);
                }

                if ( rank == 0 ) {
                    result = fabs(normdag - normlap) / (normlap * eps);

                    if ( loud > 3 ) {
                        printf( "%d: The norm %s of A is %e (LAPACK)\n",
                                rank, normsstr[i], normlap);
                    }

                    switch(norms[i]) {
                    case dplasmaMaxNorm:
                        /* result should be perfectly equal */
                        break;
                    case dplasmaInfNorm:
                        /* Sum order on the line can differ */
                        result = result / (double)N;
                        break;
                    case dplasmaOneNorm:
                        /* Sum order on the column can differ */
                        result = result / (double)M;
                        break;
                    case dplasmaFrobeniusNorm:
                        /* Sum oreder on every element can differ */
                        result = result / ((double)M * (double)N);
                        break;
                    }

                    if ( result < 1. ) {
                        printf(" ----- TESTING ZLANSY (%s, %s) ... SUCCESS !\n", uplostr[u], normsstr[i]);
                    } else {
                        printf(" ----- TESTING ZLANSY (%s, %s) ... FAILED !\n", uplostr[u], normsstr[i]);
                        printf("       | Ndag - Nlap | / Nlap = %e\n", normdag);
                        ret |= 1;
                    }
                }
            }

            parsec_data_free(dcA.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
        }
    }

#if defined(PRECISION_z) || defined(PRECISION_c)
    /*
     * Hermitian cases LANHE
     */
    {
        /* matrix generation */
        if(loud > 2) printf("+++ Generate matrices ... ");
        dplasma_zplghe( parsec, 0., dplasmaUpperLower, (parsec_tiled_matrix_t *)&dcA0, 3872);
        if(loud > 2) printf("Done\n");

        for(u=0; u<2; u++) {

            /* Computing the norm */
            PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
                parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                           rank, MB, NB, LDA, N, 0, 0,
                                           M, N, P, nodes/P, uplo[u]));

            dplasma_zplghe( parsec, 0., uplo[u], (parsec_tiled_matrix_t *)&dcA, 3872);

            for(i=0; i<4; i++) {
                if ( rank == 0 ) {
                    printf("***************************************************\n");
                }
                if(loud > 2) printf("+++ Computing norm %s ... ", normsstr[i]);
                normdag = dplasma_zlanhe(parsec, norms[i], uplo[u],
                                         (parsec_tiled_matrix_t *)&dcA);

                if ( rank == 0 ) {
                    normlap = LAPACKE_zlanhe_work(LAPACK_COL_MAJOR, normsstr[i][0], uplostr[u][0], M,
                                                  (dplasma_complex64_t*)(dcA0.mat), dcA0.super.lm, work);
                }
                if(loud > 2) printf("Done.\n");

                if ( loud > 3 ) {
                    printf( "%d: The norm %s of A is %e\n",
                            rank, normsstr[i], normdag);
                }

                if ( rank == 0 ) {
                    result = fabs(normdag - normlap) / (normlap * eps);

                    if ( loud > 3 ) {
                        printf( "%d: The norm %s of A is %e (LAPACK)\n",
                                rank, normsstr[i], normlap);
                    }
                    switch(norms[i]) {
                    case dplasmaMaxNorm:
                        /* result should be perfectly equal */
                        break;
                    case dplasmaInfNorm:
                        /* Sum order on the line can differ */
                        result = result / (double)N;
                        break;
                    case dplasmaOneNorm:
                        /* Sum order on the column can differ */
                        result = result / (double)M;
                        break;
                    case dplasmaFrobeniusNorm:
                        /* Sum oreder on every element can differ */
                        result = result / ((double)M * (double)N);
                        break;
                    }

                    if ( result < 1. ) {
                        printf(" ----- TESTING ZLANHE (%s, %s) ... SUCCESS !\n", uplostr[u], normsstr[i]);
                    } else {
                        printf(" ----- TESTING ZLANHE (%s, %s) ... FAILED !\n", uplostr[u], normsstr[i]);
                        printf("       | Ndag - Nlap | / Nlap = %e\n", normdag);
                        ret |= 1;
                    }
                }
            }

            parsec_data_free(dcA.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
        }
    }
#endif

    if ( rank == 0 ) {
        printf("***************************************************\n");
        free( work );
    }
    parsec_data_free(dcA0.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0);

    cleanup_parsec(parsec, iparam);

    return ret;
}
