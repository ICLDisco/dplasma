/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2015-2016 Inria, CNRS (LaBRI - UMR 5800), University of
 *                         Bordeaux and Bordeaux INP. All rights reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

int TS_cp1QR( int p ) {
    if ( p > 0 )
        return 4 + 6 + 12 * (p-1);
    else
        return 0;
}

int TT_cp1QR( int p ) {
    if ( p > 0 )
        return 4 + 6 + 6 * (p-1);
    else
        return 0;
}

int GD_cp1QR( int p ) {
    if ( p > 0 )
        return 4 + 6 + 6 * ceil( log2( p ) );
    else
        return 0;
}

int TS_cpBidiag( int p, int q ) {
    if ( (p > 0) && (q > 0) )
        return 12 * p * q - 6 * p +  2 * q -  4;
    else
        return 0;
}

int TT_cpBidiag( int p, int q ) {
    if ( (p > 0) && (q > 0) )
        return 6 * p * q -  4 * p + 12 * q - 10;
    else
        return 0;
}

int GD_cpBidiag( int p, int q ) {
    int k, cp;
    (void)p; (void)q;
    cp = 4 + 2 * ceil( log2( p+1-q ) );
    for (k=1; k<q; k++) {
        cp += GD_cp1QR( p + 1 - k );
    }
    for (k=2; k<q+1; k++) {
        cp += GD_cp1QR( q + 1 - k );
    }

    return cp;
}

int TS_cpQR( int p, int q ) {
    if (q == 1) {
        return 4 + (p-q) * 6;
    }
    else {
        if (p == q) {
            return 30 * (q-1) - 4;
        } else {
            return 30 * (q-1) - 2 + (p-q) * 12;
        }
    }
}

int TT_cpQR( int p, int q ) {
    if (q == 1) {
        return 2 * p + 2;
    }
    else {
        if (p == q) {
            return 22 * q - 24;
        } else {
            return 6 * p + 16 * q - 22;
        }
    }
}

int GD_cpQR( int p, int q ) {
    if (q == 1) {
        return  4 + 2 * ceil( log2(p) );
    } else {
        return 22 * q + 6 * ceil( log2(p) );
    }
}

int RunOneTest( parsec_context_t *parsec, int nodes, int cores, int rank, int loud,
                int M, int N, int LDA, int MB, int NB, int IB, int P, int Q, int hmb,
                int ltre0, int htre0, int ltree, int htree, int ts, int domino, int rbidiag,
                int nruns)
{
    int ret = 0;
    dplasma_qrtree_t qrtre0, qrtree, lqtree;
    parsec_tiled_matrix_t *subA = NULL;
    int minMN;
    int MT = (M%MB==0) ? (M/MB) : (M/MB+1);
    int NT = (N%NB==0) ? (N/NB) : (N/NB+1);
    int cp = -1;
    int rc;

    //PASTE_CODE_FLOPS(FLOPS_ZGEBRD, ((DagDouble_t)M, (DagDouble_t)N));
    DagDouble_t flops, gflops = -1.;
    DagDouble_t gflops_min = 9999999999999.;
    DagDouble_t gflops_max = 0.;
    DagDouble_t gflops_avg = 0.;
    DagDouble_t gflops_sum = 0.;
    DagDouble_t time_avg = 0.;
    flops = FLOPS_ZGEQRF( (DagDouble_t)M, (DagDouble_t)N )
        +   FLOPS_ZGELQF( (DagDouble_t)M, (DagDouble_t)(N-NB) );
    LDA = max(M, LDA);

    if ( M < N ) {
        fprintf(stderr, "This testing can only perform SVD on matrices with M >= N\n");
        return -1;
    }
    minMN = dplasma_imin(M, N);

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcTS0, rbidiag,
        parsec_matrix_block_cyclic, (&dcTS0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcTT0, rbidiag,
        parsec_matrix_block_cyclic, (&dcTT0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcTS, 1,
        parsec_matrix_block_cyclic, (&dcTS, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcTT, 1,
        parsec_matrix_block_cyclic, (&dcTT, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcBand, 1,
        parsec_matrix_block_cyclic, (&dcBand, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB+1, NB, MB+1, minMN, 0, 0,
                               MB+1, minMN, 1, 1, 1, 1, 0, 0));

    /* Initialize the matrix */
    if(loud > 3) printf("+++ Generate matrices ... ");

    if ( rbidiag ) {
        subA = parsec_tiled_matrix_submatrix( (parsec_tiled_matrix_t *)&dcA,
                                       0, 0, dcA.super.n, dcA.super.n );

        dplasma_hqr_init( &qrtre0,
                          dplasmaNoTrans, (parsec_tiled_matrix_t *)&dcA,
                          ltre0, htre0, ts, P, domino, 0 );

        /**
         * Adaptativ tree
         */
        if (ltree == 9) {
            dplasma_svd_init( &qrtree,
                              dplasmaNoTrans, subA,
                              htree, P, cores, hmb );

            dplasma_svd_init( &lqtree,
                              dplasmaTrans, subA,
                              htree, Q, cores, hmb );
        } else {
#if defined(PARSEC_SIM)
            if (ltree == DPLASMA_FLAT_TREE) {
                if (ts == 1) {
                    cp = TT_cpQR( MT, NT );
                    if (NT == 1) {
                    }
                    else if (NT == 2) {
                        cp += 6 /* LQ(1) : Update on last diagonal block */
                            + 4 /* QR(2) : Facto on last diagonal block */;
                    }
                    else if (NT == 3) {
                        cp += 12                         /* First LQ(1) = UNMQR + TTMLQ on last block */
                            + TT_cpBidiag( NT-1, NT-1 ); /* Bidiag of the remaining matrix */
                    }
                    else {
                        cp = 6 * NT * NT + 6 * NT + 6 * MT - 16;
                    }
                }
                else {
                    cp = TS_cpQR( MT, NT );
                    if (NT == 1) {
                    }
                    else if (NT == 2) {
                        cp += 6 /* LQ(1) : Update on last diagonal block */
                            + 4 /* QR(2) : Facto on last diagonal block */;
                    }
                    else {
                        if ( NT == 3 && NT == MT ) {
                            cp = 12 * NT * NT - 16 * NT + 12 * MT + 4;
                        }
                        else {
                            cp = 12 * NT * NT - 16 * NT + 12 * MT + 6;
                        }
                    }
                }
            }
            else if (ltree == DPLASMA_GREEDY1P_TREE) {
                if (ts == 1) {
                    cp = GD_cpQR( MT, NT );
                    if (NT == 1) {
                    }
                    else if (NT == 2) {
                        cp += 6 /* LQ(1) : Update on last diagonal block */
                            + 4 /* QR(2) : Facto on last diagonal block */;
                    }
                    else if (NT == 3) {
                        cp += 12                         /* First LQ(1) = UNMQR + TTMLQ on last block */
                            + GD_cpBidiag( NT-1, NT-1 ); /* Bidiag of the remaining matrix */
                    }
                    else {
                        cp += 12                            /* First LQ(1) = UNMQR + TTMLQ on last block */
                            + NT * 6 + ((NT == MT) ? 2 : 0) /* First QR(2) = unmqr( 2, q ) + (p-1) * 12 = ( 6 - (q-2) * 6 ) + (q-1) * 12 */
                            + GD_cpBidiag( NT-1, NT-1 ) - GD_cp1QR( NT-1 ); /* Bidiag of the remaining matrix starting by LQ */
                    }
                }
            }
#endif /* defined(PARSEC_SIM) */
            dplasma_hqr_init( &qrtree,
                              dplasmaNoTrans, subA,
                              ltree, htree, ts, P, 0, 0 );

            dplasma_hqr_init( &lqtree,
                              dplasmaTrans, subA,
                              ltree, htree, ts, Q, 0, 0 );
        }
    }
    else {
        /**
         * Adaptativ tree
         */
        if (ltree == 9) {
            dplasma_svd_init( &qrtree,
                              dplasmaNoTrans, (parsec_tiled_matrix_t *)&dcA,
                              htree, P, cores, hmb );

            dplasma_svd_init( &lqtree,
                              dplasmaTrans, (parsec_tiled_matrix_t *)&dcA,
                              htree, Q, cores, hmb );
        } else {
#if defined(PARSEC_SIM)
            if (ltree == 0) {
                if (ts == 1) {
                    cp = TT_cpBidiag( MT, NT );
                }
                else {
                    cp = TS_cpBidiag( MT, NT );
                }
            }
            else if (ltree == 4) {
                if (ts == 1) {
                    cp = GD_cpBidiag( MT, NT );
                }
            }
#endif /* defined(PARSEC_SIM) */
            dplasma_hqr_init( &qrtree,
                              dplasmaNoTrans, (parsec_tiled_matrix_t *)&dcA,
                              ltree, htree, ts, P, 0, 0 );

            dplasma_hqr_init( &lqtree,
                              dplasmaTrans, (parsec_tiled_matrix_t *)&dcA,
                              ltree, htree, ts, Q, 0, 0 );
        }
    }

#if defined(PARSEC_SIM)
    nruns = 1;
#endif

    for (int t=0; t<nruns+1; t++) {

        if(loud > 3) printf("+++ Generate matrices ... ");
        /* Generate the matrix on rank 0 */
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, 3872);
        if(loud > 3) printf("Done\n");

        /* Create Parsec */
        PASTE_CODE_ENQUEUE_KERNEL(parsec, zgebrd_ge2gbx,
                                  (IB, rbidiag ? &qrtre0 : &qrtree,
                                   &qrtree, &lqtree,
                                   (parsec_tiled_matrix_t*)&dcA,
                                   rbidiag ? (parsec_tiled_matrix_t*)&dcTS0 : (parsec_tiled_matrix_t*)&dcTS,
                                   rbidiag ? (parsec_tiled_matrix_t*)&dcTT0 : (parsec_tiled_matrix_t*)&dcTT,
                                   (parsec_tiled_matrix_t*)&dcTS,
                                   (parsec_tiled_matrix_t*)&dcTT,
                                   (parsec_tiled_matrix_t*)&dcBand));

        /* lets rock! */
        SYNC_TIME_START();
        rc = parsec_context_start(parsec);
        PARSEC_CHECK_ERROR(rc, "parsec_context_start");
        TIME_START();
        rc = parsec_context_wait(parsec);
        PARSEC_CHECK_ERROR(rc, "parsec_context_wait");
        SYNC_TIME_STOP();

        dplasma_zgebrd_ge2gbx_Destruct( PARSEC_zgebrd_ge2gbx );

        time_avg += sync_time_elapsed;
        gflops = (flops/1.e9)/(sync_time_elapsed);

        if(t > 0) {
            if (rank == 0){
                fprintf(stdout,
                        "zgebrd_ge2gb M= %2d N= %2d NP= %2d NC= %2d P= %2d Q= %2d NB= %2d IB= %2d R-bidiag= %2d treeh= %2d treel_rb= %2d qr_a= %2d QR(domino= %2d treel_qr= %2d ) : %.2f s %f gflops\n",
                        M, N, nodes, cores, P, Q, NB, IB,
                        rbidiag, htree, ltree, ts, domino, ltre0,
                        sync_time_elapsed, gflops);
            }
            gflops_min = (gflops_min > gflops) ? gflops : gflops_min;
            gflops_max = (gflops_max < gflops) ? gflops : gflops_max;
            gflops_avg += gflops;
            gflops_sum += (gflops * gflops);
        }
    }
    time_avg   = time_avg   / (double)nruns;
    gflops_avg = gflops_avg / (double)nruns;
    gflops_sum = sqrt( (gflops_sum/(double)nruns) - (gflops_avg * gflops_avg) );

#if defined(PARSEC_SIM)
    {
        (void)flops; (void)gflops;
        if ( rank == 0 ) {
            printf("zgebrd_ge2gb simulation M= %2d N= %2d NP= %2d NC= %2d P= %2d Q= %2d NB= %2d IB= %2d R-bidiag= %2d treeh= %2d treel_rb= %2d qr_a= %2d QR(domino= %2d treel_qr= %2d) : %3d - %3d ( error= %4d )\n",
                   M, N, nodes, cores, P, Q, NB, IB,
                   rbidiag, htree, ltree, ts, domino, ltre0,
                   parsec_getsimulationdate( parsec ), cp,
                   parsec_getsimulationdate( parsec ) - cp );
        }
    }
#else
    SYNC_TIME_PRINT(rank,
                    ("zgebrd_ge2gb computation M= %2d N= %2d NP= %2d NC= %2d P= %2d Q= %2d NB= %2d IB= %2d R-bidiag= %2d treeh= %2d treel_rb= %2d qr_a= %2d QR(domino= %2d treel_qr= %2d ) : %.2f s %f gflops (%f, %f, %f)\n",
                     M, N, nodes, cores, P, Q, NB, IB,
                     rbidiag, htree, ltree, ts, domino, ltre0,
                     time_avg, gflops_avg, gflops_min, gflops_max, gflops_sum));
#endif
    fflush(stdout);
    PASTE_CODE_PERF_LOOP_DONE();

    dplasma_hqr_finalize( &qrtree );
    dplasma_hqr_finalize( &lqtree );

    if (rbidiag) {
        dplasma_hqr_finalize( &qrtre0 );
        free(subA);
        parsec_data_free(dcTS0.mat);
        parsec_data_free(dcTT0.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcTS0);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcTT0);
    }

    parsec_data_free(dcA.mat);
    parsec_data_free(dcTS.mat);
    parsec_data_free(dcTT.mat);
    parsec_data_free(dcBand.mat);

    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcTS);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcTT);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcBand);

    (void)cp; (void)NT;

    return ret;
}


int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int m, ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 32, 200, 200);
    iparam[IPARAM_KP] = 1;
    iparam[IPARAM_KQ] = 1;
    iparam[IPARAM_HMB] = 2;
    iparam[IPARAM_HNB] = 1;
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'m';

    /* Initialize Parsec */
    parsec = setup_parsec(argc, argv, iparam);

    PASTE_CODE_IPARAM_LOCALS(iparam);

    int ltree = iparam[IPARAM_LOWLVL_TREE] == DPLASMA_GREEDY_TREE ? DPLASMA_GREEDY1P_TREE : iparam[IPARAM_LOWLVL_TREE];
    ltree = iparam[IPARAM_ASYNC] ? ltree : 9;

    /**
     * Test for varying matrix sizes m-by-n where:
     *    1) m = M .. N .. K, and n = m  (square)
     *    2) m = N .. M .. K, and n = N  (square to TS)
     *
     * -r --tsrr  : enable/disable the r-bidiagonalization
     * -b --sync  : enable/disbale the auto-adaptativ tree for the (R-)bidiagonalization step
     *
     *    --treel : Tree used for low level reduction inside nodes. (specific to xgeqrf_param)\n"
     *    --treeh : Tree used for high level reduction between nodes, only if qr_p > 1. (specific to xgeqrf_param)\n"
     *              (0: Flat, 1: Greedy, 2: Fibonacci, 3: Binary)\n"
     *    --qr_a  : Size of TS domain. (specific to xgeqrf_param)\n"
     *    --qr_p  : Size of the high level tree. (specific to xgeqrf_param)\n"
     * -d --domino: Enable/Disable the domino between upper and lower trees. (specific to xgeqrf_param) (default: 1)\n"
     */
    for (m=M; m<=N; m+=K ) {
        RunOneTest( parsec, nodes, iparam[IPARAM_NCORES], rank, loud,
                    m, m, LDA, MB, NB, IB, P, Q, iparam[IPARAM_HMB],
                    iparam[IPARAM_LOWLVL_TREE], iparam[IPARAM_HIGHLVL_TREE],
                    ltree, iparam[IPARAM_HIGHLVL_TREE],
                    iparam[IPARAM_QR_TS_SZE], iparam[IPARAM_QR_DOMINO], iparam[IPARAM_QR_TSRR],
                    iparam[IPARAM_NRUNS]);
    }

    for (m=N; m<=M; m+=K ) {
        RunOneTest( parsec, nodes, iparam[IPARAM_NCORES], rank, loud,
                    m, N, LDA, MB, NB, IB, P, Q, iparam[IPARAM_HMB],
                    iparam[IPARAM_LOWLVL_TREE], iparam[IPARAM_HIGHLVL_TREE],
                    ltree, iparam[IPARAM_HIGHLVL_TREE],
                    iparam[IPARAM_QR_TS_SZE], iparam[IPARAM_QR_DOMINO], iparam[IPARAM_QR_TSRR],
                    iparam[IPARAM_NRUNS] );
    }

    cleanup_parsec(parsec, iparam);

    return ret;
}

