/*
 * Copyright (c) 2009-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    double Cnorm, Rnorm, result;
    double threshold = 10.;
    double eps = LAPACKE_dlamch_work('e');
    char *resultstr;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;
    int s, t;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 48, 192, 192);
    iparam[IPARAM_KP] = 4;
    iparam[IPARAM_KQ] = 1;
    iparam[IPARAM_LDA] = -'m';

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);

    if (M < K) {
        printf("WARNING: M must be greater or equal to K (Set M = K)\n");
        M = K;
        MT = KT;
    }

    LDA = max(M, LDA);

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, K, 0, 0,
                               M, K, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcT, 1,
        parsec_matrix_block_cyclic, (&dcT, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, K, 0, 0,
                               MT*IB, K, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcQ, 1,
        parsec_matrix_block_cyclic, (&dcQ, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, M, 0, 0,
                               M, M, P, nodes/P, KP, KQ, IP, JQ));

    /* matrix generation */
    if(loud > 3) printf("+++ Generate matrices ... ");
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, 3872);
    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0., (parsec_tiled_matrix_t *)&dcT);
    if(loud > 3) printf("Done\n");

    if(loud > 3) printf("+++ Factorize A ... ");
    dplasma_zgeqrf(parsec,
                   (parsec_tiled_matrix_t*)&dcA,
                   (parsec_tiled_matrix_t*)&dcT);
    if(loud > 3) printf("Done\n");

    if(loud > 3) printf("+++ Generate Q ... ");
    dplasma_zungqr( parsec,
                    (parsec_tiled_matrix_t *)&dcA,
                    (parsec_tiled_matrix_t *)&dcT,
                    (parsec_tiled_matrix_t *)&dcQ);
    if(loud > 3) printf("Done\n");

    for (s=0; s<2; s++) {

        int Cm = (sides[s] == dplasmaLeft) ? M : N;
        int Cn = (sides[s] == dplasmaLeft) ? N : M;
        LDC = max(LDC, Cm);

        PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
            parsec_matrix_block_cyclic, (&dcC, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDC, Cn, 0, 0,
                                   Cm, Cn, P, nodes/P, KP, KQ, IP, JQ));
        PASTE_CODE_ALLOCATE_MATRIX(dcC0, 1,
            parsec_matrix_block_cyclic, (&dcC0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDC, Cn, 0, 0,
                                   Cm, Cn, P, nodes/P, KP, KQ, IP, JQ));

        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcC0, 2354);
        Cnorm = dplasma_zlange(parsec, dplasmaOneNorm, (parsec_tiled_matrix_t *)&dcC0);

        if (Cnorm == 0.)
            Cnorm = 1.;

        for (t=0; t<2; t++) {
#if defined(PRECISION_z) || defined(PRECISION_c)
            if (t==1) t++;
#endif

            dplasma_zlacpy( parsec, dplasmaUpperLower,
                            (parsec_tiled_matrix_t *)&dcC0,
                            (parsec_tiled_matrix_t *)&dcC);

            dplasma_zunmqr( parsec, sides[s], trans[t],
                            (parsec_tiled_matrix_t *)&dcA,
                            (parsec_tiled_matrix_t *)&dcT,
                            (parsec_tiled_matrix_t *)&dcC);

            if (sides[s] == dplasmaLeft ) {
                dplasma_zgemm( parsec, trans[t], dplasmaNoTrans,
                               -1., (parsec_tiled_matrix_t *)&dcQ,
                                    (parsec_tiled_matrix_t *)&dcC0,
                               1.,  (parsec_tiled_matrix_t *)&dcC);
            } else {
                dplasma_zgemm( parsec, dplasmaNoTrans, trans[t],
                               -1., (parsec_tiled_matrix_t *)&dcC0,
                                    (parsec_tiled_matrix_t *)&dcQ,
                               1.,  (parsec_tiled_matrix_t *)&dcC);
            }

            Rnorm = dplasma_zlange(parsec, dplasmaOneNorm, (parsec_tiled_matrix_t *)&dcC);
            result = Rnorm / ((double)M * Cnorm * eps);

            if (loud && rank == 0) {
                printf("***************************************************\n");
                if ( loud > 3 ) {
                    printf( "-- ||C||_1 = %e, ||R||_1 = %e, ||R||_1 / (M * ||C||_1 * eps) = %e\n",
                            Cnorm, Rnorm, result );
                }

                if (  isnan(Rnorm) || isinf(Rnorm) ||
                      isnan(result) || isinf(result) || (result >= threshold) ) {
                    resultstr = " FAILED";
                    ret |= 1;
                }
                else{
                    resultstr = "... PASSED";
                }
                printf(" ---- TESTING ZUNMQR (%s, %s) ...%s !\n",
                       sidestr[s], transstr[t], resultstr);
            }
        }

        parsec_data_free(dcC0.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC0);
        parsec_data_free(dcC.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcC);
    }

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    parsec_data_free(dcT.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcT);
    parsec_data_free(dcQ.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcQ);

    cleanup_parsec(parsec, iparam);

    return ret;
}
