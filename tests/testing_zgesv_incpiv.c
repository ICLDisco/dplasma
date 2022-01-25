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
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX );

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int info = 0;
    int info_solution = 0;
    int ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 40, 200, 200);
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'m';

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcA0, 1,
        parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
        parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDB, NRHS, 0, 0,
                               N, NRHS, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcX, 1,
        parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDB, NRHS, 0, 0,
                               N, NRHS, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcL, 1,
        parsec_matrix_block_cyclic, (&dcL, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcIPIV, 1,
        parsec_matrix_block_cyclic, (&dcIPIV, PARSEC_MATRIX_INTEGER, PARSEC_MATRIX_TILE,
                               rank, MB, 1, M, NT, 0, 0,
                               M, NT, P, nodes/P, KP, KQ, IP, JQ));

    /* matrix generation */
    if(loud > 2) printf("+++ Generate matrices ... ");
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA0, 3872);
    if(loud > 2) printf("Done\n");

    /*********************************************************************
     *               First Check
     */
    if ( rank == 0 ) {
        printf("***************************************************\n");
    }

    /* matrix generation */
    if( loud > 2 ) printf("Generate matrices ... ");
    dplasma_zlacpy( parsec, dplasmaUpperLower,
                    (parsec_tiled_matrix_t *)&dcA0, (parsec_tiled_matrix_t *)&dcA );
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, 3872);
    dplasma_zlacpy( parsec, dplasmaUpperLower,
                    (parsec_tiled_matrix_t *)&dcB, (parsec_tiled_matrix_t *)&dcX );
    if( loud > 2 ) printf("Done\n");

    /* Compute */
    if( loud > 2 ) printf("Compute ... ... ");
    info = dplasma_zgesv_incpiv(parsec,
                                (parsec_tiled_matrix_t *)&dcA,
                                (parsec_tiled_matrix_t *)&dcL,
                                (parsec_tiled_matrix_t *)&dcIPIV,
                                (parsec_tiled_matrix_t *)&dcX );
    if( loud > 2 ) printf("Done\n");
    if ( info != 0 ) printf("%d: Info = %d\n", rank, info);

    /* Check the solution */
    if ( info == 0 ) {
        info_solution = check_solution( parsec, rank ? 0 : loud,
                                        (parsec_tiled_matrix_t *)&dcA0,
                                        (parsec_tiled_matrix_t *)&dcB,
                                        (parsec_tiled_matrix_t *)&dcX );
    }
    if ( rank == 0 ) {
        if ( info || info_solution) {
            printf(" ----- TESTING ZGESV ... FAILED !\n");
            ret |= 1;
        }
        else {
            printf(" ----- TESTING ZGESV ....... PASSED !\n");
        }
        printf("***************************************************\n");
    }

    /*********************************************************************
     *               Second Check
     */
    if ( rank == 0 ) {
        printf("***************************************************\n");
    }
    
    /* matrix generation */
    if( loud > 2 ) printf("Generate matrices ... ");
    dplasma_zlacpy( parsec, dplasmaUpperLower,
                    (parsec_tiled_matrix_t *)&dcA0, (parsec_tiled_matrix_t *)&dcA );
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, 3872);
    dplasma_zlacpy( parsec, dplasmaUpperLower,
                    (parsec_tiled_matrix_t *)&dcB, (parsec_tiled_matrix_t *)&dcX );
    if( loud > 2 ) printf("Done\n");

    /* Compute */
    if( loud > 2 ) printf("Compute ... ... ");
    info = dplasma_zgetrf_incpiv(parsec,
                                 (parsec_tiled_matrix_t *)&dcA,
                                 (parsec_tiled_matrix_t *)&dcL,
                                 (parsec_tiled_matrix_t *)&dcIPIV );
    if ( info == 0 ) {
        dplasma_zgetrs_incpiv(parsec, dplasmaNoTrans,
                              (parsec_tiled_matrix_t *)&dcA,
                              (parsec_tiled_matrix_t *)&dcL,
                              (parsec_tiled_matrix_t *)&dcIPIV,
                              (parsec_tiled_matrix_t *)&dcX );
    }
    if ( loud > 2 ) printf("Done\n");
    if ( info != 0 ) printf("%d: Info = %d\n", rank, info);

    /* Check the solution */
    if ( info == 0 ) {
        info_solution = check_solution( parsec, rank ? 0 : loud,
                                        (parsec_tiled_matrix_t *)&dcA0,
                                        (parsec_tiled_matrix_t *)&dcB,
                                        (parsec_tiled_matrix_t *)&dcX );
    }
    if ( rank == 0 ) {
        if ( info || info_solution) {
            printf(" ----- TESTING ZGETRF + ZGETRS ... FAILED !\n");
            ret |= 1;
        }
        else {
            printf(" ----- TESTING ZGETRF + ZGETRS ....... PASSED !\n");
        }
        printf("***************************************************\n");
    }

    /*********************************************************************
     *               Third Check
     */
    if ( rank == 0 ) {
        printf("***************************************************\n");
    }
    
    /* matrix generation */
    if( loud > 2 ) printf("Generate matrices ... ");
    dplasma_zlacpy( parsec, dplasmaUpperLower,
                    (parsec_tiled_matrix_t *)&dcA0, (parsec_tiled_matrix_t *)&dcA );
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, 3872);
    dplasma_zlacpy( parsec, dplasmaUpperLower,
                    (parsec_tiled_matrix_t *)&dcB, (parsec_tiled_matrix_t *)&dcX );
    if( loud > 2 ) printf("Done\n");
    
    /* Compute */
    if( loud > 2 )printf("Compute ... ... ");
    info = dplasma_zgetrf_incpiv(parsec,
                                 (parsec_tiled_matrix_t *)&dcA,
                                 (parsec_tiled_matrix_t *)&dcL,
                                 (parsec_tiled_matrix_t *)&dcIPIV );

    if ( info == 0 ) {
        dplasma_ztrsmpl_incpiv(parsec,
                        (parsec_tiled_matrix_t *)&dcA,
                        (parsec_tiled_matrix_t *)&dcL,
                        (parsec_tiled_matrix_t *)&dcIPIV,
                        (parsec_tiled_matrix_t *)&dcX);

        dplasma_ztrsm(parsec, dplasmaLeft, dplasmaUpper,
                      dplasmaNoTrans, dplasmaNonUnit, 1.0,
                      (parsec_tiled_matrix_t *)&dcA,
                      (parsec_tiled_matrix_t *)&dcX);
    }
    if ( loud > 2 ) printf("Done\n");
    if ( info != 0 ) printf("%d: Info = %d\n", rank, info);

    /* Check the solution */
    if ( info == 0 ) {
            info_solution = check_solution( parsec, rank ? 0 : loud,
                                            (parsec_tiled_matrix_t *)&dcA0,
                                            (parsec_tiled_matrix_t *)&dcB,
                                            (parsec_tiled_matrix_t *)&dcX );
    }
    if ( rank == 0 ) {
        if ( info || info_solution ) {
            printf(" ----- TESTING ZGETRF + ZTRSMPL + ZTRSM ... FAILED !\n");
            ret |= 1;
        }
        else {
            printf(" ----- TESTING ZGETRF + ZTRSMPL + ZTRSM ....... PASSED !\n");
        }
        printf("***************************************************\n");
    }

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    parsec_data_free(dcA0.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0);
    parsec_data_free(dcB.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
    parsec_data_free(dcX.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcX);
    parsec_data_free(dcL.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcL);
    parsec_data_free(dcIPIV.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcIPIV);

    cleanup_parsec(parsec, iparam);

    return ret;
}



static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX )
{
    int info_solution;
    int N = dcB->m;
    double Rnorm, Anorm, Bnorm, Xnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    
    Anorm = dplasma_zlange(parsec, dplasmaInfNorm, dcA);
    Bnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);
    Xnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcX);

    /* Compute b - A*x */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans, -1.0, dcA, dcX, 1.0, dcB);
    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);
    result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * N * eps ) ;

    if( loud ) { 
        if( loud > 2 )
            printf( "||A||_oo = %e, ||X||_oo = %e, ||B||_oo= %e, ||A X - B||_oo = %e\n",
                    Anorm, Xnorm, Bnorm, Rnorm );

        printf("============\n");
        printf("Checking the Residual of the solution \n");
        printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);
    }
    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        info_solution = 1;
    }
    else{
        info_solution = 0;
    }

    return info_solution;
}
