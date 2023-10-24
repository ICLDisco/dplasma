/*
 * Copyright (c) 2015-2023 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "dplasma/types.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/interfaces/dtd/insert_function.h"

static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX );

static int check_inverse( parsec_context_t *parsec, int loud,
                          parsec_tiled_matrix_t *dcA,
                          parsec_tiled_matrix_t *dcInvA,
                          parsec_tiled_matrix_t *dcI );
static void warmup_zgetrf(int rank, int random_seed, parsec_context_t *parsec);

/* Global indices for the different datatypes */
static int TILE_FULL,
           TILE_RECTANGLE,
           L_TILE_RECTANGLE;

int
parsec_core_getrf_incpiv(parsec_execution_stream_t *es, parsec_task_t * this_task)
{
    (void)es;
    int m;
    int n;
    int ib;
    dplasma_complex64_t *A;
    int lda;
    int *IPIV;
    int check_info;
    int *info;

    parsec_dtd_unpack_args(this_task, &m, &n, &ib, &A, &lda, &IPIV,
                           &check_info, &info);

    CORE_zgetrf_incpiv(m, n, ib, A, lda, IPIV, info);

    if (*info != 0 && check_info)
        printf("Getrf_incpiv something is wrong\n");

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_gessm(parsec_execution_stream_t *es, parsec_task_t * this_task)
{
    (void)es;
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    dplasma_complex64_t *L;
    int ldl;
    dplasma_complex64_t *D;
    int ldd;
    dplasma_complex64_t *A;
    int lda;

    parsec_dtd_unpack_args(this_task, &m, &n, &k, &ib, &IPIV, &L, &ldl,
                           &D, &ldd, &A, &lda);

    CORE_zgessm(m, n, k, ib, IPIV, D, ldd, A, lda);

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_tstrf(parsec_execution_stream_t *es, parsec_task_t * this_task)
{
    (void)es;
    int m;
    int n;
    int ib;
    int nb;
    dplasma_complex64_t *U;
    int ldu;
    dplasma_complex64_t *A;
    int lda;
    dplasma_complex64_t *L;
    int ldl;
    int *IPIV;
    dplasma_complex64_t *WORK;
    int ldwork;
    int *check_info;
    int *info;

    parsec_dtd_unpack_args(this_task, &m, &n, &ib, &nb, &U, &ldu, &A, &lda, &L,
                           &ldl, &IPIV, &WORK, &ldwork, &check_info, &info);

    CORE_ztstrf(m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, info);

    if (*info != 0 && check_info)
        printf("Gtstrf something is wrong\n");

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_ssssm(parsec_execution_stream_t *es, parsec_task_t * this_task)
{
    (void)es;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    dplasma_complex64_t *A1;
    int lda1;
    dplasma_complex64_t *A2;
    int lda2;
    dplasma_complex64_t *L1;
    int ldl1;
    dplasma_complex64_t *L2;
    int ldl2;
    int *IPIV;

    parsec_dtd_unpack_args(this_task, &m1, &n1, &m2, &n2, &k, &ib, &A1, &lda1, &A2,
                           &lda2, &L1, &ldl1, &L2, &ldl2, &IPIV);

    CORE_zssssm(m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);

    return PARSEC_HOOK_RETURN_DONE;
}

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int info = 0;
    int ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 40, 200, 200);
    iparam[IPARAM_KP] = 4;
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'m';

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZGETRF, ((DagDouble_t)M,(DagDouble_t)N));

    LDA = max(M, LDA);

    if ( M != N && check ) {
        fprintf(stderr, "Check is impossible if M != N\n");
        check = 0;
    }
    warmup_zgetrf(rank, random_seed, parsec);

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
                               parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank, MB, NB, LDA, N, 0, 0,
                                                      M, N, P, nodes/P, KP, KQ, IP, JQ));

    /* Initializing dc for dtd */
    parsec_matrix_block_cyclic_t *__dcA = &dcA;
    parsec_dtd_data_collection_init((parsec_data_collection_t *)&dcA);

    PASTE_CODE_ALLOCATE_MATRIX(dcL, 1,
                               parsec_matrix_block_cyclic, (&dcL, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank, IB, NB, MT*IB, N, 0, 0,
                                                      MT*IB, N, P, nodes/P, KP, KQ, IP, JQ));

    /* Initializing dc for dtd */
    parsec_matrix_block_cyclic_t *__dcL = &dcL;
    parsec_dtd_data_collection_init((parsec_data_collection_t *)&dcL);

    PASTE_CODE_ALLOCATE_MATRIX(dcIPIV, 1,
                               parsec_matrix_block_cyclic, (&dcIPIV, PARSEC_MATRIX_INTEGER, PARSEC_MATRIX_TILE,
                                                      rank, MB, 1, M, NT, 0, 0,
                                                      M, NT, P, nodes/P, KP, KQ, IP, JQ));

    /* Initializing dc for dtd */
    parsec_matrix_block_cyclic_t *__dcIPIV = &dcIPIV;
    parsec_dtd_data_collection_init((parsec_data_collection_t *)&dcIPIV);

    PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
                               parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank, MB, NB, LDA, N, 0, 0,
                                                      M, N, P, nodes/P, KP, KQ, IP, JQ));
    /* Random B check */
    PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
                               parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank, MB, NB, LDB, NRHS, 0, 0,
                                                      M, NRHS, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
                               parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank, MB, NB, LDB, NRHS, 0, 0,
                                                      M, NRHS, P, nodes/P, KP, KQ, IP, JQ));
    /* Inverse check */
    PASTE_CODE_ALLOCATE_MATRIX(dcInvA, check_inv,
                               parsec_matrix_block_cyclic, (&dcInvA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank, MB, NB, LDA, N, 0, 0,
                                                      M, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcI, check_inv,
                               parsec_matrix_block_cyclic, (&dcI, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank, MB, NB, LDA, N, 0, 0,
                                                      M, N, P, nodes/P, KP, KQ, IP, JQ));

    /* matrix generation */
    if(loud > 2) printf("+++ Generate matrices ... ");
    dplasma_zpltmg( parsec, matrix_init, (parsec_tiled_matrix_t *)&dcA, random_seed );
    if ( check ) {
        dplasma_zlacpy( parsec, dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcA,
                        (parsec_tiled_matrix_t *)&dcA0 );
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, random_seed + 1 );
        dplasma_zlacpy( parsec, dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcB,
                        (parsec_tiled_matrix_t *)&dcX );
    }
    if ( check_inv ) {
        dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&dcI);
        dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&dcInvA);
    }
    if(loud > 2) printf("Done\n");

    /* Getting new parsec handle of dtd type */
    parsec_taskpool_t *dtd_tp = parsec_dtd_taskpool_new();

    /* Parameters passed on to Insert_task() */
    int k, m, n;
    int ldak, ldam;
    int tempkm, tempkn, tempmm, tempnn;
    int ib = dcL.super.mb;
    int minMNT = min(dcA.super.mt, dcA.super.nt);
    int check_info;
    int anb, nb, ldl;

    /* Allocating data arrays to be used by comm engine */
    /* A */
    parsec_arena_datatype_t *tile_full = parsec_dtd_create_arena_datatype(parsec, &TILE_FULL);
    dplasma_add2arena_tile( tile_full,
                            dcA.super.mb*dcA.super.nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, dcA.super.mb );

    /* IPIV */
    parsec_arena_datatype_t *tile_rectangle = parsec_dtd_create_arena_datatype(parsec, &TILE_RECTANGLE);
    dplasma_add2arena_rectangle( tile_rectangle,
                                 dcA.super.mb*sizeof(int),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_int_t, dcA.super.mb, 1, -1 );

    /* L */
    parsec_arena_datatype_t *l_tile_rectangle = parsec_dtd_create_arena_datatype(parsec, &L_TILE_RECTANGLE);
    dplasma_add2arena_rectangle( l_tile_rectangle,
                                 dcL.super.mb*dcL.super.nb*sizeof(dplasma_complex64_t),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_double_complex_t, dcL.super.mb, dcL.super.nb, -1);

    /* Registering the handle with parsec context */
    parsec_context_add_taskpool( parsec, dtd_tp );

    SYNC_TIME_START();

    /* #### PaRSEC context starting #### */

    /* start parsec context */
    parsec_context_start( parsec );

    /* Testing insert task function */
    for( k = 0; k < minMNT; k++ ) {
        tempkm = k == dcA.super.mt-1 ? (dcA.super.m)-k*(dcA.super.mb) : dcA.super.mb;
        tempkn = k == dcA.super.nt-1 ? (dcA.super.n)-k*(dcA.super.nb) : dcA.super.nb;
        ldak = BLKLDD((parsec_tiled_matrix_t*)&dcA, k);
        check_info = k == dcA.super.mt-1;

        parsec_dtd_insert_task( dtd_tp,     parsec_core_getrf_incpiv,             0, PARSEC_DEV_CPU, "getrf_incpiv",
                           sizeof(int),           &tempkm,                           PARSEC_VALUE,
                           sizeof(int),           &tempkn,                           PARSEC_VALUE,
                           sizeof(int),           &ib,                               PARSEC_VALUE,
                           PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, k),     PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                           sizeof(int),           &ldak,                             PARSEC_VALUE,
                           PASSED_BY_REF,         PARSEC_DTD_TILE_OF(IPIV, k, k),  PARSEC_OUTPUT | TILE_RECTANGLE,
                           sizeof(int),   &check_info,                       PARSEC_VALUE,
                           sizeof(int *),         &info,                             PARSEC_REF,
                           PARSEC_DTD_ARG_END );

        for( n = k+1; n < dcA.super.nt; n++ ) {
            tempnn = n == dcA.super.nt-1 ? (dcA.super.n)-n*(dcA.super.nb) : dcA.super.nb;
            ldl = dcL.super.mb;

            parsec_dtd_insert_task( dtd_tp,      parsec_core_gessm,           0,  PARSEC_DEV_CPU,   "gessm",
                               sizeof(int),           &tempkm,                           PARSEC_VALUE,
                               sizeof(int),           &tempnn,                           PARSEC_VALUE,
                               sizeof(int),           &tempkm,                           PARSEC_VALUE,
                               sizeof(int),           &ib,                               PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(IPIV, k, k),    PARSEC_INPUT | TILE_RECTANGLE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(L, k, k),       PARSEC_INPUT | L_TILE_RECTANGLE,
                               sizeof(int),           &ldl,                              PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, k),       PARSEC_INPUT | TILE_FULL,
                               sizeof(int),           &ldak,                             PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, n),       PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                               sizeof(int),           &ldak,                             PARSEC_VALUE,
                               PARSEC_DTD_ARG_END );
        }
        parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(L, k, k) );
        parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(IPIV, k, k) );

        for( m = k+1; m < dcA.super.mt; m++ ) {
            tempmm = m == dcA.super.mt-1 ? (dcA.super.m)-m*(dcA.super.mb) : dcA.super.mb;
            ldam = BLKLDD( (parsec_tiled_matrix_t*)&dcA, m);
            nb = dcL.super.nb;
            ldl = dcL.super.mb;
            check_info = m == dcA.super.mt-1;

            parsec_dtd_insert_task( dtd_tp,      parsec_core_tstrf,              0,      PARSEC_DEV_CPU,    "tstrf",
                               sizeof(int),           &tempmm,                           PARSEC_VALUE,
                               sizeof(int),           &tempkn,                           PARSEC_VALUE,
                               sizeof(int),           &ib,                               PARSEC_VALUE,
                               sizeof(int),           &nb,                               PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, k),     PARSEC_INOUT | TILE_FULL,
                               sizeof(int),           &ldak,                             PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, m, k),     PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                               sizeof(int),           &ldam,                             PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(L, m, k),     PARSEC_OUTPUT | L_TILE_RECTANGLE,
                               sizeof(int),           &ldl,                              PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(IPIV, m, k),  PARSEC_OUTPUT | TILE_RECTANGLE,
                               sizeof(dplasma_complex64_t)*ib*nb,    NULL,                PARSEC_SCRATCH,
                               sizeof(int),           &nb,                               PARSEC_VALUE,
                               sizeof(int),           &check_info,                       PARSEC_VALUE,
                               sizeof(int *),         &info,                             PARSEC_REF,
                               PARSEC_DTD_ARG_END );

            for( n = k+1; n < dcA.super.nt; n++ ) {
                tempnn = n == dcA.super.nt-1 ? (dcA.super.n)-n*(dcA.super.nb) : dcA.super.nb;
                anb = dcA.super.nb;
                ldl = dcL.super.mb;

                parsec_dtd_insert_task( dtd_tp,      parsec_core_ssssm,            0,         PARSEC_DEV_CPU,   "ssssm",
                                   sizeof(int),           &anb,                               PARSEC_VALUE,
                                   sizeof(int),           &tempnn,                            PARSEC_VALUE,
                                   sizeof(int),           &tempmm,                            PARSEC_VALUE,
                                   sizeof(int),           &tempnn,                            PARSEC_VALUE,
                                   sizeof(int),           &anb,                               PARSEC_VALUE,
                                   sizeof(int),           &ib,                                PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, n),     PARSEC_INOUT | TILE_FULL,
                                   sizeof(int),           &ldak,                              PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, m, n),     PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                                   sizeof(int),           &ldam,                              PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(L, m, k),     PARSEC_INPUT | L_TILE_RECTANGLE,
                                   sizeof(int),           &ldl,                               PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, m, k),     PARSEC_INPUT | TILE_FULL,
                                   sizeof(int),           &ldam,                              PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(IPIV, m, k),  PARSEC_INPUT | TILE_RECTANGLE,
                                   PARSEC_DTD_ARG_END );
            }
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(L, m, k) );
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(IPIV, m, k) );
        }
        for( n = k+1; n < dcA.super.nt; n++ ) {
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, k, n));
        }
        for( m = k+1; m < dcA.super.mt; m++ ) {
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, m, k) );
        }
        parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, k, k) );
    }

    parsec_dtd_data_flush_all( dtd_tp, (parsec_data_collection_t *)&dcA );
    parsec_dtd_data_flush_all( dtd_tp, (parsec_data_collection_t *)&dcL );
    parsec_dtd_data_flush_all( dtd_tp, (parsec_data_collection_t *)&dcIPIV );

    /* finishing all the tasks inserted, but not finishing the handle */
    parsec_taskpool_wait( dtd_tp );

    /* Waiting on all handle and turning everything off for this context */
    parsec_context_wait( parsec );

    /* #### PaRSEC context is done #### */

    SYNC_TIME_PRINT(rank, ("\tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops\n",
                           P, Q, NB, N,
                           gflops=(flops/1e9)/sync_time_elapsed));

    /* Cleaning up the parsec handle */
    parsec_taskpool_free( dtd_tp );

    if ( info != 0 ) {
        if( rank == 0 && loud ) printf("-- Factorization is suspicious (info = %d) ! \n", info );
        ret |= 1;
    }
    else if ( check ) {
        /*
         * First check with a right hand side
         */
        dplasma_zgetrs_incpiv( parsec, dplasmaNoTrans,
                               (parsec_tiled_matrix_t *)&dcA,
                               (parsec_tiled_matrix_t *)&dcL,
                               (parsec_tiled_matrix_t *)&dcIPIV,
                               (parsec_tiled_matrix_t *)&dcX );

        /* Check the solution */
        ret |= check_solution( parsec, (rank == 0) ? loud : 0,
                               (parsec_tiled_matrix_t *)&dcA0,
                               (parsec_tiled_matrix_t *)&dcB,
                               (parsec_tiled_matrix_t *)&dcX);

        /*
         * Second check with inverse
         */
        if ( check_inv ) {
            dplasma_zgetrs_incpiv( parsec, dplasmaNoTrans,
                                   (parsec_tiled_matrix_t *)&dcA,
                                   (parsec_tiled_matrix_t *)&dcL,
                                   (parsec_tiled_matrix_t *)&dcIPIV,
                                   (parsec_tiled_matrix_t *)&dcInvA );

            /* Check the solution */
            ret |= check_inverse(parsec, (rank == 0) ? loud : 0,
                                 (parsec_tiled_matrix_t *)&dcA0,
                                 (parsec_tiled_matrix_t *)&dcInvA,
                                 (parsec_tiled_matrix_t *)&dcI);
        }
    }

    if ( check ) {
        parsec_data_free(dcA0.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0);
        parsec_data_free(dcB.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
        parsec_data_free(dcX.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcX);
        if ( check_inv ) {
            parsec_data_free(dcInvA.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcInvA);
            parsec_data_free(dcI.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcI);
        }
    }

    /* Cleaning data arrays we allocated for communication */
    dplasma_matrix_del2arena( tile_full );
    parsec_dtd_destroy_arena_datatype(parsec, TILE_FULL);
    dplasma_matrix_del2arena( tile_rectangle );
    parsec_dtd_destroy_arena_datatype(parsec, TILE_RECTANGLE);
    dplasma_matrix_del2arena( l_tile_rectangle );
    parsec_dtd_destroy_arena_datatype(parsec, L_TILE_RECTANGLE);

    parsec_dtd_data_collection_fini( (parsec_data_collection_t *)&dcA );
    parsec_dtd_data_collection_fini( (parsec_data_collection_t *)&dcL );
    parsec_dtd_data_collection_fini( (parsec_data_collection_t *)&dcIPIV );

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
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
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double Bnorm = 0.0;
    double Xnorm, result;
    int m = dcB->m;
    double eps = LAPACKE_dlamch_work('e');

    Anorm = dplasma_zlange(parsec, dplasmaInfNorm, dcA);
    Bnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);
    Xnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcX);

    /* Compute b - A*x */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans, -1.0, dcA, dcX, 1.0, dcB);

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);

    result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * m * eps ) ;

    if ( loud > 2 ) {
        printf("============\n");
        printf("Checking the Residual of the solution \n");
        if ( loud > 3 )
            printf( "-- ||A||_oo = %e, ||X||_oo = %e, ||B||_oo= %e, ||A X - B||_oo = %e\n",
                    Anorm, Xnorm, Bnorm, Rnorm );
        printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);
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

static int check_inverse( parsec_context_t *parsec, int loud,
                          parsec_tiled_matrix_t *dcA,
                          parsec_tiled_matrix_t *dcInvA,
                          parsec_tiled_matrix_t *dcI )
{
    int info_solution;
    double Anorm    = 0.0;
    double InvAnorm = 0.0;
    double Rnorm, result;
    int m = dcA->m;
    double eps = LAPACKE_dlamch_work('e');

    Anorm    = dplasma_zlange(parsec, dplasmaInfNorm, dcA   );
    InvAnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcInvA);

    /* Compute I - A*A^{-1} */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans, -1.0, dcA, dcInvA, 1.0, dcI);

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcI);

    result = Rnorm / ( ( Anorm * InvAnorm ) * m * eps ) ;

    if ( loud > 2 ) {
        printf("============\n");
        printf("Checking the Residual of the solution \n");
        if ( loud > 3 )
            printf( "-- ||A||_oo = %e, ||A^{-1}||_oo = %e, ||A A^{-1} - I||_oo = %e\n",
                    Anorm, InvAnorm, Rnorm );
        printf("-- ||AA^{-1}-I||_oo/((||A||_oo||A^{-1}||_oo).N.eps) = %e \n", result);
    }

    if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        if( loud ) printf("-- Solution is suspicious ! \n");
        info_solution = 1;
    }
    else{
        if( loud ) printf("-- Solution is CORRECT ! \n");
        info_solution = 0;
    }

    return info_solution;
}

static uint32_t always_local_rank_of(parsec_data_collection_t * desc, ...)
{
    return desc->myrank;
}

static uint32_t always_local_rank_of_key(parsec_data_collection_t * desc, parsec_data_key_t key)
{
    (void)key;
    return desc->myrank;
}

static void warmup_zgetrf(int rank, int random_seed, parsec_context_t *parsec)
{
    int MB = 64;
    int IB = 40;
    int NB = 64;
    int MT = 4;
    int NT = 4;
    int N = NB*NT;
    int M = MB*MT;
    int matrix_init = dplasmaMatrixRandom;
    int info;

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, M, N, 0, 0,
                               M, N, 1, 1, 1, 1, 0, 0));
    dcA.super.super.rank_of = always_local_rank_of;
    dcA.super.super.rank_of_key = always_local_rank_of_key;
    PASTE_CODE_ALLOCATE_MATRIX(dcL, 1,
        parsec_matrix_block_cyclic, (&dcL, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, 1, 1, 1, 1, 0, 0));
    dcL.super.super.rank_of = always_local_rank_of;
    dcL.super.super.rank_of_key = always_local_rank_of_key;
    PASTE_CODE_ALLOCATE_MATRIX(dcIPIV, 1,
        parsec_matrix_block_cyclic, (&dcIPIV, PARSEC_MATRIX_INTEGER, PARSEC_MATRIX_TILE,
                               rank, MB, 1, M, NT, 0, 0,
                               M, NT, 1, 1, 1, 1, 0, 0));
    dcIPIV.super.super.rank_of = always_local_rank_of;
    dcIPIV.super.super.rank_of_key = always_local_rank_of_key;

    /* Do the CPU warmup first */
    dplasma_zpltmg( parsec, matrix_init, (parsec_tiled_matrix_t *)&dcA, random_seed );
    parsec_taskpool_t *zgetrf_incpiv = dplasma_zgetrf_incpiv_New((parsec_tiled_matrix_t*)&dcA,
                            (parsec_tiled_matrix_t*)&dcL,
                            (parsec_tiled_matrix_t*)&dcIPIV,
                            &info);
    zgetrf_incpiv->devices_index_mask = 1<<0; /* Only CPU ! */
    parsec_context_add_taskpool(parsec, zgetrf_incpiv);
    parsec_context_start(parsec);
    parsec_context_wait(parsec);

    /* Check for which device type (skipping RECURSIVE), we need to warmup this operation */
    for(int dtype = PARSEC_DEV_RECURSIVE+1; dtype < PARSEC_DEV_MAX_NB_TYPE; dtype++) {
        for(int i = 0; i < (int)zgetrf_incpiv->nb_task_classes; i++) {
            for(int j = 0; NULL != zgetrf_incpiv->task_classes_array[i]->incarnations[j].hook; j++) {
                if( zgetrf_incpiv->task_classes_array[i]->incarnations[j].type == dtype ) {
                    goto do_run; /* We found one class that was on that device, no need to try more incarnations or task classes */
                }
            }
        }
        continue; /* No incarnation of this device type on any task class; try another type */
    do_run:
        for(int did = 0; did < (int)parsec_nb_devices; did++) {
            parsec_device_module_t *dev = parsec_mca_device_get(did);
            if(dev->type != dtype)
                continue;
            /* This should work, right? Unfortunately, we can't test until there is a <dev>-enabled implementation for this test */
            for(int m = 0; m < MT; m++) {
                for(int n = 0; n < NT; n++) {
                    parsec_data_t *dta = dcA.super.super.data_of(&dcA.super.super, m, n);
                    parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                    dta = dcL.super.super.data_of(&dcL.super.super, m, n);
                    parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                    dta = dcIPIV.super.super.data_of(&dcIPIV.super.super, m, n);
                    parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                }
            }
            dplasma_zpltmg( parsec, matrix_init, (parsec_tiled_matrix_t *)&dcA, random_seed );
            parsec_taskpool_t *zgetrf_incpiv_device = dplasma_zgetrf_incpiv_New((parsec_tiled_matrix_t*)&dcA,
                                    (parsec_tiled_matrix_t*)&dcL,
                                    (parsec_tiled_matrix_t*)&dcIPIV,
                                    &info);
            parsec_context_add_taskpool(parsec, zgetrf_incpiv_device);
            parsec_context_start(parsec);
            parsec_context_wait(parsec);

            dplasma_zgetrf_incpiv_Destruct(zgetrf_incpiv_device);
        }
    }

    dplasma_zgetrf_incpiv_Destruct(zgetrf_incpiv);

}
