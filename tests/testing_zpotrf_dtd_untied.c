/*
 * Copyright (c) 2015-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

/*
 *
 * This version of potrf lets any thread insert task in PaRSEC.
 * In the other approach the main thread inserts task and the others
 * execute them. The main thread maintains a window of insertion, reaching
 * which it joins the worker thread. This version lets all thread have
 * the same kind of behavior. Here, we have a task that inserts other tasks.
 * They all insert task depending on which thread acquires the right to execute
 * the task that inserts task.
 *
 */

#include "common.h"
#include "flops.h"
#include "dplasma/types.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/interfaces/dtd/insert_function.h"

/* Global index for the full tile datatype */
static int TILE_FULL;

int
parsec_core_potrf(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int uplo;
    int m, lda, *info;
    dplasma_complex64_t *A;

    parsec_dtd_unpack_args(this_task, &uplo, &m, &A, &lda, &info);

    CORE_zpotrf(uplo, m, A, lda, info);

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_trsm(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int side, uplo, trans, diag;
    int  m, n, lda, ldc;
    dplasma_complex64_t alpha;
    dplasma_complex64_t *A, *C;

    parsec_dtd_unpack_args(this_task, &side, &uplo, &trans, &diag, &m, &n,
                           &alpha, &A, &lda, &C, &ldc);

    CORE_ztrsm(side, uplo, trans, diag,
               m, n, alpha,
               A, lda,
               C, ldc);

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_herk(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int uplo, trans;
    int m, n, lda, ldc;
    dplasma_complex64_t alpha;
    dplasma_complex64_t beta;
    dplasma_complex64_t *A;
    dplasma_complex64_t *C;

    parsec_dtd_unpack_args(this_task, &uplo, &trans, &m, &n, &alpha, &A,
                           &lda, &beta, &C, &ldc);

    CORE_zherk(uplo, trans, m, n,
               alpha, A, lda,
               beta,  C, ldc);

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_gemm(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int transA, transB;
    int m, n, k, lda, ldb, ldc;
    dplasma_complex64_t alpha, beta;
    dplasma_complex64_t *A;
    dplasma_complex64_t *B;
    dplasma_complex64_t *C;

    parsec_dtd_unpack_args(this_task, &transA, &transB, &m, &n, &k, &alpha,
                           &A, &lda, &B, &ldb, &beta, &C, &ldc);

    CORE_zgemm(transA, transB,
               m, n, k,
               alpha, A, lda,
                      B, ldb,
               beta,  C, ldc);

    return PARSEC_HOOK_RETURN_DONE;
}

int
insert_task_lower(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    /* Parameters passed on to Insert_task() */
    int tempkm, tempmm, ldak, ldam, side, transA_p, transA_g, diag, trans, transB, ldan;
    dplasma_complex64_t alpha_trsm, alpha_herk, beta;
    parsec_matrix_sym_block_cyclic_t *dcA;

    int total, *iteration, uplo, *info, m, n, k, count = 0;

    parsec_taskpool_t *dtd_tp = (parsec_taskpool_t *)this_task->taskpool;

    parsec_dtd_unpack_args(this_task, &total, &iteration, &uplo,
                           &info, &dcA);

    parsec_matrix_sym_block_cyclic_t *__dcA = dcA;

    side = dplasmaRight;
    transA_p = dplasmaConjTrans;
    diag = dplasmaNonUnit;
    alpha_trsm = 1.0;
    trans = dplasmaNoTrans;
    alpha_herk = -1.0;
    beta = 1.0;
    transB = dplasmaConjTrans;
    transA_g = dplasmaNoTrans;

    /* Testing Insert Function */
    for( k = *iteration; k < total; k++, *iteration += 1, count++ ) {
        if( count > parsec_dtd_window_size-1000 ) {
            return PARSEC_HOOK_RETURN_AGAIN;
        }

        tempkm = (k == (dcA->super.mt - 1)) ? dcA->super.m - k * dcA->super.mb : dcA->super.mb;
        ldak = BLKLDD(&dcA->super, k);

        parsec_dtd_insert_task(dtd_tp, parsec_core_potrf,
                         (total - k) * (total-k) * (total - k)/*priority*/, PARSEC_DEV_CPU, "Potrf",
                         sizeof(int),      &uplo,              PARSEC_VALUE,
                         sizeof(int),      &tempkm,            PARSEC_VALUE,
                         PASSED_BY_REF,    PARSEC_DTD_TILE_OF(A, k, k), PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                         sizeof(int),      &ldak,              PARSEC_VALUE,
                         sizeof(int *),    &info,              PARSEC_SCRATCH,
                         PARSEC_DTD_ARG_END);
        for( m = k+1; m < total; m++, count++ ) {
            tempmm = m == dcA->super.mt - 1 ? dcA->super.m - m * dcA->super.mb : dcA->super.mb;
            ldam = BLKLDD(&dcA->super, m);
            parsec_dtd_insert_task(dtd_tp, parsec_core_trsm,
                             (total - m) * (total-m) * (total - m) + 3 * ((2 * total) - k - m - 1) * (m - k)/*priority*/,
                             PARSEC_DEV_CPU, "Trsm",
                             sizeof(int),      &side,               PARSEC_VALUE,
                             sizeof(int),      &uplo,               PARSEC_VALUE,
                             sizeof(int),      &transA_p,           PARSEC_VALUE,
                             sizeof(int),      &diag,               PARSEC_VALUE,
                             sizeof(int),      &tempmm,             PARSEC_VALUE,
                             sizeof(int),      &dcA->super.nb,   PARSEC_VALUE,
                             sizeof(dplasma_complex64_t),      &alpha_trsm,         PARSEC_VALUE,
                             PASSED_BY_REF,    PARSEC_DTD_TILE_OF(A, k, k), PARSEC_INPUT | TILE_FULL,
                             sizeof(int),      &ldak,               PARSEC_VALUE,
                             PASSED_BY_REF,    PARSEC_DTD_TILE_OF(A, m, k), PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                             sizeof(int),      &ldam,               PARSEC_VALUE,
                             PARSEC_DTD_ARG_END);
        }

        parsec_dtd_data_flush(dtd_tp, PARSEC_DTD_TILE_OF(A, k, k));
        for( m = k+1; m < dcA->super.nt; m++, count++ ){
            tempmm = m == dcA->super.mt - 1 ? dcA->super.m - m * dcA->super.mb : dcA->super.mb;
            ldam = BLKLDD(&dcA->super, m);
            parsec_dtd_insert_task(dtd_tp, parsec_core_herk,
                            (total - m) * (total - m) * (total - m) + 3 * (m - k)/*priority*/,
                            PARSEC_DEV_CPU, "Herk",
                            sizeof(int),       &uplo,               PARSEC_VALUE,
                            sizeof(int),       &trans,              PARSEC_VALUE,
                            sizeof(int),       &tempmm,             PARSEC_VALUE,
                            sizeof(int),       &dcA->super.mb,   PARSEC_VALUE,
                            sizeof(dplasma_complex64_t),       &alpha_herk,         PARSEC_VALUE,
                            PASSED_BY_REF,     PARSEC_DTD_TILE_OF(A, m, k), PARSEC_INPUT | TILE_FULL,
                            sizeof(int),       &ldam,               PARSEC_VALUE,
                            sizeof(dplasma_complex64_t),       &beta,               PARSEC_VALUE,
                            PASSED_BY_REF,     PARSEC_DTD_TILE_OF(A, m, m), PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                            sizeof(int),       &ldam,               PARSEC_VALUE,
                            PARSEC_DTD_ARG_END);
            for( n = m+1; n < total; n++, count++ ){
                ldan = BLKLDD(&dcA->super, n);
                parsec_dtd_insert_task(dtd_tp,  parsec_core_gemm,
                               (total - m) * (total - m) * (total - m) + 3 * ((2 * total) - m - n - 3) * (m - n) + 6 * (m - k) /*priority*/,
                               PARSEC_DEV_CPU, "Gemm",
                               sizeof(int),        &transA_g,           PARSEC_VALUE,
                               sizeof(int),        &transB,             PARSEC_VALUE,
                               sizeof(int),        &tempmm,             PARSEC_VALUE,
                               sizeof(int),        &dcA->super.mb,   PARSEC_VALUE,
                               sizeof(int),        &dcA->super.mb,   PARSEC_VALUE,
                               sizeof(dplasma_complex64_t),        &alpha_herk,         PARSEC_VALUE,
                               PASSED_BY_REF,      PARSEC_DTD_TILE_OF(A, n, k), PARSEC_INPUT | TILE_FULL,
                               sizeof(int),        &ldan,               PARSEC_VALUE,
                               PASSED_BY_REF,      PARSEC_DTD_TILE_OF(A, m, k), PARSEC_INPUT | TILE_FULL,
                               sizeof(int),        &ldam,               PARSEC_VALUE,
                               sizeof(dplasma_complex64_t),        &beta,               PARSEC_VALUE,
                               PASSED_BY_REF,      PARSEC_DTD_TILE_OF(A, n, m), PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                               sizeof(int),        &ldan,               PARSEC_VALUE,
                               PARSEC_DTD_ARG_END);
            }
            parsec_dtd_data_flush(dtd_tp, PARSEC_DTD_TILE_OF(A, m, k));
        }
    }

    parsec_dtd_data_flush_all(dtd_tp, (parsec_data_collection_t *)dcA);
    return PARSEC_HOOK_RETURN_DONE;
}

int
insert_task_upper(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    /* Parameters passed on to Insert_task() */
    int tempkm, tempmm, ldak, ldam, side, transA_p, transA_g, diag, trans, transB, ldan;
    dplasma_complex64_t alpha_trsm, alpha_herk, beta;
    parsec_matrix_sym_block_cyclic_t *dcA;

    int total, *iteration, uplo, *info, m, n, k, count = 0;

    parsec_taskpool_t *dtd_tp = (parsec_taskpool_t *)this_task->taskpool;

    parsec_dtd_unpack_args(this_task, &total, &iteration, &uplo, &info, &dcA);

    parsec_matrix_sym_block_cyclic_t *__dcA = dcA;

    side = dplasmaLeft;
    transA_p = dplasmaConjTrans;
    diag = dplasmaNonUnit;
    alpha_trsm = 1.0;
    trans = dplasmaConjTrans;
    alpha_herk = -1.0;
    beta = 1.0;
    transB = dplasmaNoTrans;
    transA_g = dplasmaConjTrans;

    /* Testing Insert Function */
    for( k = *iteration; k < total; k++, *iteration += 1, count++ ) {
        if( count > parsec_dtd_window_size-1000 ) {
            return PARSEC_HOOK_RETURN_AGAIN;
        }

        tempkm = k == dcA->super.nt-1 ? dcA->super.n-k*dcA->super.nb : dcA->super.nb;
        ldak = BLKLDD(&dcA->super, k);
        parsec_dtd_insert_task(dtd_tp, parsec_core_potrf,
                         (total - k) * (total-k) * (total - k)/*priority*/, PARSEC_DEV_CPU, "Potrf",
                         sizeof(int),      &uplo,              PARSEC_VALUE,
                         sizeof(int),      &tempkm,            PARSEC_VALUE,
                         PASSED_BY_REF,    PARSEC_DTD_TILE_OF(A, k, k), PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                         sizeof(int),      &ldak,              PARSEC_VALUE,
                         sizeof(int *),    info,               PARSEC_SCRATCH,
                         PARSEC_DTD_ARG_END);
        for (m = k+1; m < total; m++, count++) {
            tempmm = m == dcA->super.nt-1 ? dcA->super.n-m*dcA->super.nb : dcA->super.nb;
            parsec_dtd_insert_task(dtd_tp, parsec_core_trsm,
                             (total - m) * (total-m) * (total - m) + 3 * ((2 * total) - k - m - 1) * (m - k)/*priority*/,
                             PARSEC_DEV_CPU, "Trsm",
                             sizeof(int),      &side,               PARSEC_VALUE,
                             sizeof(int),      &uplo,               PARSEC_VALUE,
                             sizeof(int),      &transA_p,           PARSEC_VALUE,
                             sizeof(int),      &diag,               PARSEC_VALUE,
                             sizeof(int),      &dcA->super.nb,   PARSEC_VALUE,
                             sizeof(int),      &tempmm,             PARSEC_VALUE,
                             sizeof(dplasma_complex64_t),      &alpha_trsm,         PARSEC_VALUE,
                             PASSED_BY_REF,    PARSEC_DTD_TILE_OF(A, k, k), PARSEC_INPUT | TILE_FULL,
                             sizeof(int),      &ldak,               PARSEC_VALUE,
                             PASSED_BY_REF,    PARSEC_DTD_TILE_OF(A, k, m), PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                             sizeof(int),      &ldak,               PARSEC_VALUE,
                             PARSEC_DTD_ARG_END);
        }
        parsec_dtd_data_flush(dtd_tp, PARSEC_DTD_TILE_OF(A, k, k));

        for (m = k+1; m < dcA->super.mt; m++, count++) {
            tempmm = m == dcA->super.nt-1 ? dcA->super.n-m*dcA->super.nb : dcA->super.nb;
            ldam = BLKLDD(&dcA->super, m);
            parsec_dtd_insert_task(dtd_tp, parsec_core_herk,
                            (total - m) * (total - m) * (total - m) + 3 * (m - k)/*priority*/,
                            PARSEC_DEV_CPU, "Herk",
                            sizeof(int),       &uplo,               PARSEC_VALUE,
                            sizeof(int),       &trans,              PARSEC_VALUE,
                            sizeof(int),       &tempmm,             PARSEC_VALUE,
                            sizeof(int),       &dcA->super.mb,   PARSEC_VALUE,
                            sizeof(dplasma_complex64_t),       &alpha_herk,         PARSEC_VALUE,
                            PASSED_BY_REF,     PARSEC_DTD_TILE_OF(A, k, m), PARSEC_INPUT | TILE_FULL,
                            sizeof(int),       &ldak,               PARSEC_VALUE,
                            sizeof(dplasma_complex64_t),    &beta,                  PARSEC_VALUE,
                            PASSED_BY_REF,     PARSEC_DTD_TILE_OF(A, m, m), PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                            sizeof(int),       &ldam,               PARSEC_VALUE,
                            PARSEC_DTD_ARG_END);
            for (n = m+1; n < total; n++, count++) {
               ldan = BLKLDD(&dcA->super, n);
               parsec_dtd_insert_task(dtd_tp,  parsec_core_gemm,
                               (total - m) * (total - m) * (total - m) + 3 * ((2 * total) - m - n - 3) * (m - n) + 6 * (m - k) /*priority*/,
                               PARSEC_DEV_CPU, "Gemm",
                               sizeof(int),        &transA_g,           PARSEC_VALUE,
                               sizeof(int),        &transB,             PARSEC_VALUE,
                               sizeof(int),        &dcA->super.mb,   PARSEC_VALUE,
                               sizeof(int),        &tempmm,             PARSEC_VALUE,
                               sizeof(int),        &dcA->super.mb,   PARSEC_VALUE,
                               sizeof(dplasma_complex64_t),        &alpha_herk,         PARSEC_VALUE,
                               PASSED_BY_REF,      PARSEC_DTD_TILE_OF(A, k, m), PARSEC_INPUT | TILE_FULL,
                               sizeof(int),        &ldak,               PARSEC_VALUE,
                               PASSED_BY_REF,      PARSEC_DTD_TILE_OF(A, k, n), PARSEC_INPUT | TILE_FULL,
                               sizeof(int),        &ldak,               PARSEC_VALUE,
                               sizeof(dplasma_complex64_t),        &beta,               PARSEC_VALUE,
                               PASSED_BY_REF,      PARSEC_DTD_TILE_OF(A, m, n), PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                               sizeof(int),        &ldan,               PARSEC_VALUE,
                               PARSEC_DTD_ARG_END);
            }
            parsec_dtd_data_flush(dtd_tp, PARSEC_DTD_TILE_OF(A, k, m));
        }
    }

    parsec_dtd_data_flush_all(dtd_tp, (parsec_data_collection_t *)dcA);
    return PARSEC_HOOK_RETURN_DONE;
}

int main(int argc, char **argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int uplo = dplasmaUpper;
    int info = 0;
    int ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 0, 180, 180);
    iparam[IPARAM_NGPUS] = DPLASMA_ERR_NOT_SUPPORTED;

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZPOTRF, ((DagDouble_t)N));

    /* Initializing matrix structure */
    LDA = dplasma_imax( LDA, N );
    LDB = dplasma_imax( LDB, N );
    KP = 1;
    KQ = 1;

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   N, N, P, nodes/P, uplo));

    parsec_dtd_data_collection_init((parsec_data_collection_t *)&dcA);

    parsec_taskpool_t *dtd_tp = parsec_dtd_taskpool_new( );

    /* Default type */
    parsec_arena_datatype_t *tile_full = parsec_dtd_create_arena_datatype(parsec, &TILE_FULL);
    dplasma_add2arena_tile( tile_full,
                            dcA.super.mb*dcA.super.nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, dcA.super.mb );

    for(int t = 0; t < nruns+1; t++) {
        /* matrix generation */
        if(loud > 3) printf("+++ Generate matrices ... ");
        dplasma_zplghe( parsec, (double)(N), uplo,
                        (parsec_tiled_matrix_t *)&dcA, random_seed);
        if(loud > 3) printf("Done\n");

        /* Registering the handle with parsec context */
        parsec_context_add_taskpool( parsec, dtd_tp );

        SYNC_TIME_START();

        /* #### PaRSEC context starting #### */

        /* Start parsec context */
        parsec_context_start(parsec);

        int *iteration = malloc(sizeof(int));
        *iteration = 0;
        int total;
        if( dplasmaLower == uplo ) {
            total = dcA.super.mt;
            parsec_dtd_insert_task( dtd_tp,       insert_task_lower, 0, PARSEC_DEV_CPU, "insert_task_lower",
                                    sizeof(int),           &total,              PARSEC_VALUE,
                                    sizeof(int),           iteration,           PARSEC_REF,
                                    sizeof(int),           &uplo,               PARSEC_VALUE,
                                    sizeof(int *),         &info,               PARSEC_REF,
                                    sizeof(parsec_matrix_sym_block_cyclic_t *), &dcA, PARSEC_REF,
                                    PARSEC_DTD_ARG_END );

        } else {
            total = dcA.super.nt;
            parsec_dtd_insert_task( dtd_tp,       insert_task_upper, 0, PARSEC_DEV_CPU, "insert_task_upper",
                                    sizeof(int),           &total,              PARSEC_VALUE,
                                    sizeof(int),           iteration,           PARSEC_REF,
                                    sizeof(int),           &uplo,               PARSEC_VALUE,
                                    sizeof(int *),         &info,               PARSEC_REF,
                                    sizeof(parsec_matrix_sym_block_cyclic_t *), &dcA, PARSEC_REF,
                                    PARSEC_DTD_ARG_END );

        }

        /* Finishing all the tasks inserted, but not finishing the handle */
        parsec_dtd_taskpool_wait( dtd_tp );

        /* Waiting on all handle and turning everything off for this context */
        parsec_context_wait(parsec);

        /* #### PaRSEC context is done #### */

        if(t > 0) {
            SYNC_TIME_PRINT(rank, ("\tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops\n",
                    P, Q, NB, N,
                    gflops=(flops/1e9)/sync_time_elapsed));
            gflops_avg += gflops/nruns;
        }
    }
    /* Cleaning up the parsec handle */
    parsec_taskpool_free( dtd_tp );
    PASTE_CODE_PERF_LOOP_DONE();

    if( 0 == rank && info != 0 ) {
        printf("-- Factorization is suspicious (info = %d) ! \n", info);
        ret |= 1;
    }
    if( !info && check ) {
        /* Check the factorization */
        PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
            parsec_matrix_sym_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                       rank, MB, NB, LDA, N, 0, 0,
                                       N, N, P, nodes/P, uplo));
        dplasma_zplghe( parsec, (double)(N), uplo,
                        (parsec_tiled_matrix_t *)&dcA0, random_seed);

        ret |= check_zpotrf( parsec, (rank == 0) ? loud : 0, uplo,
                             (parsec_tiled_matrix_t *)&dcA,
                             (parsec_tiled_matrix_t *)&dcA0);

        /* Check the solution */
        PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
            parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, NRHS, 0, 0,
                                   N, NRHS, P, nodes/P, KP, KQ, IP, JQ));
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcB, random_seed+1);

        PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
            parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                   rank, MB, NB, LDB, NRHS, 0, 0,
                                   N, NRHS, P, nodes/P, KP, KQ, IP, JQ));
        dplasma_zlacpy( parsec, dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcB, (parsec_tiled_matrix_t *)&dcX );

        dplasma_zpotrs(parsec, uplo,
                       (parsec_tiled_matrix_t *)&dcA,
                       (parsec_tiled_matrix_t *)&dcX );

        ret |= check_zaxmb( parsec, (rank == 0) ? loud : 0, uplo,
                            (parsec_tiled_matrix_t *)&dcA0,
                            (parsec_tiled_matrix_t *)&dcB,
                            (parsec_tiled_matrix_t *)&dcX);

        /* Cleanup */
        parsec_data_free(dcA0.mat); dcA0.mat = NULL;
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0 );
        parsec_data_free(dcB.mat); dcB.mat = NULL;
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB );
        parsec_data_free(dcX.mat); dcX.mat = NULL;
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcX );
    }

    parsec_dtd_destroy_arena_datatype(parsec, TILE_FULL);
    dplasma_matrix_del2arena( tile_full );
    parsec_dtd_data_collection_fini( (parsec_data_collection_t *)&dcA );
    parsec_data_free(dcA.mat); dcA.mat = NULL;
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);

    cleanup_parsec(parsec, iparam);
    return ret;
}
