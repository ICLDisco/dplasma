/*
 * Copyright (c) 2015-2024 The University of Tennessee and The University
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

static int check_orthogonality(parsec_context_t *parsec, int loud,
                               parsec_tiled_matrix_t *Q);
static int check_factorization(parsec_context_t *parsec, int loud,
                               parsec_tiled_matrix_t *Aorig,
                               parsec_tiled_matrix_t *A,
                               parsec_tiled_matrix_t *Q);
static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX );
static void warmup_zgeqrf(int rank, int random_seed, parsec_context_t *parsec);

/* Global indices for the different datatypes */
static int TILE_FULL,
           TILE_RECTANGLE;

int
parsec_core_geqrt(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int m;
    int n;
    int ib;
    dplasma_complex64_t *A;
    int lda;
    dplasma_complex64_t *T;
    int ldt;
    dplasma_complex64_t *TAU;
    dplasma_complex64_t *WORK;

    parsec_dtd_unpack_args(this_task, &m, &n, &ib, &A, &lda, &T, &ldt, &TAU, &WORK);

    CORE_zgeqrt(m, n, ib, A, lda, T, ldt, TAU, WORK);

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_unmqr(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int side;
    int trans;
    int m;
    int n;
    int k;
    int ib;
    dplasma_complex64_t *A;
    int lda;
    dplasma_complex64_t *T;
    int ldt;
    dplasma_complex64_t *C;
    int ldc;
    dplasma_complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(this_task, &side, &trans, &m, &n, &k, &ib, &A,
                           &lda, &T, &ldt, &C, &ldc, &WORK, &ldwork);

    CORE_zunmqr(side, trans, m, n, k, ib,
                A, lda, T, ldt, C, ldc, WORK, ldwork);

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_tsqrt(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int m;
    int n;
    int ib;
    dplasma_complex64_t *A1;
    int lda1;
    dplasma_complex64_t *A2;
    int lda2;
    dplasma_complex64_t *T;
    int ldt;
    dplasma_complex64_t *TAU;
    dplasma_complex64_t *WORK;

    parsec_dtd_unpack_args(this_task, &m, &n, &ib, &A1, &lda1, &A2, &lda2, &T, &ldt, &TAU, &WORK);

    CORE_ztsqrt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);

    return PARSEC_HOOK_RETURN_DONE;
}

int
parsec_core_tsmqr(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int side;
    int trans;
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
    dplasma_complex64_t *V;
    int ldv;
    dplasma_complex64_t *T;
    int ldt;
    dplasma_complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(this_task, &side, &trans, &m1, &n1, &m2, &n2, &k,
                           &ib, &A1, &lda1, &A2, &lda2, &V, &ldv, &T, &ldt, &WORK, &ldwork);

    CORE_ztsmqr(side, trans, m1, n1, m2, n2, k, ib,
                A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);

    return PARSEC_HOOK_RETURN_DONE;
}

int main(int argc, char **argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 32, 200, 200);
    iparam[IPARAM_KP] = 4;
    iparam[IPARAM_KQ] = 1;
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'m';

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZGEQRF, ((DagDouble_t)M, (DagDouble_t)N));

    LDA = max(M, LDA);
    LDB = max(M, LDB);

    warmup_zgeqrf(rank, random_seed, parsec);

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    /* Initializing dc for dtd */
    parsec_matrix_block_cyclic_t *__dcA = &dcA;
    parsec_dtd_data_collection_init((parsec_data_collection_t *)&dcA);

    PASTE_CODE_ALLOCATE_MATRIX(dcT, 1,
        parsec_matrix_block_cyclic, (&dcT, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, P, nodes/P, KP, KQ, IP, JQ));

    /* Initializing dc for dtd */
    parsec_matrix_block_cyclic_t *__dcT = &dcT;
    parsec_dtd_data_collection_init((parsec_data_collection_t *)&dcT);

    PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
        parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcQ, check,
        parsec_matrix_block_cyclic, (&dcQ, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));

    /* Check the solution */
    PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
        parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDB, NRHS, 0, 0,
                               M, NRHS, P, nodes/P, KP, KQ, IP, JQ));

    PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
        parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDB, NRHS, 0, 0,
                               M, NRHS, P, nodes/P, KP, KQ, IP, JQ));

    /* matrix generation */
    if(loud > 3) printf("+++ Generate matrices ... ");
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, 3872);
    if( check )
        dplasma_zlacpy( parsec, dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcA, (parsec_tiled_matrix_t *)&dcA0 );
    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0., (parsec_tiled_matrix_t *)&dcT);
    if(loud > 3) printf("Done\n");

    /* Getting new parsec handle of dtd type */
    parsec_taskpool_t *dtd_tp = parsec_dtd_taskpool_new();

    /* Parameters passed on to Insert_task() */
    int k, m, n;
    int ldak, ldam;
    int tempkm, tempkn, tempnn, tempmm;
    int ib = dcT.super.mb;
    int minMNT = min(dcA.super.mt, dcA.super.nt);
    int side = dplasmaLeft,
        trans = dplasmaConjTrans;

    /* Allocating data arrays to be used by comm engine */
    /* Default type */
    parsec_arena_datatype_t *tile_full = PARSEC_OBJ_NEW(parsec_arena_datatype_t);
    dplasma_add2arena_tile( tile_full,
                            dcA.super.mb*dcA.super.nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, dcA.super.mb );
    parsec_dtd_attach_arena_datatype(parsec, tile_full, &TILE_FULL);

    parsec_arena_datatype_t *tile_rectangle = PARSEC_OBJ_NEW(parsec_arena_datatype_t);
    dplasma_add2arena_rectangle( tile_rectangle,
                                 dcT.super.mb*dcT.super.nb*sizeof(dplasma_complex64_t),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_double_complex_t, dcT.super.mb, dcT.super.nb, -1);
    parsec_dtd_attach_arena_datatype(parsec, tile_rectangle, &TILE_RECTANGLE);

    /* Registering the handle with parsec context */
    parsec_context_add_taskpool(parsec, dtd_tp);

    SYNC_TIME_START();
    /* start parsec context */
    parsec_context_start(parsec);

    /* Testing Insert Function */
    for( k = 0; k < minMNT; k++ ) {
        tempkm = k == dcA.super.mt-1 ? dcA.super.m-(k*dcA.super.mb) : dcA.super.mb;
        tempkn = k == dcA.super.nt-1 ? dcA.super.n-(k*dcA.super.nb) : dcA.super.nb;
        ldak = BLKLDD( (parsec_tiled_matrix_t*)&dcA, k);

        parsec_dtd_insert_task( dtd_tp,      parsec_core_geqrt,
                          (dcA.super.nt-k)*(dcA.super.nt-k)*(dcA.super.nt-k),
                          PARSEC_DEV_CPU, "geqrt",
                           sizeof(int),           &tempkm,                        PARSEC_VALUE,
                           sizeof(int),           &tempkn,                        PARSEC_VALUE,
                           sizeof(int),           &ib,                            PARSEC_VALUE,
                           PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, k),     PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                           sizeof(int),           &ldak,                          PARSEC_VALUE,
                           PASSED_BY_REF,         PARSEC_DTD_TILE_OF(T, k, k),     PARSEC_OUTPUT | TILE_RECTANGLE,
                           sizeof(int),           &dcT.super.mb,                  PARSEC_VALUE,
                           sizeof(dplasma_complex64_t)*dcT.super.nb,       NULL,   PARSEC_SCRATCH,
                           sizeof(dplasma_complex64_t)*ib*dcT.super.nb,    NULL,   PARSEC_SCRATCH,
                           PARSEC_DTD_ARG_END );

        for( n = k+1; n < dcA.super.nt; n++ ) {
            tempnn = n == dcA.super.nt-1 ? dcA.super.n-(n*dcA.super.nb) : dcA.super.nb;

            parsec_dtd_insert_task( dtd_tp,      parsec_core_unmqr,          0,    PARSEC_DEV_CPU,      "unmqr",
                               sizeof(int),           &side,                              PARSEC_VALUE,
                               sizeof(int),           &trans,                             PARSEC_VALUE,
                               sizeof(int),           &tempkm,                            PARSEC_VALUE,
                               sizeof(int),           &tempnn,                            PARSEC_VALUE,
                               sizeof(int),           &tempkm,                            PARSEC_VALUE,
                               sizeof(int),           &ib,                                PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, k),      PARSEC_INPUT | TILE_FULL,
                               sizeof(int),           &ldak,                              PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(T, k, k),      PARSEC_INPUT | TILE_RECTANGLE,
                               sizeof(int),           &dcT.super.mb,                   PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, n),      PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                               sizeof(int),           &ldak,                              PARSEC_VALUE,
                               sizeof(dplasma_complex64_t)*ib*dcT.super.nb,   NULL,        PARSEC_SCRATCH,
                               sizeof(int),           &dcT.super.nb,                      PARSEC_VALUE,
                               PARSEC_DTD_ARG_END );
        }
        parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(T, k, k) );

        for( m = k+1; m < dcA.super.mt; m++ ) {
            tempmm = m == dcA.super.mt-1 ? dcA.super.m-(m*dcA.super.mb) : dcA.super.mb;
            ldam = BLKLDD( (parsec_tiled_matrix_t*)&dcA, m);

            parsec_dtd_insert_task( dtd_tp,      parsec_core_tsqrt,
                              (dcA.super.mt-k)*(dcA.super.mt-k)*(dcA.super.mt-k),  PARSEC_DEV_CPU, "tsqrt",
                               sizeof(int),           &tempmm,                            PARSEC_VALUE,
                               sizeof(int),           &tempkn,                            PARSEC_VALUE,
                               sizeof(int),           &ib,                                PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, k),     PARSEC_INOUT | TILE_FULL,
                               sizeof(int),           &ldak,                              PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, m, k),     PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                               sizeof(int),           &ldam,                              PARSEC_VALUE,
                               PASSED_BY_REF,         PARSEC_DTD_TILE_OF(T, m, k),     PARSEC_OUTPUT | TILE_RECTANGLE,
                               sizeof(int),           &dcT.super.mb,                     PARSEC_VALUE,
                               sizeof(dplasma_complex64_t)*dcT.super.nb,       NULL,      PARSEC_SCRATCH,
                               sizeof(dplasma_complex64_t)*ib*dcT.super.nb,    NULL,      PARSEC_SCRATCH,
                               PARSEC_DTD_ARG_END );

            for( n = k+1; n < dcA.super.nt; n++ ) {
                tempnn = n == dcA.super.nt-1 ? dcA.super.n-(n*dcA.super.nb) : dcA.super.nb;
                int ldwork = dplasmaLeft == dplasmaLeft ? ib : dcT.super.nb;

                parsec_dtd_insert_task( dtd_tp,      parsec_core_tsmqr,
                                  (dcA.super.mt-k)*(dcA.super.mt-n)*(dcA.super.mt-n),        PARSEC_DEV_CPU,     "tsmqr",
                                   sizeof(int),           &side,                             PARSEC_VALUE,
                                   sizeof(int),           &trans,                            PARSEC_VALUE,
                                   sizeof(int),           &dcA.super.mb,                     PARSEC_VALUE,
                                   sizeof(int),           &tempnn,                           PARSEC_VALUE,
                                   sizeof(int),           &tempmm,                           PARSEC_VALUE,
                                   sizeof(int),           &tempnn,                           PARSEC_VALUE,
                                   sizeof(int),           &dcA.super.nb,                     PARSEC_VALUE,
                                   sizeof(int),           &ib,                               PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, k, n),     PARSEC_INOUT | TILE_FULL,
                                   sizeof(int),           &ldak,                             PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, m, n),     PARSEC_INOUT | TILE_FULL | PARSEC_AFFINITY,
                                   sizeof(int),           &ldam,                             PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(A, m, k),     PARSEC_INPUT | TILE_FULL,
                                   sizeof(int),           &ldam,                             PARSEC_VALUE,
                                   PASSED_BY_REF,         PARSEC_DTD_TILE_OF(T, m, k),     PARSEC_INPUT | TILE_RECTANGLE,
                                   sizeof(int),           &dcT.super.mb,                     PARSEC_VALUE,
                                   sizeof(dplasma_complex64_t)*ib*dcT.super.nb,    NULL,      PARSEC_SCRATCH,
                                   sizeof(int),           &ldwork,                           PARSEC_VALUE,
                                   PARSEC_DTD_ARG_END );
            }
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(T, m, k) );
        }
        for( n = k+1; n < dcA.super.nt; n++ ) {
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, k, n) );
        }
        for( m = k+1; m < dcA.super.mt; m++ ) {
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, m, k) );
        }
        parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, k, k) );
    }

    parsec_dtd_data_flush_all( dtd_tp, (parsec_data_collection_t *)&dcA );
    parsec_dtd_data_flush_all( dtd_tp, (parsec_data_collection_t *)&dcT );

    /* finishing all the tasks inserted, but not finishing the handle */
    parsec_taskpool_wait( dtd_tp );

    /* Waiting on all handle and turning everything off for this context */
    parsec_context_wait( parsec );

    SYNC_TIME_PRINT(rank, ("\tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops\n",
                           P, Q, NB, N,
                           gflops=(flops/1e9)/sync_time_elapsed));

    /* Cleaning up the parsec handle */
    parsec_taskpool_free( dtd_tp );

    if( check ) {
        if (M >= N) {
            if(loud > 2) printf("+++ Generate the Q ...");
            dplasma_zungqr( parsec,
                            (parsec_tiled_matrix_t *)&dcA,
                            (parsec_tiled_matrix_t *)&dcT,
                            (parsec_tiled_matrix_t *)&dcQ);
            if(loud > 2) printf("Done\n");

            if(loud > 2) printf("+++ Solve the system ...");
            dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcX, 2354);
            dplasma_zlacpy( parsec, dplasmaUpperLower,
                            (parsec_tiled_matrix_t *)&dcX,
                            (parsec_tiled_matrix_t *)&dcB );
            dplasma_zgeqrs( parsec,
                            (parsec_tiled_matrix_t *)&dcA,
                            (parsec_tiled_matrix_t *)&dcT,
                            (parsec_tiled_matrix_t *)&dcX );
            if(loud > 2) printf("Done\n");

            /* Check the orthogonality, factorization and the solution */
            ret |= check_orthogonality( parsec, (rank == 0) ? loud : 0,
                                        (parsec_tiled_matrix_t *)&dcQ);
            ret |= check_factorization( parsec, (rank == 0) ? loud : 0,
                                        (parsec_tiled_matrix_t *)&dcA0,
                                        (parsec_tiled_matrix_t *)&dcA,
                                        (parsec_tiled_matrix_t *)&dcQ );
            ret |= check_solution( parsec, (rank == 0) ? loud : 0,
                                   (parsec_tiled_matrix_t *)&dcA0,
                                   (parsec_tiled_matrix_t *)&dcB,
                                   (parsec_tiled_matrix_t *)&dcX );

        } else {
            printf("Check cannot be performed when N > M\n");
        }

        parsec_data_free(dcA0.mat);
        parsec_data_free(dcQ.mat);
        parsec_data_free(dcB.mat);
        parsec_data_free(dcX.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcQ);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcX);
    }

    /* Cleaning data arrays we allocated for communication */
    dplasma_matrix_del2arena( tile_full );
    parsec_dtd_free_arena_datatype(parsec, TILE_FULL);
    dplasma_matrix_del2arena( tile_rectangle );
    parsec_dtd_free_arena_datatype(parsec, TILE_RECTANGLE);

    parsec_dtd_data_collection_fini( (parsec_data_collection_t *)&dcA );
    parsec_dtd_data_collection_fini( (parsec_data_collection_t *)&dcT );

    parsec_data_free(dcA.mat);
    parsec_data_free(dcT.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcT);

    cleanup_parsec(parsec, iparam);

    return ret;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int check_orthogonality(parsec_context_t *parsec, int loud, parsec_tiled_matrix_t *Q)
{
    parsec_matrix_block_cyclic_t *twodQ = (parsec_matrix_block_cyclic_t *)Q;
    double normQ = 999999.0;
    double result;
    double eps = LAPACKE_dlamch_work('e');
    int info_ortho;
    int M = Q->m;
    int N = Q->n;
    int minMN = min(M, N);

    PASTE_CODE_ALLOCATE_MATRIX(Id, 1,
        parsec_matrix_block_cyclic, (&Id, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               twodQ->grid.rank,
                               Q->mb, Q->nb, minMN, minMN, 0, 0,
                               minMN, minMN, twodQ->grid.rows, twodQ->grid.cols, twodQ->grid.krows, twodQ->grid.kcols, twodQ->grid.ip, twodQ->grid.jq));

    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&Id);

    /* Perform Id - Q'Q */
    if ( M >= N ) {
        dplasma_zherk( parsec, dplasmaUpper, dplasmaConjTrans,
                       1.0, Q, -1.0, (parsec_tiled_matrix_t*)&Id );
    } else {
        dplasma_zherk( parsec, dplasmaUpper, dplasmaNoTrans,
                       1.0, Q, -1.0, (parsec_tiled_matrix_t*)&Id );
    }

    normQ = dplasma_zlanhe(parsec, dplasmaInfNorm, dplasmaUpper, (parsec_tiled_matrix_t*)&Id);

    result = normQ / (minMN * eps);
    if ( loud ) {
        printf("============\n");
        printf("Checking the orthogonality of Q \n");
        printf("||Id-Q'*Q||_oo / (N*eps) = %e \n", result);
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        if ( loud ) printf("-- Orthogonality is suspicious ! \n");
        info_ortho=1;
    }
    else {
        if ( loud ) printf("-- Orthogonality is CORRECT ! \n");
        info_ortho=0;
    }

    parsec_data_free(Id.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&Id);
    return info_ortho;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int
check_factorization(parsec_context_t *parsec, int loud,
                    parsec_tiled_matrix_t *Aorig,
                    parsec_tiled_matrix_t *A,
                    parsec_tiled_matrix_t *Q)
{
    parsec_tiled_matrix_t *subA;
    parsec_matrix_block_cyclic_t *twodA = (parsec_matrix_block_cyclic_t *)A;
    double Anorm, Rnorm;
    double result;
    double eps = LAPACKE_dlamch_work('e');
    int info_factorization;
    int M = A->m;
    int N = A->n;
    int minMN = min(M, N);

    PASTE_CODE_ALLOCATE_MATRIX(Residual, 1,
        parsec_matrix_block_cyclic, (&Residual, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               twodA->grid.rank,
                               A->mb, A->nb, M, N, 0, 0,
                               M, N, twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols, twodA->grid.ip, twodA->grid.jq));

    PASTE_CODE_ALLOCATE_MATRIX(R, 1,
        parsec_matrix_block_cyclic, (&R, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               twodA->grid.rank,
                               A->mb, A->nb, N, N, 0, 0,
                               N, N, twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols, twodA->grid.ip, twodA->grid.jq));

    /* Copy the original A in Residual */
    dplasma_zlacpy( parsec, dplasmaUpperLower, Aorig, (parsec_tiled_matrix_t *)&Residual );

    /* Extract the R */
    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0., (parsec_tiled_matrix_t *)&R);

    subA = parsec_tiled_matrix_submatrix( A, 0, 0, N, N );
    dplasma_zlacpy( parsec, dplasmaUpper, subA, (parsec_tiled_matrix_t *)&R );
    free(subA);

    /* Perform Residual = Aorig - Q*R */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans,
                   -1.0, Q, (parsec_tiled_matrix_t *)&R,
                    1.0, (parsec_tiled_matrix_t *)&Residual);

    /* Free R */
    parsec_data_free(R.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&R);

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, (parsec_tiled_matrix_t*)&Residual);
    Anorm = dplasma_zlange(parsec, dplasmaInfNorm, Aorig);

    result = Rnorm / ( Anorm * minMN * eps);

    if ( loud ) {
        printf("============\n");
        printf("Checking the QR Factorization \n");
        printf("-- ||A-QR||_oo/(||A||_oo.N.eps) = %e \n", result );
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        if ( loud ) printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else {
        if ( loud ) printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    parsec_data_free(Residual.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&Residual);
    return info_factorization;
}

static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcB,
                           parsec_tiled_matrix_t *dcX )
{
    parsec_tiled_matrix_t *subX;
    int info_solution;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double Bnorm = 0.0;
    double Xnorm, result;
    double eps = LAPACKE_dlamch_work('e');

    subX = parsec_tiled_matrix_submatrix( dcX, 0, 0, dcA->n, dcX->n );

    Anorm = dplasma_zlange(parsec, dplasmaInfNorm, dcA);
    Bnorm = dplasma_zlange(parsec, dplasmaInfNorm, dcB);
    Xnorm = dplasma_zlange(parsec, dplasmaInfNorm, subX);

    /* Compute A*x-b */
    dplasma_zgemm( parsec, dplasmaNoTrans, dplasmaNoTrans, 1.0, dcA, subX, -1.0, dcB);

    /* Compute A' * ( A*x - b ) */
    dplasma_zgemm( parsec, dplasmaConjTrans, dplasmaNoTrans,
                   1.0, dcA, dcB, 0., subX );

    Rnorm = dplasma_zlange( parsec, dplasmaInfNorm, subX );
    free(subX);

    result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * dcA->n * eps ) ;

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

static uint32_t always_local_rank_of(parsec_data_collection_t * desc, ...)
{
    return desc->myrank;
}

static uint32_t always_local_rank_of_key(parsec_data_collection_t * desc, parsec_data_key_t key)
{
    (void)key;
    return desc->myrank;
}

static void warmup_zgeqrf(int rank, int random_seed, parsec_context_t *parsec)
{
    int IB = 32;
    int MB = 64;
    int NB = 64;
    int MT = 4;
    int NT = 4;
    int N = NB*NT;
    int M = MB*MT;
    int did;

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, M, N, 0, 0,
                               M, N, 1, 1, 1, 1, 0, 0));
    dcA.super.super.rank_of = always_local_rank_of;
    dcA.super.super.rank_of_key = always_local_rank_of_key;
    PASTE_CODE_ALLOCATE_MATRIX(dcT, 1,
        parsec_matrix_block_cyclic, (&dcT, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, IB, NB, MT*IB, N, 0, 0,
                               MT*IB, N, 1, 1, 1, 1, 0, 0));
    dcT.super.super.rank_of = always_local_rank_of;
    dcT.super.super.rank_of_key = always_local_rank_of_key;

    dplasma_zpltmg( parsec, dplasmaMatrixRandom, (parsec_tiled_matrix_t *)&dcA, random_seed );
    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0., (parsec_tiled_matrix_t *)&dcT);
    parsec_taskpool_t *zgeqrf = dplasma_zgeqrf_New( (parsec_tiled_matrix_t*)&dcA, (parsec_tiled_matrix_t*)&dcT);
    zgeqrf->devices_index_mask = 1<<0; /* Only CPU ! */
    parsec_context_add_taskpool(parsec, zgeqrf);
    parsec_context_start(parsec);
    parsec_context_wait(parsec);
    dplasma_zgeqrf_Destruct(zgeqrf);

    /* We know that there is a GPU-enabled version of this operation, so warm it up if some device is enabled */
    for(did = 1; did < (int)parsec_nb_devices; did++) {
        if(!parsec_mca_device_is_gpu(did)) continue;
        for(int i = 0; i < MT; i++) {
            for(int j = 0; j < NT; j++) {
                parsec_data_t *dta = dcA.super.super.data_of(&dcA.super.super, i, j);
                parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                dta = dcT.super.super.data_of(&dcT.super.super, i, j);
                parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
            }
        }
        dplasma_zpltmg( parsec, dplasmaMatrixRandom, (parsec_tiled_matrix_t *)&dcA, random_seed );
        dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0., (parsec_tiled_matrix_t *)&dcT);
        dplasma_zgeqrf( parsec, (parsec_tiled_matrix_t*)&dcA, (parsec_tiled_matrix_t*)&dcT);
        parsec_devices_release_memory();
    }

    parsec_data_free(dcA.mat); dcA.mat = NULL;
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA );
    parsec_data_free(dcT.mat); dcT.mat = NULL;
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcT );
    parsec_devices_reset_load(parsec);
}
