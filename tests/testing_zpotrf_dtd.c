/*
 * Copyright (c) 2013-2023 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "flops.h"
#include "dplasma/types.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/interfaces/dtd/insert_function.h"

#include "dtd_wrappers/dplasma_z_dtd.h"

static void warmup_zpotrf(int rank, dplasma_enum_t uplo, int random_seed, parsec_context_t *parsec);

/* Global index for the full tile datatype */
static int TILE_FULL;

int main(int argc, char **argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int uplo = dplasmaUpper;
    int info = 0;
    int ret = 0;

    int m, n, k, total; /* loop counter */
    /* Parameters passed on to Insert_task() */
    int tempkm, tempmm, ldak, ldam, side, transA_p, transA_g, diag, trans, transB, ldan;
    dplasma_complex64_t alpha_trsm, beta;
    double alpha_herk, beta_herk;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 0, 180, 180);

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZPOTRF, ((DagDouble_t)N));

    /* initializing matrix structure */
    LDA = dplasma_imax( LDA, N );
    LDB = dplasma_imax( LDB, N );
    KP = 1;
    KQ = 1;

    /* warm up */
    if(loud > 3) printf("+++ warm up ... ");
    warmup_zpotrf(rank, uplo, random_seed, parsec);
    if(loud > 3) printf("Done\n");

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   N, N, P, nodes/P, uplo));

    /* Initializing dc for dtd */
    parsec_matrix_sym_block_cyclic_t *__dcA = &dcA;
    parsec_dtd_data_collection_init((parsec_data_collection_t *)&dcA);

    /* matrix generation */
    if(loud > 3) printf("+++ Generate matrices ... ");
    dplasma_zplghe( parsec, (double)(N), uplo,
                    (parsec_tiled_matrix_t *)&dcA, random_seed);
    if(loud > 3) printf("Done\n");

    /* Getting new parsec handle of dtd type */
    parsec_taskpool_t *dtd_tp = parsec_dtd_taskpool_new();

    /* Allocating data arrays to be used by comm engine */
    parsec_arena_datatype_t *tile_full = parsec_dtd_create_arena_datatype(parsec, &TILE_FULL);
    dplasma_add2arena_tile( tile_full,
                            dcA.super.mb*dcA.super.nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, dcA.super.mb );

    /* Registering the handle with parsec context */
    parsec_context_add_taskpool( parsec, dtd_tp );

    SYNC_TIME_START();

    /* #### parsec context Starting #### */

    /* start parsec context */
    parsec_context_start( parsec );
    parsec_task_class_t *zpotrf_tc = parsec_dtd_create_zpotrf_task_class(dtd_tp, TILE_FULL, PARSEC_DEV_ALL);
    parsec_task_class_t *ztrsm_tc  = parsec_dtd_create_ztrsm_task_class( dtd_tp, TILE_FULL, PARSEC_DEV_ALL);
    parsec_task_class_t *zherk_tc  = parsec_dtd_create_zherk_task_class( dtd_tp, TILE_FULL, PARSEC_DEV_ALL);
    parsec_task_class_t *zgemm_tc  = parsec_dtd_create_zgemm_task_class( dtd_tp, TILE_FULL, PARSEC_DEV_ALL);

#if defined(DPLASMA_HAVE_CUDA)
    zpotrf_dtd_workspace_info_t* infos = (zpotrf_dtd_workspace_info_t*) malloc(sizeof(zpotrf_dtd_workspace_info_t));
    if( gpus > 0 ){
        CuHI = parsec_info_lookup(&parsec_per_stream_infos, "DPLASMA::CUDA::HANDLES", NULL);
        assert(CuHI != -1);
    }
#endif

    if( dplasmaLower == uplo ) {

        side = dplasmaRight;
        transA_p = dplasmaConjTrans;
        diag = dplasmaNonUnit;
        alpha_trsm = 1.0;
        trans = dplasmaNoTrans;
        alpha_herk = -1.0;
        beta_herk = 1.0;
        beta = 1.0;
        transB = dplasmaConjTrans;
        transA_g = dplasmaNoTrans;

        total = dcA.super.mt;

#if defined(DPLASMA_HAVE_CUDA)
        infos->mb   = dcA.super.mb;
        infos->nb   = dcA.super.nb;
        infos->uplo = uplo;
        WoSI = parsec_info_register( &parsec_per_stream_infos, "DPLASMA::ZPOTRF::WS",
                                     zpotrf_dtd_destroy_workspace, NULL,
                                     zpotrf_dtd_create_workspace, infos,
                                     NULL);
#endif
        /* Testing Insert Function */
        for( k = 0; k < total; k++ ) {
            tempkm = (k == (dcA.super.mt - 1)) ? dcA.super.m - k * dcA.super.mb : dcA.super.mb;
            ldak = BLKLDD(&dcA.super, k);

            parsec_dtd_insert_task_with_task_class( dtd_tp, zpotrf_tc,
                               (total - k) * (total-k) * (total - k)/*priority*/,
                               PARSEC_DEV_ALL,
                               PARSEC_DTD_EMPTY_FLAG, &uplo,
                               PARSEC_DTD_EMPTY_FLAG, &tempkm,
                               PARSEC_PUSHOUT,        PARSEC_DTD_TILE_OF(A, k, k),
                               PARSEC_DTD_EMPTY_FLAG, &ldak,
                               PARSEC_DTD_EMPTY_FLAG, &info,
                               PARSEC_DTD_ARG_END );

            for( m = k+1; m < total; m++ ) {
                tempmm = m == dcA.super.mt - 1 ? dcA.super.m - m * dcA.super.mb : dcA.super.mb;
                ldam = BLKLDD(&dcA.super, m);
                parsec_dtd_insert_task_with_task_class( dtd_tp, ztrsm_tc,
                                   (total - m) * (total-m) * (total - m) + 3 * ((2 * total) - k - m - 1) * (m - k)/*priority*/,
                                   PARSEC_DEV_ALL,
                                   PARSEC_DTD_EMPTY_FLAG, &side,
                                   PARSEC_DTD_EMPTY_FLAG, &uplo,
                                   PARSEC_DTD_EMPTY_FLAG, &transA_p,
                                   PARSEC_DTD_EMPTY_FLAG, &diag,
                                   PARSEC_DTD_EMPTY_FLAG, &tempmm,
                                   PARSEC_DTD_EMPTY_FLAG, &dcA.super.nb,
                                   PARSEC_DTD_EMPTY_FLAG, &alpha_trsm,
                                   PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, k, k),
                                   PARSEC_DTD_EMPTY_FLAG, &ldak,
                                   PARSEC_PUSHOUT,        PARSEC_DTD_TILE_OF(A, m, k),
                                   PARSEC_DTD_EMPTY_FLAG, &ldam,
                                   PARSEC_DTD_ARG_END );
            }
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, k, k) );

            for( m = k+1; m < dcA.super.nt; m++ ) {
                tempmm = m == dcA.super.mt - 1 ? dcA.super.m - m * dcA.super.mb : dcA.super.mb;
                ldam = BLKLDD(&dcA.super, m);
                parsec_dtd_insert_task_with_task_class( dtd_tp, zherk_tc,
                                   (total - m) * (total - m) * (total - m) + 3 * (m - k)/*priority*/,
                                   PARSEC_DEV_ALL,
                                   PARSEC_DTD_EMPTY_FLAG, &uplo,
                                   PARSEC_DTD_EMPTY_FLAG, &trans,
                                   PARSEC_DTD_EMPTY_FLAG, &tempmm,
                                   PARSEC_DTD_EMPTY_FLAG, &dcA.super.mb,
                                   PARSEC_DTD_EMPTY_FLAG, &alpha_herk,
                                   PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, m, k),
                                   PARSEC_DTD_EMPTY_FLAG, &ldam,
                                   PARSEC_DTD_EMPTY_FLAG, &beta_herk,
                                   PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, m, m),
                                   PARSEC_DTD_EMPTY_FLAG, &ldam,
                                   PARSEC_DTD_ARG_END );

                for( n = m+1; n < total; n++ ) {
                    ldan = BLKLDD(&dcA.super, n);
                    parsec_dtd_insert_task_with_task_class( dtd_tp,  zgemm_tc,
                                       (total - m) * (total - m) * (total - m) + 3 * ((2 * total) - m - n - 3) * (m - n) + 6 * (m - k) /*priority*/,
                                       PARSEC_DEV_ALL,
                                       PARSEC_DTD_EMPTY_FLAG, &transA_g,
                                       PARSEC_DTD_EMPTY_FLAG, &transB,
                                       PARSEC_DTD_EMPTY_FLAG, &tempmm,
                                       PARSEC_DTD_EMPTY_FLAG, &dcA.super.mb,
                                       PARSEC_DTD_EMPTY_FLAG, &dcA.super.mb,
                                       PARSEC_DTD_EMPTY_FLAG, &alpha_herk,
                                       PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, n, k),
                                       PARSEC_DTD_EMPTY_FLAG, &ldan,
                                       PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, m, k),
                                       PARSEC_DTD_EMPTY_FLAG, &ldam,
                                       PARSEC_DTD_EMPTY_FLAG, &beta,
                                       PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, n, m),
                                       PARSEC_DTD_EMPTY_FLAG, &ldan,
                                       PARSEC_DTD_ARG_END );
                }
                parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, m, k) );
            }
        }
    } else {
        side = dplasmaLeft;
        transA_p = dplasmaConjTrans;
        diag = dplasmaNonUnit;
        alpha_trsm = 1.0;
        trans = dplasmaConjTrans;
        alpha_herk = -1.0;
        beta_herk = 1.0;
        beta = 1.0;
        transB = dplasmaNoTrans;
        transA_g = dplasmaConjTrans;

        total = dcA.super.nt;

#if defined(DPLASMA_HAVE_CUDA)
        infos->mb   = dcA.super.mb;
        infos->nb   = dcA.super.nb;
        infos->uplo = uplo;
        WoSI = parsec_info_register( &parsec_per_device_infos, "DPLASMA::ZPOTRF::WS",
                                     zpotrf_dtd_destroy_workspace, NULL,
                                     zpotrf_dtd_create_workspace, infos,
                                     NULL);
#endif
        for( k = 0; k < total; k++ ) {
            tempkm = k == dcA.super.nt-1 ? dcA.super.n-k*dcA.super.nb : dcA.super.nb;
            ldak = BLKLDD(&dcA.super, k);
            parsec_dtd_insert_task_with_task_class( dtd_tp, zpotrf_tc,
                               (total - k) * (total-k) * (total - k)/*priority*/,
                               PARSEC_DEV_ALL,
                               PARSEC_DTD_EMPTY_FLAG, &uplo,
                               PARSEC_DTD_EMPTY_FLAG, &tempkm,
                               PARSEC_PUSHOUT,        PARSEC_DTD_TILE_OF(A, k, k),
                               PARSEC_DTD_EMPTY_FLAG, &ldak,
                               PARSEC_DTD_EMPTY_FLAG, &info,
                               PARSEC_DTD_ARG_END );

            for( m = k+1; m < total; m++ ) {
                tempmm = m == dcA.super.nt-1 ? dcA.super.n-m*dcA.super.nb : dcA.super.nb;
                parsec_dtd_insert_task_with_task_class( dtd_tp, ztrsm_tc,
                                   (total - m) * (total-m) * (total - m) + 3 * ((2 * total) - k - m - 1) * (m - k)/*priority*/,
                                   PARSEC_DEV_ALL,
                                   PARSEC_DTD_EMPTY_FLAG, &side,
                                   PARSEC_DTD_EMPTY_FLAG, &uplo,
                                   PARSEC_DTD_EMPTY_FLAG, &transA_p,
                                   PARSEC_DTD_EMPTY_FLAG, &diag,
                                   PARSEC_DTD_EMPTY_FLAG, &dcA.super.nb,
                                   PARSEC_DTD_EMPTY_FLAG, &tempmm,
                                   PARSEC_DTD_EMPTY_FLAG, &alpha_trsm,
                                   PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, k, k),
                                   PARSEC_DTD_EMPTY_FLAG, &ldak,
                                   PARSEC_PUSHOUT,        PARSEC_DTD_TILE_OF(A, k, m),
                                   PARSEC_DTD_EMPTY_FLAG, &ldak,
                                   PARSEC_DTD_ARG_END );
            }
            parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, k, k) );

            for( m = k+1; m < dcA.super.mt; m++ ) {
                tempmm = m == dcA.super.nt-1 ? dcA.super.n-m*dcA.super.nb : dcA.super.nb;
                ldam = BLKLDD(&dcA.super, m);
                parsec_dtd_insert_task_with_task_class( dtd_tp, zherk_tc,
                                   (total - m) * (total - m) * (total - m) + 3 * (m - k)/*priority*/,
                                   PARSEC_DEV_ALL,
                                   PARSEC_DTD_EMPTY_FLAG, &uplo,
                                   PARSEC_DTD_EMPTY_FLAG, &trans,
                                   PARSEC_DTD_EMPTY_FLAG, &tempmm,
                                   PARSEC_DTD_EMPTY_FLAG, &dcA.super.mb,
                                   PARSEC_DTD_EMPTY_FLAG, &alpha_herk,
                                   PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, k, m),
                                   PARSEC_DTD_EMPTY_FLAG, &ldak,
                                   PARSEC_DTD_EMPTY_FLAG, &beta_herk,
                                   PARSEC_DTD_EMPTY_FLAG, PARSEC_DTD_TILE_OF(A, m, m),
                                   PARSEC_DTD_EMPTY_FLAG, &ldam,
                                   PARSEC_DTD_ARG_END );

                for( n = m+1; n < total; n++ ) {
                   ldan = BLKLDD(&dcA.super, n);
                   parsec_dtd_insert_task_with_task_class( dtd_tp,  zgemm_tc,
                                      (total - m) * (total - m) * (total - m) + 3 * ((2 * total) - m - n - 3) * (m - n) + 6 * (m - k) /*priority*/,
                                      PARSEC_DEV_ALL,
                                      PARSEC_DTD_EMPTY_FLAG,        &transA_g,
                                      PARSEC_DTD_EMPTY_FLAG,        &transB,
                                      PARSEC_DTD_EMPTY_FLAG,        &dcA.super.mb,
                                      PARSEC_DTD_EMPTY_FLAG,        &tempmm,
                                      PARSEC_DTD_EMPTY_FLAG,        &dcA.super.mb,
                                      PARSEC_DTD_EMPTY_FLAG,        &alpha_herk,
                                      PARSEC_DTD_EMPTY_FLAG,        PARSEC_DTD_TILE_OF(A, k, m),
                                      PARSEC_DTD_EMPTY_FLAG,        &ldak,
                                      PARSEC_DTD_EMPTY_FLAG,        PARSEC_DTD_TILE_OF(A, k, n),
                                      PARSEC_DTD_EMPTY_FLAG,        &ldak,
                                      PARSEC_DTD_EMPTY_FLAG,        &beta,
                                      PARSEC_DTD_EMPTY_FLAG,        PARSEC_DTD_TILE_OF(A, m, n),
                                      PARSEC_DTD_EMPTY_FLAG,        &ldan,
                                      PARSEC_DTD_ARG_END );
                }
                parsec_dtd_data_flush( dtd_tp, PARSEC_DTD_TILE_OF(A, k, m) );
            }
        }
    }

    parsec_dtd_data_flush_all( dtd_tp, (parsec_data_collection_t *)&dcA );

    /* finishing all the tasks inserted, but not finishing the handle */
    parsec_taskpool_wait( dtd_tp );

    /* Waiting on all handle and turning everything off for this context */
    parsec_context_wait( parsec );

    /* #### PaRSEC context is done #### */

    SYNC_TIME_PRINT(rank, ("\tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops\n",
                           P, Q, NB, N,
                           gflops=(flops/1e9)/sync_time_elapsed));

    parsec_dtd_task_class_release(dtd_tp, zpotrf_tc );
    parsec_dtd_task_class_release(dtd_tp, ztrsm_tc );
    parsec_dtd_task_class_release(dtd_tp, zherk_tc );
    parsec_dtd_task_class_release(dtd_tp, zgemm_tc );

    /* Cleaning up the parsec handle */
    parsec_taskpool_free( dtd_tp );
#if defined(DPLASMA_HAVE_CUDA)
    parsec_info_unregister(&parsec_per_device_infos, WoSI, NULL);
    free(infos);
#endif

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

    /* Cleaning data arrays we allocated for communication */
    dplasma_matrix_del2arena( tile_full );
    parsec_dtd_destroy_arena_datatype(parsec, TILE_FULL);
    parsec_dtd_data_collection_fini( (parsec_data_collection_t *)&dcA );

    parsec_data_free(dcA.mat); dcA.mat = NULL;
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);

    cleanup_parsec(parsec, iparam);
    return ret;
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

static void warmup_zpotrf(int rank, dplasma_enum_t uplo, int random_seed, parsec_context_t *parsec)
{
    int MB = 64;
    int NB = 64;
    int MT = 4;
    int NT = 4;
    int N = NB*NT;
    int M = MB*MT;
    int did;
    int info;

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                   rank, MB, NB, M, N, 0, 0,
                                   M, N, 1, 1, uplo));
    dcA.super.super.rank_of = always_local_rank_of;
    dcA.super.super.rank_of_key = always_local_rank_of_key;

    /* Do the CPU warmup first */
    dplasma_zplghe(parsec, (double)(N), uplo, &dcA.super, random_seed);
    parsec_taskpool_t *zpotrf = dplasma_zpotrf_New(uplo, &dcA.super, &info );
    zpotrf->devices_index_mask = 1<<0; /* Only CPU ! */
    parsec_context_add_taskpool(parsec, zpotrf);
    parsec_context_start(parsec);
    parsec_context_wait(parsec);
    dplasma_zpotrf_Destruct(zpotrf);

    /* Now do the other devices, skipping RECURSIVE */
    /* We know that there is a GPU-enabled version of this operation, so warm it up if some device is enabled */
    for(did = 2; did < (int)parsec_nb_devices; did++) {
        if(PARSEC_MATRIX_LOWER == uplo) {
            for(int i = 0; i < MT; i++) {
                for(int j = 0; j <= i; j++) {
                    parsec_data_t *dta = dcA.super.super.data_of(&dcA.super.super, i, j);
                    parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                }
            }
        } else {
            for(int i = 0; i < MT; i++) {
                for(int j = i; j < NT; j++) {
                    parsec_data_t *dta = dcA.super.super.data_of(&dcA.super.super, i, j);
                    parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                }
            }
        }
        dplasma_zplghe( parsec, (double)(N), uplo,
                        (parsec_tiled_matrix_t *)&dcA, random_seed);
        dplasma_zpotrf( parsec, uplo, &dcA.super );
        parsec_devices_release_memory();
    }

    parsec_data_free(dcA.mat); dcA.mat = NULL;
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA );
}
