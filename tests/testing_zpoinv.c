/*
 * Copyright (c) 2009-2023 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "flops.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

static void warmup_zpoinv(int rank, dplasma_enum_t uplo, int random_seed, parsec_context_t *parsec);

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    dplasma_enum_t uplo = dplasmaLower;
    int info = 0;
    int ret = 0;

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 0, 180, 180);
    iparam[IPARAM_NGPUS] = 0;

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZPOTRF, ((DagDouble_t)N));
    flops += FLOPS_ZPOTRI((DagDouble_t)N);

    /* initializing matrix structure */
    LDA = dplasma_imax( LDA, N );
    KP = 1;
    KQ = 1;

    warmup_zpoinv(rank, uplo, random_seed, parsec);

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
                                   rank, MB, NB, LDA, N, 0, 0,
                                   N, N, P, nodes/P, uplo));

    for(int t = 0; t < nruns; t++) {
        /* matrix generation */
        if(loud > 3) printf("+++ Generate matrices ... ");
        dplasma_zplghe( parsec, (double)(N), uplo,
                        (parsec_tiled_matrix_t *)&dcA, random_seed);
        if(loud > 3) printf("Done\n");

        if (async) {
            PASTE_CODE_ENQUEUE_KERNEL(parsec, zpoinv,
                                    (uplo, (parsec_tiled_matrix_t*)&dcA, &info));
            PASTE_CODE_PROGRESS_KERNEL(parsec, zpoinv);
            dplasma_zpoinv_Destruct( PARSEC_zpoinv );
        }
        else {
            SYNC_TIME_START();
            info = dplasma_zpoinv_sync( parsec, uplo, (parsec_tiled_matrix_t*)&dcA );
            SYNC_TIME_PRINT(rank, ("zpoinv\tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops\n",
                                P, Q, NB, N,
                                gflops=(flops/1e9)/sync_time_elapsed));
        }
    }

    if( 0 == rank && info != 0 ) {
        printf("-- Factorization is suspicious (info = %d) ! \n", info);
        ret |= 1;
    }
    if( !info && check ) {
        /* Check the factorization */
        PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
            parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE, rank,
                                   MB, NB, LDA, N, 0, 0, N, N, P, nodes/P, 1, 1, IP, JQ));
        dplasma_zplghe( parsec, (double)(N), dplasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcA0, random_seed);

        ret |= check_zpoinv( parsec, (rank == 0) ? loud : 0, uplo,
                             (parsec_tiled_matrix_t *)&dcA0,
                             (parsec_tiled_matrix_t *)&dcA );

        if (ret) {
            printf("-- Innversion is suspicious ! \n");
        }
        else
        {
            printf("-- Inversion is CORRECT ! \n");
        }

        /* Cleanup */
        parsec_data_free(dcA0.mat); dcA0.mat = NULL;
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0 );
    }

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

static void warmup_zpoinv(int rank, dplasma_enum_t uplo, int random_seed, parsec_context_t *parsec)
{
    int MB = 64;
    int NB = 64;
    int MT = 4;
    int NT = 4;
    int M = MB*MT;
    int N = NB*NT;
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
    parsec_taskpool_t *zpoinv = dplasma_zpoinv_New(uplo, &dcA.super, &info );
    zpoinv->devices_index_mask = 1<<0; /* Only CPU ! */
    parsec_context_add_taskpool(parsec, zpoinv);
    parsec_context_start(parsec);
    parsec_context_wait(parsec);
    dplasma_zpoinv_Destruct(zpoinv);

    /* Now do the other devices, skipping RECURSIVE */
    /* We know that there is a GPU-enabled version of this operation, so warm it up if some device is enabled */
    for(did = 1; did < (int)parsec_nb_devices; did++) {
        if(!parsec_mca_device_is_gpu(did)) continue;
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
        dplasma_zplghe(parsec, (double)(N), uplo, &dcA.super, random_seed);
        dplasma_zpoinv( parsec, uplo, &dcA.super );
        parsec_devices_release_memory();
    }

    parsec_data_free(dcA.mat); dcA.mat = NULL;
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA );
}
