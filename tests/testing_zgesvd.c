/*
 * Copyright (c) 2011-2023 The University of Tennessee and The University
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

static int check_solution(int N, const double *E1, const double *E2);
static void warmup_zgesvd(int rank, int random_seed, parsec_context_t *parsec);

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF];
    int ret = 0;
    double *s0 = NULL;
    double *s1 = NULL;
    double *e  = NULL;
    int minMN;
    int info_solution;
    double time_ge2gb, time_gb2bd, time_solve = -1.;
    int rc;

    /* Ensure BLAS are sequential and set thread affinity for the master */
/* #if defined(__ICC) || defined(__INTEL_COMPILER) */
/*     kmp_set_defaults("KMP_AFFINITY=disabled"); */
/*     mkl_set_num_threads( 1 ); */
/* #endif */

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 32, 200, 200);
    iparam[IPARAM_KP] = 1;
    iparam[IPARAM_KQ] = 1;
    iparam[IPARAM_LDA] = -'m';
    iparam[IPARAM_LDB] = -'m';

    /* Initialize Parsec */
    parsec = setup_parsec(argc, argv, iparam);

    /* Make sure KP and KQ are set to 1, since it conflicts with HQR */
    iparam[IPARAM_KP] = 1;
    iparam[IPARAM_KQ] = 1;

    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_FLOPS(FLOPS_ZGEBRD, ((DagDouble_t)M, (DagDouble_t)N));

    LDA = max(M, LDA);

    warmup_zgesvd(rank, random_seed, parsec);

    if ( M < N ) {
        fprintf(stderr, "This testing can only perform SVD on matrices with M >= N\n");
        return EXIT_FAILURE;
    }
    minMN = dplasma_imin(M, N);

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, LDA, N, 0, 0,
                               M, N, P, nodes/P, KP, KQ, IP, JQ));
    PASTE_CODE_ALLOCATE_MATRIX(dcBand, 1,
        parsec_matrix_block_cyclic, (&dcBand, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB+1, NB, MB+1, minMN, 0, 0,
                               MB+1, minMN, 1, 1, 1, 1, IP, JQ));

    s1 = (double*)malloc( minMN * sizeof(double));
    e  = (double*)malloc( minMN * sizeof(double));

    for(int t = 0; t < nruns; t++) {
        /* Initialize the matrix */
        if(loud > 3) printf("+++ Generate matrices ... ");

        /* Generate the matrix on rank 0 */
        if ( check ) {

            /* Generate the singular values vector as in latms routines for check purpose */
            if (rank == 0 && NULL ==s0 )
            {
                double tmp = 1. / (double)N;
                double alp = ( 1. - tmp ) / ((double)( N - 1 ));
                int i;
                s0 = (double *) malloc(minMN * sizeof(double));

                s0[0] = 1.;
                for(i=1; i < minMN; i++){
                    s0[i] = (double)(N-i-1) * alp + tmp;
                }
            }

            dplasma_zlatms( parsec, dplasmaGeneral, (double)N, (parsec_tiled_matrix_t *)&dcA, 3872);
        }
        else {
            dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, 3872);
        }

        /* Create Parsec */
        PASTE_CODE_ENQUEUE_KERNEL(parsec, zgebrd_ge2gb,
                                (IB,
                                (parsec_tiled_matrix_t*)&dcA,
                                (parsec_tiled_matrix_t*)&dcBand));

        /* lets rock! */
        SYNC_TIME_START();
        rc = parsec_context_start(parsec);
        PARSEC_CHECK_ERROR(rc, "parsec_context_start");
        TIME_START();
        rc = parsec_context_wait(parsec);
        PARSEC_CHECK_ERROR(rc, "parsec_context_wait");
        SYNC_TIME_STOP();
        time_ge2gb = sync_time_elapsed;

        if( rank == 0 ) {
    /* #if defined(__ICC) || defined(__INTEL_COMPILER) */
    /*         mkl_set_num_threads( iparam[IPARAM_NCORES] ); */
    /* #endif */
            /* Reduce the band */
            TIME_START();
            info_solution = LAPACKE_zgbbrd( LAPACK_COL_MAJOR,
                                            'N',
                                            M, N,
                                            0, 0, NB,
                                            dcBand.mat, MB+1,
                                            s1, e,
                                            NULL, 1,
                                            NULL, 1,
                                            NULL, 1 );
            TIME_STOP();
            time_gb2bd = time_elapsed;

            /* Solve the bidiagonal SVD problem */
            if (info_solution == 0){
                TIME_START();
                info_solution = LAPACKE_zbdsqr( LAPACK_COL_MAJOR, 'U',
                                                minMN, 0, 0, 0,
                                                s1, e,
                                                NULL, 1, NULL, 1, NULL, 1 );
                TIME_STOP();
                time_solve = time_elapsed;
            }

    /* #if defined(__ICC) || defined(__INTEL_COMPILER) */
    /*         mkl_set_num_threads( 1 ); */
    /* #endif */
            fprintf(stderr, "WARNING: This code is using the non optimized Lapack zbdsqr subroutine to reduce the band to bi-diagonal form. Please replace this call by the multi-threaded PLASMA implementation in order to get performance\n");
            printf("zgeqrf GESVD computation NP= %d NC= %d P= %d IB= %d MB= %d NB= %d qr_a= %d qr_p = %d treel= %d treeh= %d domino= %d R-bidiag= %d M= %d N= %d : %e %e %e / %f gflops\n",
                iparam[IPARAM_NNODES],
                iparam[IPARAM_NCORES],
                iparam[IPARAM_P],
                iparam[IPARAM_IB],
                iparam[IPARAM_MB],
                iparam[IPARAM_NB],
                iparam[IPARAM_QR_TS_SZE],
                iparam[IPARAM_QR_HLVL_SZE],
                iparam[IPARAM_LOWLVL_TREE],
                iparam[IPARAM_HIGHLVL_TREE],
                iparam[IPARAM_QR_DOMINO],
                iparam[IPARAM_QR_TSRR],
                iparam[IPARAM_M],
                iparam[IPARAM_N],
                time_ge2gb, time_gb2bd, time_solve,
                gflops = (flops/1e9)/(time_ge2gb+time_gb2bd+time_solve));

    #if defined(PARSEC_SIM)
            printf("zgeqrf GESVD simulation NP= %d NC= %d P= %d qr_a= %d qr_p = %d treel= %d treeh= %d domino= %d RR= %d MT= %d NT= %d : %d \n",
                iparam[IPARAM_NNODES],
                iparam[IPARAM_NCORES],
                iparam[IPARAM_P],
                iparam[IPARAM_QR_TS_SZE],
                iparam[IPARAM_QR_HLVL_SZE],
                iparam[IPARAM_LOWLVL_TREE],
                iparam[IPARAM_HIGHLVL_TREE],
                iparam[IPARAM_QR_DOMINO],
                iparam[IPARAM_QR_TSRR],
                MT, NT,
                parsec_getsimulationdate( parsec ));
    #endif
        }

        dplasma_zgebrd_ge2gb_Destruct( PARSEC_zgebrd_ge2gb );
    }

    if( check && (rank==0) ) {
        if (info_solution == 0 ) {
            info_solution = check_solution(minMN, s0, s1);
        }

        if (info_solution == 0) {
            printf("***************************************************\n"
                   " ---- TESTING ZGESVD .. M >= N ........... PASSED !\n"
                   "***************************************************\n");
        }
        else {
            printf("***************************************************\n"
                   " ---- TESTING ZGESVD .. M >= N .. FAILED !\n"
                   "***************************************************\n");
        }
    }

    free(s1);
    free(s0);
    free(e);

    parsec_data_free(dcA.mat);
    parsec_data_free(dcBand.mat);

    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcBand);

    cleanup_parsec(parsec, iparam);

    return ret;
}

static double dplasma_dmax( double a, double b ) {
    if ( a < b ) {
        return b;
    } else {
        return a;
    }
}

/*------------------------------------------------------------
 *  Check the eigenvalues
 */
static int check_solution(int N, const double *E1, const double *E2)
{
    int info_solution, i;
    double resid;
    double maxtmp;
    double maxel = fabs( fabs(E1[0]) - fabs(E2[0]) );
    double maxeig = dplasma_dmax( fabs(E1[0]), fabs(E2[0]) );
    double eps = LAPACKE_dlamch_work('e');

    for (i = 1; i < N; i++){
        resid   = fabs(fabs(E1[i])-fabs(E2[i]));
        maxtmp  = dplasma_dmax(fabs(E1[i]), fabs(E2[i]));

        /* Update */
        maxeig = dplasma_dmax(maxtmp, maxeig);
        maxel  = dplasma_dmax(resid,  maxel );
    }

    maxel = maxel / (maxeig * N * eps);
    printf(" ======================================================\n");
    printf(" | S - singularcomputed | / (|S| * N * eps) : %e \n",  maxel );
    printf(" ======================================================\n");

    if ( isnan(maxel) || isinf(maxel) || (maxel > 100) ) {
        printf("-- The singular values are suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The singular values are CORRECT ! \n");
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

static void warmup_zgesvd(int rank, int random_seed, parsec_context_t *parsec)
{
    int MB = 64;
    int IB = 40;
    int NB = 64;
    int MT = 4;
    int NT = 4;
    int N = NB*NT;
    int M = MB*MT;
    double *s1, *e;

    /* initializing matrix structure */
    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               rank, MB, NB, M, N, 0, 0,
                               M, N, 1, 1, 1, 1, 0, 0));
    dcA.super.super.rank_of = always_local_rank_of;
    dcA.super.super.rank_of_key = always_local_rank_of_key;
    PASTE_CODE_ALLOCATE_MATRIX(dcBand, 1,
        parsec_matrix_block_cyclic, (&dcBand, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_LAPACK,
                               rank, MB+1, NB, MB+1, M, 0, 0,
                               MB+1, M, 1, 1, 1, 1, 0, 0));
    dcBand.super.super.rank_of = always_local_rank_of;
    dcBand.super.super.rank_of_key = always_local_rank_of_key;
    s1 = (double*)malloc( M * sizeof(double));
    e  = (double*)malloc( M * sizeof(double));

    /* Do the CPU warmup first */
    dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, random_seed);
    parsec_taskpool_t *zgesvd = dplasma_zgebrd_ge2gb_New(IB,
                               (parsec_tiled_matrix_t*)&dcA,
                               (parsec_tiled_matrix_t*)&dcBand);
    zgesvd->devices_index_mask = 1<<0; /* Only CPU ! */
    parsec_context_add_taskpool(parsec, zgesvd);
    parsec_context_start(parsec);
    parsec_context_wait(parsec);
    (void)LAPACKE_zgbbrd( LAPACK_COL_MAJOR,
                                        'N',
                                        M, N,
                                        0, 0, NB,
                                        dcBand.mat, MB+1,
                                        s1, e,
                                        NULL, 1,
                                        NULL, 1,
                                        NULL, 1 );
    (void)LAPACKE_zbdsqr( LAPACK_COL_MAJOR, 'U',
                                            M, 0, 0, 0,
                                            s1, e,
                                            NULL, 1, NULL, 1, NULL, 1 );

    /* Check for which device type (skipping RECURSIVE), we need to warmup this operation */
    for(int dtype = PARSEC_DEV_RECURSIVE+1; dtype < PARSEC_DEV_MAX_NB_TYPE; dtype++) {
        for(int i = 0; i < (int)zgesvd->nb_task_classes; i++) {
            for(int j = 0; NULL != zgesvd->task_classes_array[i]->incarnations[j].hook; j++) {
                if( zgesvd->task_classes_array[i]->incarnations[j].type & dtype ) {
                    goto do_run; /* We found one class that was on that device, no need to try more incarnations or task classes */
                }
            }
        }
        continue; /* No incarnation of this device type on any task class; try another type */
    do_run:
        for(int did = 0; did < (int)parsec_nb_devices; did++) {
            parsec_device_module_t *dev = parsec_mca_device_get(did);
            if( !(dev->type & dtype) )
                continue;
            /* This should work, right? Unfortunately, we can't test until there is a <dev>-enabled implementation for this test */
            for(int m = 0; m < MT; m++) {
                for(int n = 0; n < NT; n++) {
                    parsec_data_t *dta = dcA.super.super.data_of(&dcA.super.super, m, n);
                    parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                    dta = dcA.super.super.data_of(&dcA.super.super, m, n);
                    parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                    if(m == 0) {
                        dta = dcBand.super.super.data_of(&dcBand.super.super, m, n);
                        parsec_advise_data_on_device( dta, did, PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
                    }
                }
            }
            dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, random_seed);
            parsec_taskpool_t *zgesvd_device = dplasma_zgebrd_ge2gb_New(IB,
                                    (parsec_tiled_matrix_t*)&dcA,
                                    (parsec_tiled_matrix_t*)&dcBand);
            zgesvd->devices_index_mask = 1<<0; /* Only CPU ! */
            parsec_context_add_taskpool(parsec, zgesvd_device);
            parsec_context_start(parsec);
            parsec_context_wait(parsec);
            dplasma_zgebrd_ge2gb_Destruct( zgesvd_device );
            /* No need to redo zgbbrd and zbdsqr as those are LAPACK / CPU-only */
        }
    }

    free(e);
    free(s1);
    dplasma_zgebrd_ge2gb_Destruct( zgesvd );

}
