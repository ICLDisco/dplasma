/*
 * Copyright (c) 2009-2010 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */


#include "dague.h"

#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

/* Plasma and math libs */
#include <cblas.h>
#include <plasma.h>
#include <lapacke.h>
#include <core_blas.h>

#include "scheduling.h"
#include "profiling.h"
#include "data_dist/matrix/two_dim_rectangle_cyclic/two_dim_rectangle_cyclic.h"
#include "dplasma.h"

#include "common.h"
#include "common_timing.h"

#define _FMULS(M, N, K) ( (DagDouble_t)(M) * (DagDouble_t)(N) * (DagDouble_t)(K) )
#define _FADDS(M, N, K) ( (DagDouble_t)(M) * (DagDouble_t)(N) * (DagDouble_t)(K) )

/**********************************
 * static functions
 **********************************/

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
static int check_solution(PLASMA_enum transA, PLASMA_enum transB,
                          Dague_Complex64_t alpha, two_dim_block_cyclic_t *ddescA, two_dim_block_cyclic_t *ddescB, 
                          Dague_Complex64_t beta, two_dim_block_cyclic_t *ddescC, two_dim_block_cyclic_t *ddescCfinal )
{
    int info_solution;
    double Anorm, Bnorm, Cinitnorm, Cdplasmanorm, Clapacknorm, Rnorm;
    double eps, result;
    double *work;
    int Am, An, Bm, Bn;
    Dague_Complex64_t mzone = (Dague_Complex64_t)-1.0;
    Dague_Complex64_t *A, *B, *Cinit, *Cfinal;

    int M   = ddescC->super.m;
    int N   = ddescC->super.n;
    int K   = ( transA == PlasmaNoTrans ) ? ddescA->super.n : ddescA->super.m ;
    int LDA = ddescA->super.lm;
    int LDB = ddescB->super.lm;
    int LDC = ddescC->super.lm;

    eps = LAPACKE_dlamch_work('e');

    if (transA == PlasmaNoTrans) {
        Am = M; An = K;
    } else {
        Am = K; An = M;
    }
    if (transB == PlasmaNoTrans) {
        Bm = K; Bn = N;
    } else {
        Bm = N; Bn = K;
    }

    work  = (double *)malloc(max(K,max(M, N))* sizeof(double));
    A     = (Dague_Complex64_t *)malloc((ddescA->super.mt)*(ddescA->super.nt)*(ddescA->super.bsiz)*sizeof(Dague_Complex64_t));
    B     = (Dague_Complex64_t *)malloc((ddescB->super.mt)*(ddescB->super.nt)*(ddescB->super.bsiz)*sizeof(Dague_Complex64_t));
    Cinit = (Dague_Complex64_t *)malloc((ddescC->super.mt)*(ddescC->super.nt)*(ddescC->super.bsiz)*sizeof(Dague_Complex64_t));
    Cfinal= (Dague_Complex64_t *)malloc((ddescC->super.mt)*(ddescC->super.nt)*(ddescC->super.bsiz)*sizeof(Dague_Complex64_t));

    twoDBC_to_lapack( ddescA, A,     LDA );
    twoDBC_to_lapack( ddescB, B,     LDB );
    twoDBC_to_lapack( ddescC,      Cinit,  LDC );
    twoDBC_to_lapack( ddescCfinal, Cfinal, LDC );

    Anorm        = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', Am, An, A,      LDA, work);
    Bnorm        = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', Bm, Bn, B,      LDB, work);
    Cinitnorm    = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', M,  N,  Cinit,  LDC, work);
    Cdplasmanorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', M,  N,  Cfinal, LDC, work);

    cblas_zgemm(CblasColMajor,
                (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                M, K, N, CBLAS_SADDR(alpha), A, LDA, B, LDB,
                CBLAS_SADDR(beta), Cinit, LDC);

    Clapacknorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', M, N, Cinit, LDC, work);

    cblas_zaxpy(LDC * N, CBLAS_SADDR(mzone), Cinit, 1, Cfinal, 1);
    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'i', M, N, Cfinal, LDC, work);

    if (getenv("DPLASMA_TESTING_VERBOSE"))
        printf("Rnorm %e, Anorm %e, Bnorm %e, Cinit %e, Cdplasmanorm %e, Clapacknorm %e\n",
               Rnorm, Anorm, Bnorm, Cinitnorm, Cdplasmanorm, Clapacknorm);
    
    result = Rnorm / ((Anorm + Bnorm + Cinitnorm) * max(M,N) * eps);
    if (  isinf(Clapacknorm) || isinf(Cdplasmanorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }

    free(work);
    free(A);
    free(B);
    free(Cinit);
    free(Cfinal);

    return info_solution;
}

int main(int argc, char ** argv)
{
    dague_context_t* dague;
    int iparam[IPARAM_SIZEOF];

    /* Set defaults for non argv iparams */
    iparam_default_gemm(iparam);
    /* Initialize DAGuE */
    dague = setup_dague(argc, argv, iparam);
    
    int rank  = iparam[IPARAM_RANK];
    int nodes = iparam[IPARAM_NNODES];
    int cores = iparam[IPARAM_NCORES];
    int gpus  = iparam[IPARAM_NGPUS];
    int prio  = iparam[IPARAM_PRIO];
    int P     = iparam[IPARAM_P];
    int Q     = iparam[IPARAM_Q];
    int M     = iparam[IPARAM_M];
    int N     = iparam[IPARAM_N];
    int K     = iparam[IPARAM_K];
    int LDA   = iparam[IPARAM_LDA];
    int LDB   = iparam[IPARAM_LDB];
    int LDC   = iparam[IPARAM_LDC];
    int MB    = iparam[IPARAM_MB];
    int NB    = iparam[IPARAM_NB];
    int SMB   = iparam[IPARAM_SMB];
    int SNB   = iparam[IPARAM_SNB];
    int check = iparam[IPARAM_CHECK];
    int loud  = iparam[IPARAM_VERBOSE];
    
    int tA    = PlasmaNoTrans;
    int tB    = PlasmaNoTrans;
    Dague_Complex64_t alpha =  0.51;
    Dague_Complex64_t beta  = -0.42;

    DagDouble_t flops, gflops;
#if defined(PRECISIONS_z) || defined(PRECISIONS_c)
    flops = 2.*_FADDS(M, N, K) + 6.*_FMULS(M, N, K);
#else
    flops = _FADDS(M, N, K) + _FMULS(M, N, K);
#endif

    /* initializing matrix structure */
    two_dim_block_cyclic_t ddescA;
    two_dim_block_cyclic_init(&ddescA, matrix_ComplexDouble, nodes, cores, rank, MB, NB, M, K, 0, 0, LDA, K, SMB, SNB, P);
    ddescA.mat = dague_data_allocate((size_t)ddescA.super.nb_local_tiles * (size_t)ddescA.super.bsiz * (size_t)ddescA.super.mtype);
    
    two_dim_block_cyclic_t ddescB;
    two_dim_block_cyclic_init(&ddescB, matrix_ComplexDouble, nodes, cores, rank, MB, NB, K, N, 0, 0, LDB, N, SMB, SNB, P);
    ddescB.mat = dague_data_allocate((size_t)ddescB.super.nb_local_tiles * (size_t)ddescB.super.bsiz * (size_t)ddescB.super.mtype);
    
    two_dim_block_cyclic_t ddescC;
    two_dim_block_cyclic_init(&ddescC, matrix_ComplexDouble, nodes, cores, rank, MB, NB, M, N, 0, 0, LDC, N, SMB, SNB, P);
    ddescC.mat = dague_data_allocate((size_t)ddescC.super.nb_local_tiles * (size_t)ddescC.super.bsiz * (size_t)ddescC.super.mtype);

    two_dim_block_cyclic_t ddescC2;
    if(check) {
        two_dim_block_cyclic_init(&ddescC2, matrix_ComplexDouble, nodes, cores, rank, MB, NB, M, N, 0, 0, LDC, N, SMB, SNB, P);
        ddescC2.mat = dague_data_allocate((size_t)ddescC2.super.nb_local_tiles * (size_t)ddescC2.super.bsiz * (size_t)ddescC2.super.mtype);
    }
    
    if(!check) 
    {
        /* matrix generation */
        if(loud) printf("Generate matrices ... ");
        generate_tiled_random_mat((tiled_matrix_desc_t *) &ddescA, 100);
        generate_tiled_random_mat((tiled_matrix_desc_t *) &ddescB, 200);
        generate_tiled_random_mat((tiled_matrix_desc_t *) &ddescC, 300);
        if(loud) printf("Done\n");

        /* Create GEMM DAGuE */
        if(loud) printf("Generate GEMM DAG ... ");
        TIME_START();
        dague_object_t* dague_gemm = 
            dplasma_zgemm_New(tA, tB, 
                              (Dague_Complex64_t)alpha,
                              (tiled_matrix_desc_t *)&ddescA, 
                              (tiled_matrix_desc_t *)&ddescB,
                              (Dague_Complex64_t)beta,
                              (tiled_matrix_desc_t *)&ddescC);
        dague_enqueue(dague, dague_gemm);
        if(loud) printf("Done\n");
        if(loud) TIME_PRINT(("DAG creation: %u total tasks enqueued\n", dague->taskstodo));

        /* lets rock! */
        SYNC_TIME_START();
        TIME_START();
        dague_progress(dague);
        if(loud) TIME_PRINT(("Dague proc %d:\tcomputed %u tasks,\t%f task/s\n",
                    rank, dague_gemm->nb_local_tasks,
                    dague_gemm->nb_local_tasks/time_elapsed));
        SYNC_TIME_PRINT(("Dague progress:\t%d %d %f gflops\n", N, NB,
                         gflops = (flops/1e9)/(sync_time_elapsed)));
    }
    else
    { 
/* Iterate on the transpose forms. TODO: LDB is set incorrecly for T and H */
#if defined(PRECISIONS_z) || defined(PRECISIONS_c)
        for(tA=0; tA<3; tA++) {
            for(tB=0; tB<3; tB++) {
#else
        for(tA=0; tA<2; tA++) {
            for(tB=0; tB<2; tB++) {
#endif
                printf("***************************************************\n");
                printf(" ----- TESTING DGEMM (%s, %s) -------- \n",
                       transstr[tA], transstr[tB]);
                
                /* matrix generation */
                if(loud) printf("Generate matrices ... ");
                generate_tiled_random_mat((tiled_matrix_desc_t *) &ddescA,  100);
                generate_tiled_random_mat((tiled_matrix_desc_t *) &ddescB,  200);
                generate_tiled_random_mat((tiled_matrix_desc_t *) &ddescC,  300);
                generate_tiled_random_mat((tiled_matrix_desc_t *) &ddescC2, 300);
                if(loud) printf("Done\n");

                /* Create GEMM DAGuE */
                if(loud) printf("Compute ... ... ");
                dplasma_dgemm(dague, trans[tA], trans[tB],
                              (Dague_Complex64_t)alpha,
                              (tiled_matrix_desc_t *)&ddescA, 
                              (tiled_matrix_desc_t *)&ddescB, 
                              (Dague_Complex64_t)beta, 
                              (tiled_matrix_desc_t *)&ddescC);
                if(loud) printf("Done\n");
                
                /* Check the solution */
                int info_solution = check_solution( trans[tA], trans[tB], 
                                                    alpha, &ddescA,  &ddescB, 
                                                    beta,  &ddescC2, &ddescC);
                if (info_solution == 0) {
                    printf(" ---- TESTING DGEMM (%s, %s) ...... PASSED !\n",
                           transstr[tA], transstr[tB]);
                }
                else {
                    printf(" ---- TESTING DGEMM (%s, %s) ... FAILED !\n",
                           transstr[tA], transstr[tB]);
                }
                printf("***************************************************\n");
            }
        }

        dague_data_free(ddescC2.mat);
    }

    dague_data_free(ddescA.mat);
    dague_data_free(ddescB.mat);
    dague_data_free(ddescC.mat);

    cleanup_dague(dague);
    return 0;
}