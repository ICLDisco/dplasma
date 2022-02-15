/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "common.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/data_internal.h"

/* Including the bulge chassing */
#define FADDS_ZHBRDT(__n) (-1)
#define FMULS_ZHBRDT(__n) (-1)

int main(int argc, char *argv[])
{
    parsec_context_t *parsec;
    int iparam[IPARAM_SIZEOF];

    /* Set defaults for non argv iparams */
    iparam_default_facto(iparam);
    iparam_default_ibnbmb(iparam, 48, 144, 144);
    iparam[IPARAM_NGPUS] = DPLASMA_ERR_NOT_SUPPORTED;

    parsec = setup_parsec(argc, argv, iparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    if(P != 1)
        fprintf(stderr, "!!! This algorithm works on a band 1D matrix. The value of P=%d has been overriden, the actual grid is %dx%d\n", P, 1, nodes);

    PASTE_CODE_FLOPS_COUNT(FADDS_ZHBRDT, FMULS_ZHBRDT, ((DagDouble_t)N));

    /*
     PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
     parsec_matrix_sym_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE,
     rank, MB, NB, LDA, N, 0, 0,
     N, N, P, nodes/P, MatrixLower))
     */

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
                               parsec_matrix_block_cyclic, (&dcA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank, MB+1, NB+2, MB+1, (NB+2)*NT, 0, 0,
                                                      MB+1, (NB+2)*NT, 1, nodes, 1, KQ, IP, JQ /* 1D cyclic */ ));

    double time_avg = 0.0;

    for(int t = 0; t < nruns+1; t++) {
        if(loud > 2) printf("+++ Generate matrices ... ");
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcA, 3872);
        if(loud > 2) printf("Done\n");

        PASTE_CODE_ENQUEUE_KERNEL(parsec, zhbrdt,
                                  ((parsec_tiled_matrix_t*)&dcA));

        PASTE_CODE_PROGRESS_KERNEL(parsec, zhbrdt, t);
        if(t > 0) {
            time_avg += sync_time_elapsed/nruns;
        }
        dplasma_zhbrdt_Destruct( PARSEC_zhbrdt );
    }
    PASTE_PROF_INFO;
    if(loud >= 5 && rank == 0) {
        printf("<DartMeasurement name=\"duration\" type=\"numeric/double\"\n"
               "                 encoding=\"none\" compression=\"none\">\n"
               "%g\n"
               "</DartMeasurement>\n",
               time_avg);
    }

    if( check ) {
        printf( "No check implemented yet.\n" );

#if defined(PARSEC_HAVE_MPI)
        /* Regenerate A, distributed so that the random generators are doing
         * the same things */
        PASTE_CODE_ALLOCATE_MATRIX(dcAcpy, 1,
                                   parsec_matrix_block_cyclic, (&dcAcpy, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                          rank, MB+1, NB+2, MB+1, (NB+2)*NT,
                                                          0, 0, MB+1, (NB+2)*NT, 1, nodes, 1, KQ, IP, JQ));
        dplasma_zplrnt( parsec, 0, (parsec_tiled_matrix_t *)&dcAcpy, 3872);

        /* Gather Acpy on rank 0 */
        PASTE_CODE_ALLOCATE_MATRIX(dcLAcpy, 1,
                                   parsec_matrix_block_cyclic, (&dcLAcpy, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                          rank, MB+1, NB+2, MB+1, (NB+2)*NT,
                                                          0, 0, MB+1, (NB+2)*NT, 1, 1, 1, 1, IP, JQ));

        /* Gather A diagonal and subdiagonal on rank 0 */
        PASTE_CODE_ALLOCATE_MATRIX(dcLA, 1,
                                   parsec_matrix_block_cyclic, (&dcLA, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                                                          rank, 2, NB, 2, NB*NT,
                                                          0, 0, 2, NB*NT, 1, 1, 1, 1, IP, JQ));
        if(rank == 0) {
            for(int t = 0; t < NT; t++)
            {
                int rsrc = dcA.super.super.rank_of(&dcA.super.super, 0, t);
                if(rsrc == 0)
                {
                    dplasma_complex64_t* datain = parsec_data_copy_get_ptr(parsec_data_get_copy(dcA.super.super.data_of(&dcA.super.super, 0, t), 0));
                    dplasma_complex64_t* dataout = parsec_data_copy_get_ptr(parsec_data_get_copy(dcLA.super.super.data_of(&dcA.super.super, 0, t), 0));
                    for(int n = 0; n < NB; n++)
                        for(int m = 0; m < 2; m++)
                            {
                                dataout[m+n*2] = datain[m+n*(MB+1)];
                            }
                }
                else
                {
                    dplasma_complex64_t* dataout = parsec_data_copy_get_ptr(parsec_data_get_copy(dcLA.super.super.data_of(&dcLA.super.super, 0, t), 0));
                    MPI_Recv(dataout, 2*NB, parsec_datatype_double_complex_t, rsrc, t, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        else
        {
            MPI_Datatype bidiagband_dtt;
            MPI_Type_vector(NB, 2, MB+1, parsec_datatype_double_complex_t, &bidiagband_dtt);

            for(int t = 0; t < NT; t++) {
                if(dcA.super.super.rank_of(&dcA.super.super, 0, t) == (uint32_t)rank)
                {
                    dplasma_complex64_t* datain = parsec_data_copy_get_ptr(parsec_data_get_copy(dcA.super.super.data_of(&dcA.super.super, 0, t), 0));
                    MPI_Send(datain, 1, bidiagband_dtt, 0, t, MPI_COMM_WORLD);
                }
            }
        }
#endif  /* defined(PARSEC_HAVE_MPI) */
    }

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA);

    cleanup_parsec(parsec, iparam);

    return EXIT_SUCCESS;
}

/*--------------------------------------------------------------
 * Check the solution
 */
#if 0
static int check_solution(int N, double *E1, double *E2, double eps)
{
    int info_solution, i;
    double *Residual = (double *)malloc(N*sizeof(double));
    double maxtmp;
    double maxel = fabs(fabs(E1[0])-fabs(E2[0]));
    double maxeig = dplasma_fmax(fabs(E1[0]), fabs(E2[0]));
    for (i = 1; i < N; i++){
        Residual[i] = fabs(fabs(E1[i])-fabs(E2[i]));
        maxtmp      = dplasma_fmax(fabs(E1[i]), fabs(E2[i]));
        maxeig      = dplasma_fmax(maxtmp, maxeig);
        //printf("Residu: %f E1: %f E2: %f\n", Residual[i], E1[i], E2[i] );
        if (maxel < Residual[i])
            maxel =  Residual[i];
    }

    //printf("maxel: %.16f maxeig: %.16f \n", maxel, maxeig );

    printf(" ======================================================\n");
    printf(" | D -  eigcomputed | / (|D| ulp)      : %15.3E \n",  maxel/(maxeig*eps) );
    printf(" ======================================================\n");


    printf("============\n");
    printf("Checking the eigenvalues of A\n");
    if (isnan(maxel / eps) || isinf(maxel / eps) || ((maxel / (maxeig*eps)) > 1000.0) ) {
        //printf("isnan: %d %f %e\n", isnan(maxel / eps), maxel, eps );
        //printf("isinf: %d %f %e\n", isinf(maxel / eps), maxel, eps );
        printf("-- The eigenvalues are suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The eigenvalues are CORRECT ! \n");
        info_solution = 0;
    }

    free(Residual);
    return info_solution;
}
#endif
