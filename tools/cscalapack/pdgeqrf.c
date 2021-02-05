/*
 * Copyright (c) 2009-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2010      University of Denver, Colorado.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include "myscalapack.h"
#include "common.h"

static double check_solution( int params[], double *Alu, double *tau );

int main( int argc, char **argv ) {
    int params[PARAMS_SIZE];
    int info;
    int ictxt, nprow, npcol, myrow, mycol, iam;
    int number_runs;
    int m, n, nb, mb, s, mloc, nloc, verif, iseed;
    int descA[9];
    double *A = NULL;
    double *tau = NULL;
    double resid, telapsed, gflops, pgflops;

    setup_params( params, argc, argv );
    ictxt = params[PARAM_BLACS_CTX];
    iam   = params[PARAM_RANK];
    m     = params[PARAM_M];
    n     = params[PARAM_N];
    mb    = params[PARAM_MB];
    nb    = params[PARAM_NB];
    s     = params[PARAM_NRHS];
    iseed = params[PARAM_SEED];
    verif = params[PARAM_VALIDATE];
    number_runs = params[PARAM_NRUNS];

#ifdef DPLASMA_WRAPPER_ON
    parsec_init_wrapper_();
#endif

    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
    mloc = numroc_( &m, &mb, &myrow, &i0, &nprow );
    nloc = numroc_( &n, &nb, &mycol, &i0, &npcol );
    descinit_( descA, &m, &n, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );

    A = malloc( sizeof(double)*mloc*nloc );
    tau = malloc( sizeof(double)*imin(m,n) );

    double *work=NULL; int lwork=-1; double getlwork;
    double t1, t2;
#ifdef DPLASMA_WRAPPER_ON
    pdlatsqr_( &m, &n, A, &i1, &i1, descA, tau, &getlwork, &lwork, &info );
#else
    pdgeqrf_( &m, &n, A, &i1, &i1, descA, tau, &getlwork, &lwork, &info );
#endif
    assert( 0 == info );
    lwork = (int)getlwork;
    work = malloc( sizeof(double)*lwork );

    int t;
    for(t = 0; t < number_runs; t++) {
#ifdef DPLASMA_WRAPPER_ON
        parsec_wrapper_devices_release_memory_();
#endif
        scalapack_pdplrnt( A,
                           m, n,
                           mb, nb,
                           myrow, mycol,
                           nprow, npcol,
                           mloc,
                           iseed );
        t1 = MPI_Wtime();
#ifdef DPLASMA_WRAPPER_ON
        pdlatsqr_( &m, &n, A, &i1, &i1, descA, tau, work, &lwork, &info );
#else
        pdgeqrf_( &m, &n, A, &i1, &i1, descA, tau, work, &lwork, &info );
#endif
        assert( 0 == info );
        t2 = MPI_Wtime();
        telapsed = t2-t1;

        if( 0 != iam ) {
            MPI_Reduce( &telapsed, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
        }
        else {
            MPI_Reduce( MPI_IN_PLACE, &telapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
            gflops = FLOPS_DGEQRF((double)m, (double)n)/1e+9/telapsed;
            pgflops = gflops/(((double)nprow)*((double)npcol));
        }

        if( 0 == iam ) {
            printf("[****] TIMEHL(s) %12.5f : PDGEQRF \tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops"
                   " - ENQ&PROG&DEST %12.5f : %14f gflops"
                   " - ENQ %12.5f - DEST %12.5f\n",
                   telapsed, nprow, npcol, nb, n,
                   gflops,
                   telapsed,
                   gflops,
                   0.0, 0.0);
        }
#ifdef DPLASMA_WRAPPER_ON
        parsec_wrapper_devices_reset_load_();
#endif
    }

    if ( verif ) {
        resid = check_solution( params, A, tau );
    } else {
        resid = -1;
    }
    if( 0 == iam ) {
        printf( "### PDGEQRF ###\n"
                "#%4sx%-4s %7s %7s %4s %4s # %10s \n", "P", "Q", "M", "N", "NB", "NRHS", "resid");
        printf( " %4d %-4d %7d %7d %4d %4d   %10.3e \n", nprow, npcol, m, n, nb, s, resid );
    }
#ifdef DPLASMA_WRAPPER_ON
    parsec_fini_wrapper_();
#endif

    free( A ); A = NULL;
    free( tau ); tau = NULL;
    free(work); work = NULL;

    Cblacs_exit( 0 );
    return 0;
}

static double check_solution( int params[], double* Aqr, double *tau ) {
    double resid = NAN;
    int info;
    int ictxt = params[PARAM_BLACS_CTX];
    int iam   = params[PARAM_RANK];
    int m     = params[PARAM_M];
    int n     = params[PARAM_N];
    int mb    = params[PARAM_MB];
    int nb    = params[PARAM_NB];
    int s     = params[PARAM_NRHS];
    int iseed = params[PARAM_SEED];
    int nprow, npcol, myrow, mycol;
    int mloc, nloc, sloc;
    double *A=NULL; int descA[9];
    double *B=NULL; int descB[9];
    double *X=NULL;
    double eps, Anorm, Bnorm, Xnorm, Rnorm;

    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
    mloc = numroc_( &m, &mb, &myrow, &i0, &nprow );
    nloc = numroc_( &n, &nb, &mycol, &i0, &npcol );
    sloc = numroc_( &s, &nb, &mycol, &i0, &npcol );

    /* recreate A */
    descinit_( descA, &m, &n, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
    assert( 0 == info );
    A = malloc( sizeof(double)*mloc*nloc );
    scalapack_pdplrnt( A,
                       m, n,
                       mb, nb,
                       myrow, mycol,
                       nprow, npcol,
                       mloc,
                       iseed );

    /* create B and copy it to X */
    descinit_( descB, &n, &s, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
    assert( 0 == info );
    B = malloc( sizeof(double)*mloc*sloc );
    X = malloc( sizeof(double)*mloc*sloc );

    scalapack_pdplrnt( B,
                       n, s,
                       mb, nb,
                       myrow, mycol,
                       nprow, npcol,
                       mloc,
                       iseed + 1 );

    pdlacpy_( "All", &n, &s, B, &i1, &i1, descB, X, &i1, &i1, descB );
    {
        double *work = malloc( sizeof(double)*mloc );
        Anorm = pdlange_( "I", &n, &n, A, &i1, &i1, descA, work );
        Bnorm = pdlange_( "I", &n, &s, B, &i1, &i1, descB, work );
        free( work );
    }

    /* Compute X from Aqr */
    {
        double *work=NULL; int lwork=-1; double getlwork;
        /* Get workspace size */
        pdormqr_( "L", "T", &n, &s, &n, Aqr, &i1, &i1,
                  descA, tau, X, &i1, &i1, descB,
                  &getlwork, &lwork, &info );
        lwork = (int)getlwork;
        work = malloc( sizeof(double)*lwork );

        /* Actual call to ormqr */
        pdormqr_( "L", "T", &n, &s, &n, Aqr, &i1, &i1,
                  descA, tau, X, &i1, &i1, descB,
                  work, &lwork, &info );
        free(work);
    }
    pdtrsm_( "L", "U", "N", "N", &n, &s, &p1, Aqr, &i1, &i1, descA, X, &i1, &i1, descB );

    /* Compute B-AX */
    pdgemm_( "N", "N", &n, &s, &n, &m1, A, &i1, &i1, descA, X, &i1, &i1, descB,
             &p1, B, &i1, &i1, descB );
    {
        double *work = malloc( sizeof(double)*mloc );
        Xnorm = pdlange_( "I", &n, &s, X, &i1, &i1, descB, work );
        Rnorm = pdlange_( "I", &n, &s, B, &i1, &i1, descB, work );
        free( work );
    }

    if(iam == 0)
        printf("||A||oo = %e, ||X||oo = %e, ||B||oo = %e, ||R||oo = %e\n",
               Anorm, Xnorm, Bnorm, Rnorm );

    eps = pdlamch_( &ictxt, "Epsilon" );
    resid = Rnorm / ( (Bnorm + Anorm * Xnorm) * fmax(m, n) * eps );
    free( A ); free( B ); free( X );

    return resid;
}
