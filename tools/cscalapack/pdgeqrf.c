/*
 * Copyright (c) 2009-2020 The University of Tennessee and The University
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

// #define ASSERT_INFO
// #define USE_SUBMATRIX_A
// #define VERBOSE
// #define EXTRA_TOP
int main( int argc, char **argv ) {
    int params[PARAMS_SIZE];
    int info;
    int ictxt, nprow, npcol, myrow, mycol, iam;
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


#ifdef USE_SUBMATRIX_A
#ifdef EXTRA_TOP
    int extra_top = params[PARAM_EXTRA_ROWS];
#else
    int extra_top = 0;
#endif
    int extra_bottom = params[PARAM_EXTRA_ROWS];
    m = m + extra_top + extra_bottom;
#endif

#ifdef DPLASMA_WRAPPER_ON
    parsec_init_wrapper_();
#endif

    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
    mloc = numroc_( &m, &mb, &myrow, &i0, &nprow );
    nloc = numroc_( &n, &nb, &mycol, &i0, &npcol );
    descinit_( descA, &m, &n, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
#ifdef ASSERT_INFO
//    if((nb>=n)||(mb>=m)){
//        /*don't assert because some processes have not work*/
//    }else{
//        assert( 0 == info );
//    }
#endif

#ifdef USE_SUBMATRIX_A
    double * A_base = malloc( sizeof(double)*mloc*nloc );
    A = &(A_base[extra_top]);
    {
        int ii;
        int jj;
        for(jj=0; jj<extra_top; jj++){
            for (ii=0; ii<nloc; ii++){
                A_base[jj +  ii*(mloc)] = -111.0;//-t*rank-(a+1);
            }
        }

        for(jj=0; jj<extra_bottom; jj++){
            for (ii=0; ii<nloc; ii++){
                A_base[(mloc-jj-1) +  ii*(mloc)] = -222.0;//-t*rank-(a+1);
            }
        }
    }
#else
    A = malloc( sizeof(double)*mloc*nloc );
#endif

#ifdef USE_SUBMATRIX_A
    m = m - extra_top - extra_bottom;
#endif
    scalapack_pdplrnt( A,
                       m, n,
                       mb, nb,
                       myrow, mycol,
                       nprow, npcol,
                       mloc,
                       iseed );
#ifdef USE_SUBMATRIX_A
    m = m + extra_top + extra_bottom;
#endif

    tau = malloc( sizeof(double)*imin(m,n) );

#ifdef VERBOSE
    char name[10]="A";
    double * pwork = malloc( sizeof(double)*mloc*nloc);
    int nout=6;
#ifdef USE_SUBMATRIX_A
    pdlaprnt_(&m, &n, A_base, &i1, &i1, descA, &i0, &i0, name, &nout, pwork);
#else
    pdlaprnt_(&m, &n, A, &i1, &i1, descA, &i0, &i0, name, &nout, pwork);
#endif
#endif

    {
        double *work=NULL; int lwork=-1; double getlwork;
        double t1, t2;

#ifdef USE_SUBMATRIX_A
        descA[2] = descA[2] - extra_top - extra_bottom;
        m = m - extra_top - extra_bottom;
        pdgeqrf_( &m, &n, A, &i1, &i1, descA, tau, &getlwork, &lwork, &info );
        m = m + extra_top + extra_bottom;
        descA[2] = descA[2] + extra_top + extra_bottom;
#else
        pdgeqrf_( &m, &n, A, &i1, &i1, descA, tau, &getlwork, &lwork, &info );
#endif

#ifdef ASSERT_INFO
        assert( 0 == info );
#endif
        lwork = (int)getlwork;
        work = malloc( sizeof(double)*lwork );
        t1 = MPI_Wtime();

#ifdef USE_SUBMATRIX_A
        descA[2] = descA[2] - extra_top - extra_bottom;
        m = m - extra_top - extra_bottom;
        pdgeqrf_( &m, &n, A, &i1, &i1, descA, tau, work, &lwork, &info );
        m = m + extra_top + extra_bottom;
        descA[2] = descA[2] + extra_top + extra_bottom;
#else
        pdgeqrf_( &m, &n, A, &i1, &i1, descA, tau, work, &lwork, &info );
#endif
        assert( 0 == info );
        t2 = MPI_Wtime();
        telapsed = t2-t1;
        free(work); work = NULL;
    }
    if ( verif ) {
        resid = check_solution( params, A, tau );
    } else {
        resid = -1;
    }

    if( 0 != iam ) {
        MPI_Reduce( &telapsed, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    }
    else {
        MPI_Reduce( MPI_IN_PLACE, &telapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
        gflops = FLOPS_DGEQRF((double)m, (double)n)/1e+9/telapsed;
        pgflops = gflops/(((double)nprow)*((double)npcol));
        printf( "### PDGEQRF ###\n"
                "#%4sx%-4s %7s %7s %4s %4s %4s # %10s %10s %10s %11s\n", "P", "Q", "M", "N", "MB", "NB", "NRHS", "resid", "time(s)", "gflops", "gflops/pxq" );
        printf( " %4d %-4d %7d %7d %4d %4d %4d   %10.3e %10.3g %10.3g %11.3g\n", nprow, npcol, m, n, mb, nb, s, resid, telapsed, gflops, pgflops );
    }

#ifdef VERBOSE
#ifdef USE_SUBMATRIX_A
    pdlaprnt_(&m, &n, A_base, &i1, &i1, descA, &i0, &i0, name, &nout, pwork);
#else
    pdlaprnt_(&m, &n, A, &i1, &i1, descA, &i0, &i0, name, &nout, pwork);
#endif
#endif

#ifdef DPLASMA_WRAPPER_ON
    parsec_fini_wrapper_();
#endif

#ifdef EXTRA_TOP
    free( A_base ); A_base = NULL;
#else
    free( A ); A = NULL;
#endif
    free( tau ); tau = NULL;

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
#ifdef ASSERT_INFO
//    if((nb>=n)||(mb>=m)){
//        /*don't assert because some processes have not work*/
//    }else{
//        assert( 0 == info );
//    }
#endif

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
#ifdef ASSERT_INFO
//    if((nb>=s)||(mb>=n)){
//        /*don't assert because some processes have not work*/
//    }else{
//        assert( 0 == info );
//    }
#endif

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
