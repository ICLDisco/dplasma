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


int main( int argc, char **argv ) {
    int params[PARAMS_SIZE];
    int info;
    int ictxt, nprow, npcol, myrow, mycol, iam;
    int number_runs;
    int m, n, k, mb, nb, s, mloc, nloc, verif, iseed;
    double *A=NULL, *B=NULL, *C=NULL; int descA[9], descB[9], descC[9];
    double resid = NAN;
    double telapsed, gflops, pgflops;

    setup_params( params, argc, argv );
    ictxt = params[PARAM_BLACS_CTX];
    iam   = params[PARAM_RANK];
    m     = params[PARAM_M];
    n     = params[PARAM_N];
    k     = params[PARAM_K];
    mb    = params[PARAM_MB];
    nb    = params[PARAM_NB];
    s     = params[PARAM_NRHS];
    iseed = params[PARAM_SEED];
    verif = params[PARAM_VALIDATE];
    number_runs = params[PARAM_NRUNS];

    int Aseed = 3872;
    int Bseed = 4674;
    int Cseed = 2873;
    double alpha =  0.51;
    double beta  = -0.42;

#ifdef DPLASMA_WRAPPER_ON
    parsec_init_wrapper_();
#endif

    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

    if(verif){
        int tA; int tB;
        char aTA[3] = "nt";
        char aTB[3] = "nt";

        for(tA=0; tA<2; tA++) {
          for(tB=0; tB<2; tB++) {
            char TA = aTA[tA];
            char TB = aTB[tB];
            int Am, An, Ai, Aj, Amb, Anb;
            int Bm, Bn, Bi, Bj, Bmb, Bnb;
            int Cm, Cn;
            if ( TA == 'n') {
                Am  = m;
                An  = k;
            } else {
                Am  = k;
                An  = m;
            }
            if ( TB == 'n') {
                Bm  = k;
                Bn  = n;
            } else {
                Bm  = k;
                Bn  = n;
            }

            Cm = m;
            Cn = n;

            mloc = numroc_( &Am, &mb, &myrow, &i0, &nprow );
            nloc = numroc_( &An, &nb, &mycol, &i0, &npcol );
            descinit_( descA, &Am, &An, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
            assert( 0 == info );
            A = malloc( sizeof(double)*mloc*nloc );
            scalapack_pdplrnt(A,
                               Am, An,
                               mb, nb,
                               myrow, mycol,
                               nprow, npcol,
                               mloc,
                               Aseed);

            mloc = numroc_( &Bm, &mb, &myrow, &i0, &nprow );
            nloc = numroc_( &Bn, &nb, &mycol, &i0, &npcol );
            descinit_( descB, &Bm, &Bn, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
            assert( 0 == info );
            B = malloc( sizeof(double)*mloc*nloc );
            scalapack_pdplrnt(B,
                               Bm, Bn,
                               mb, nb,
                               myrow, mycol,
                               nprow, npcol,
                               mloc,
                               Bseed);

            mloc = numroc_( &Cm, &mb, &myrow, &i0, &nprow );
            nloc = numroc_( &Cn, &nb, &mycol, &i0, &npcol );
            descinit_( descC, &Cm, &Cn, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
            C = malloc( sizeof(double)*mloc*nloc );
            scalapack_pdplrnt(C,
                               Cm, Cn,
                               mb, nb,
                               myrow, mycol,
                               nprow, npcol,
                               mloc,
                               Cseed);
            pdgemm_( &TA, &TB, &m, &n, &k, &alpha, A, &i1, &i1, descA,
                                            B, &i1, &i1, descB,
                                     &beta, C, &i1, &i1, descC );

            MPI_Barrier(MPI_COMM_WORLD);
            if( 0 == iam ) {
              printf("[****] completed dgemm [%c %c] \tPxQ= %3d %-3d NB= %4d N= %7d\n",
                        TA, TB, nprow, npcol, nb, n);
            }
            free( A ); A = NULL;
            free( B ); B = NULL;
            free( C ); C = NULL;
          }
        }

    }else{// NOT VERIF
      char TA = 'n';
      char TB = 'n';
      int Am, An, Ai, Aj, Amb, Anb;
      int Bm, Bn, Bi, Bj, Bmb, Bnb;
      int Cm, Cn;
      if ( TA == 'n') {
          Am  = m;
          An  = k;
      } else {
          Am  = k;
          An  = m;
      }
      if ( TB == 'n') {
          Bm  = k;
          Bn  = n;
      } else {
          Bm  = k;
          Bn  = n;
      }

      Cm = m;
      Cn = n;

      Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
      mloc = numroc_( &Am, &mb, &myrow, &i0, &nprow );
      nloc = numroc_( &An, &nb, &mycol, &i0, &npcol );
      descinit_( descA, &Am, &An, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
      assert( 0 == info );
      A = malloc( sizeof(double)*mloc*nloc );
      scalapack_pdplrnt(A,
                         Am, An,
                         mb, nb,
                         myrow, mycol,
                         nprow, npcol,
                         mloc,
                         Aseed);

      mloc = numroc_( &Bm, &mb, &myrow, &i0, &nprow );
      nloc = numroc_( &Bn, &nb, &mycol, &i0, &npcol );
      descinit_( descB, &Bm, &Bn, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
      assert( 0 == info );
      B = malloc( sizeof(double)*mloc*nloc );
      scalapack_pdplrnt(B,
                         Bm, Bn,
                         mb, nb,
                         myrow, mycol,
                         nprow, npcol,
                         mloc,
                         Bseed);

      mloc = numroc_( &Cm, &mb, &myrow, &i0, &nprow );
      nloc = numroc_( &Cn, &nb, &mycol, &i0, &npcol );
      descinit_( descC, &Cm, &Cn, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
      C = malloc( sizeof(double)*mloc*nloc );
      scalapack_pdplrnt(C,
                         Cm, Cn,
                         mb, nb,
                         myrow, mycol,
                         nprow, npcol,
                         mloc,
                         Cseed);

      int t;
      for(t = 0; t < number_runs; t++) {

  #ifdef DPLASMA_WRAPPER_ON
          parsec_wrapper_devices_release_memory_();
  #endif
          double t1, t2;
          t1 = MPI_Wtime();

          pdgemm_( &TA, &TB, &m, &n, &k, &alpha, A, &i1, &i1, descA,
                                          B, &i1, &i1, descB,
                                   &beta, C, &i1, &i1, descC );

          t2 = MPI_Wtime();
          telapsed = t2-t1;
          if( 0 != iam ) {
              MPI_Reduce( &telapsed, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
          }
          else {
              MPI_Reduce( MPI_IN_PLACE, &telapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
              gflops = FLOPS_DGEMM((double)m, (double)n, (double)n)/1e+9/telapsed;
              pgflops = gflops/(((double)nprow)*((double)npcol));
          }

          if( 0 == iam ) {
              printf("[****] TIMEHL(s) %12.5f : dgemm \tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops"
                    " - ENQ&PROG&DEST %12.5f : %14f gflops"
                    " - ENQ %12.5f - DEST %12.5f\n",
                            telapsed, nprow, npcol, nb, n,
                            gflops,
                            telapsed,
                            gflops,
                            0.0,0.0);
          }

  #ifdef DPLASMA_WRAPPER_ON
          parsec_wrapper_devices_reset_load_();
  #endif

      }
      free( A ); A = NULL;
      free( B ); B = NULL;
      free( C ); C = NULL;
    }

#ifdef DPLASMA_WRAPPER_ON
    parsec_fini_wrapper_();
#endif

    Cblacs_exit( 0 );
    return 0;
}
