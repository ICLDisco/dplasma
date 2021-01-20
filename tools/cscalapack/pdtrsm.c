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
//workaround dplasmaLeft
#include "../../src/include/dplasma/constants.h"


int main( int argc, char **argv ) {
    int params[PARAMS_SIZE];
    int info;
    int ictxt, nprow, npcol, myrow, mycol, iam;
    int m, n, mb, nb, s, mloc, nloc, verif, iseed;
    double *A=NULL, *B=NULL; int descA[9], descB[9];
    double resid = NAN;
    double telapsed, gflops, pgflops;

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
    int number_runs = params[PARAM_NRUNS];

    int Aseed = 3872;
    int Bseed = 2873;
    double alpha =  3.5;

#ifdef DPLASMA_WRAPPER_ON
    parsec_init_wrapper_();
#endif

    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

    if( verif ) {
        /* TODO: check devel on wrapper  */
    } else { // NOT VERIF
      int Am, An, Ai, Aj, Amb, Anb;
      int Bm, Bn, Bi, Bj, Bmb, Bnb;

      char side  = 'R';
      char uplo  = 'U';
      char trans = 'N';
      char diag  = 'N';

      // char side  = 'L';
      // char uplo  = 'L';
      // char trans = 'N';
      // char diag  = 'U';

      // char side  = 'L';
      // char uplo  = 'U';
      // char trans = 'N';
      // char diag  = 'N';

      printf("TRSM side %c uplo %c trans %c diag %c \n", side, uplo, trans, diag);
      if ( side == 'L' ) {
          Am = m; An = m;
      } else {
          Am = n; An = n;
      }
      Bm = m; Bn = n;
      Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
      mloc = numroc_( &Am, &mb, &myrow, &i0, &nprow );
      nloc = numroc_( &An, &nb, &mycol, &i0, &npcol );
      descinit_( descA, &Am, &An, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
      assert( 0 == info );
      A = malloc( sizeof(double)*mloc*nloc );
      scalapack_pdplrnt( A,
                         Am, An,
                         mb, nb,
                         myrow, mycol,
                         nprow, npcol,
                         mloc,
                         Aseed);
    // dplasma_dplgsy( parsec, 0., dplasmaUpperLower, (parsec_tiled_matrix_dc_t *)&dcA, Aseed);
    // dplasma_dplrnt( parsec, 0, (parsec_tiled_matrix_dc_t *)&dcC, Cseed );

      A = malloc( sizeof(double)*mloc*nloc );
      pdlacpy_( "All", &Am, &An, A, &i1, &i1, descA, A, &i1, &i1, descA );

      mloc = numroc_( &Bm, &mb, &myrow, &i0, &nprow );
      nloc = numroc_( &Bn, &nb, &mycol, &i0, &npcol );
      descinit_( descB, &Bm, &Bn, &mb, &nb, &i0, &i0, &ictxt, &mloc, &info );
      assert( 0 == info );
      B = malloc( sizeof(double)*mloc*nloc );
      scalapack_pdplrnt( B,
                         Bm, Bn,
                         mb, nb,
                         myrow, mycol,
                         nprow, npcol,
                         mloc,
                         Bseed);

      int t;
      for( t = 0; t < number_runs; t++ ) {

  #ifdef DPLASMA_WRAPPER_ON
          parsec_wrapper_devices_release_memory_();
  #endif
          double t1, t2;
          t1 = MPI_Wtime();

          pdtrsm_(&side, &uplo, &trans, &diag,
                   &m, &n, &alpha,
                   A, &i1, &i1, descA,
                  B, &i1, &i1, descB);

          t2 = MPI_Wtime();
          telapsed = t2-t1;
          if( 0 != iam ) {
              MPI_Reduce( &telapsed, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
          }
          else {
              MPI_Reduce( MPI_IN_PLACE, &telapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
              gflops = FLOPS_DTRSM((side=='L'? dplasmaLeft : dplasmaRight), (double)m, (double)n)/1e+9/telapsed;
              pgflops = gflops/(((double)nprow)*((double)npcol));
          }

          if( 0 == iam ) {
              printf("[****] TIMEHL(s) %12.5f : dtrsm \tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops"
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
    }

#ifdef DPLASMA_WRAPPER_ON
    parsec_fini_wrapper_();
#endif

    Cblacs_exit( 0 );
    return 0;
}
