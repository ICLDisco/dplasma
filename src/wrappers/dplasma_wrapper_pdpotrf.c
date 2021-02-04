#include "common.h"

/*      SUBROUTINE PDPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
 *
 *     .. Scalar Arguments ..
 *       CHARACTER          UPLO
 *       INTEGER            IA, INFO, JA, N
 *     ..
 *     .. Array Arguments ..
 *      INTEGER            DESCA( * )
 *      DOUBLE PRECISION   A( * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  PDPOTRF computes the Cholesky factorization of an N-by-N real
 *  symmetric positive definite distributed matrix sub( A ) denoting
 *  A(IA:IA+N-1, JA:JA+N-1).
 *
 *  The factorization has the form
 *
 *            sub( A ) = U' * U ,  if UPLO = 'U', or
 *
 *            sub( A ) = L  * L',  if UPLO = 'L',
 *
 *  where U is an upper triangular matrix and L is lower triangular.
 *
 *  Notes
 *  =====
 *
 *  Each global data object is described by an associated description
 *  vector.  This vector stores the information required to establish
 *  the mapping between an object element and its corresponding process
 *  and memory location.
 *
 *  Let A be a generic term for any 2D block cyclicly distributed array.
 *  Such a global array has an associated description vector DESCA.
 *  In the following comments, the character _ should be read as
 *  "of the global array".
 *
 *  NOTATION        STORED IN      EXPLANATION
 *  --------------- -------------- --------------------------------------
 *  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
 *                                 DTYPE_A = 1.
 *  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
 *                                 the BLACS process grid A is distribu-
 *                                 ted over. The context itself is glo-
 *                                 bal, but the handle (the integer
 *                                 value) may vary.
 *  M_A    (global) DESCA( M_ )    The number of rows in the global
 *                                 array A.
 *  N_A    (global) DESCA( N_ )    The number of columns in the global
 *                                 array A.
 *  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
 *                                 the rows of the array.
 *  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
 *                                 the columns of the array.
 *  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
 *                                 row of the array A is distributed.
 *  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
 *                                 first column of the array A is
 *                                 distributed.
 *  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
 *                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
 *
 *  Let K be the number of rows or columns of a distributed matrix,
 *  and assume that its process grid has dimension p x q.
 *  LOCr( K ) denotes the number of elements of K that a process
 *  would receive if K were distributed over the p processes of its
 *  process column.
 *  Similarly, LOCc( K ) denotes the number of elements of K that a
 *  process would receive if K were distributed over the q processes of
 *  its process row.
 *  The values of LOCr() and LOCc() may be determined via a call to the
 *  ScaLAPACK tool function, NUMROC:
 *          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
 *          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
 *  An upper bound for these quantities may be computed by:
 *          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
 *          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
 *
 *  This routine requires square block decomposition ( MB_A = NB_A ).
 *
 *  Arguments
 *  =========
 *
 *  UPLO    (global input) CHARACTER
 *          = 'U':  Upper triangle of sub( A ) is stored;
 *          = 'L':  Lower triangle of sub( A ) is stored.
 *
 *  N       (global input) INTEGER
 *          The number of rows and columns to be operated on, i.e. the
 *          order of the distributed submatrix sub( A ). N >= 0.
 *
 *  A       (local input/local output) DOUBLE PRECISION pointer into the
 *          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
 *          On entry, this array contains the local pieces of the
 *          N-by-N symmetric distributed matrix sub( A ) to be factored.
 *          If UPLO = 'U', the leading N-by-N upper triangular part of
 *          sub( A ) contains the upper triangular part of the matrix,
 *          and its strictly lower triangular part is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          sub( A ) contains the lower triangular part of the distribu-
 *          ted matrix, and its strictly upper triangular part is not
 *          referenced. On exit, if UPLO = 'U', the upper triangular
 *          part of the distributed matrix contains the Cholesky factor
 *          U, if UPLO = 'L', the lower triangular part of the distribu-
 *          ted matrix contains the Cholesky factor L.
 *
 *  IA      (global input) INTEGER
 *          The row index in the global array A indicating the first
 *          row of sub( A ).
 *
 *  JA      (global input) INTEGER
 *          The column index in the global array A indicating the
 *          first column of sub( A ).
 *
 *  DESCA   (global and local input) INTEGER array of dimension DLEN_.
 *          The array descriptor for the distributed matrix A.
 *
 *  INFO    (global output) INTEGER
 *          = 0:  successful exit
 *          < 0:  If the i-th argument is an array and the j-entry had
 *                an illegal value, then INFO = -(i*100+j), if the i-th
 *                argument is a scalar and had an illegal value, then
 *                INFO = -i.
 *          > 0:  If INFO = K, the leading minor of order K,
 *                A(IA:IA+K-1,JA:JA+K-1) is not positive definite, and
 *                the factorization could not be completed.
 */

void pdpotrf_w(char * UPLO,
              int * N,
              double * A,
              int * IA,
              int * JA,
              int * DESCA,
              int * info){


#ifdef COUNT_WRAPPED_CALLS
    count_PDPOTRF++;
#endif

    *info=0;
    if(*N == 0){
      /* NOP */
      return;
    }

    int KP = 1;
    int KQ = 1;

    PASTE_SETUP(A);

#ifdef WRAPPER_VERBOSE_CALLS
    if(rank_A == 0){
      printf("V-PDPOTRF N%d "
             "IA%d JA%d A%p MBA%d NBA%d %c \n",
             *N,
             *IA, *JA, A, DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_], *UPLO);
    }
#endif

    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  " N%d IA%d JA%d (ictxt)DESCA[WRAPPER_CTXT1_] %d, "
          "(gM)DESCA[WRAPPER_M1_] %d, (gN)DESCA[WRAPPER_N1_] %d, (MB)DESCA[WRAPPER_MB1_] %d, (NB)DESCA[WRAPPER_NB1_] %d, "
          "DESCA[WRAPPER_RSRC1_] %d, DESCA[WRAPPER_CSRC1_] %d, (LDD)DESCA[WRAPPER_LLD1_] %d\n",
          *N, *IA, *JA, DESCA[WRAPPER_CTXT1_],
          DESCA[WRAPPER_M1_], DESCA[WRAPPER_N1_], DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_],
          DESCA[WRAPPER_RSRC1_], DESCA[WRAPPER_CSRC1_], DESCA[WRAPPER_LLD1_]);

    parsec_init_wrapped_call((void*)comm_A);

    dplasma_enum_t uplo_parsec = OP_UPLO(*UPLO);

    two_dim_block_cyclic_t dcA_lapack;
    two_dim_block_cyclic_lapack_init(&dcA_lapack, matrix_RealDouble, matrix_Lapack,
                                      rank_A,
                                      MB_A, NB_A,
                                      gM_A, gN_A,
                                      cIA, cJA,
                                      *N, *N,
                                      P_A, Q_A,
                                      KP, KQ,
                                      iP_A, jQ_A,
                                      LLD_A, nloc_A);
    dcA_lapack.mat = A;
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcA_lapack, "dcA_lapack");

#ifdef CHECK_RESULTS
    int check=1;
    int loud=5;
    PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
        two_dim_block_cyclic, (&dcA0, matrix_RealDouble, matrix_Tile,
                               rank_A, MB_A, NB_A, *N, *N, 0, 0,
                               *N, *N,
                               P_A, Q_A,
                               KP, KQ,
                               0, 0));

    if( check ) {
        dcopy_lapack_tile(parsec_ctx, &dcA_lapack, &dcA0, mloc_A, nloc_A);
    }
#endif

    PRINT(parsec_ctx, comm_A, uplo_parsec, "dcA", (&dcA_lapack));

    int redisA = 0;
    if( (cIA % MB_A != 0) || ( cJA % NB_A != 0)) redisA = 1;
    assert(redisA == 0); /* not aligned offsets are not supported for POTRF */

    two_dim_block_cyclic_t *dcA = redistribute_lapack_input(&dcA_lapack, redisA, comm_A, rank_A, "redisA");

#ifdef MEASURE_INTERNAL_TIMES
    PASTE_CODE_FLOPS(FLOPS_DPOTRF, ((DagDouble_t)*N));
#endif

    WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec_ctx, dpotrf,
                              (uplo_parsec, (parsec_tiled_matrix_dc_t*)dcA, info),
                              dplasma_dpotrf_Destruct( PARSEC_dpotrf ),
                              rank_A, P_A, Q_A, NB_A, gN_A, comm_A);

    dcA = redistribute_lapack_output_cleanup(&dcA_lapack, dcA, 1, comm_A, rank_A, "redisA");

    PRINT(parsec_ctx, comm_A, uplo_parsec, "dcA", dcA);

    if( 0 == rank_A && *info != 0 ) {
      printf("-- Factorization is suspicious (info = %d) ! \n", *info);
    }

#ifdef CHECK_RESULTS
    if( check ) {
        /* Check the factorization */
        PASTE_CODE_ALLOCATE_MATRIX(dcA_out, check,
            two_dim_block_cyclic, (&dcA_out, matrix_RealDouble, matrix_Tile,
                                   rank_A, MB_A, NB_A, *N, *N, 0, 0,
                                   *N, *N,
                                   P_A, Q_A,
                                   KP, KQ,
                                   0, 0));
        dcopy_lapack_tile(parsec_ctx, dcA, &dcA_out, mloc_A, nloc_A);

        int ret = 0;
        ret |= check_dpotrf( parsec_ctx, (rank_A == 0) ? loud : 0, uplo_parsec,
                             (parsec_tiled_matrix_dc_t *)&dcA_out,
                             (parsec_tiled_matrix_dc_t *)&dcA0);

        int NRHS = 1;
        int LDB   = *N;
        int random_seed = 3872;

        /* Check the solution */
        PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
            two_dim_block_cyclic, (&dcB, matrix_RealDouble, matrix_Tile,
                                   rank_A, MB_A, NB_A, LDB, NRHS, 0, 0,
                                   *N, NRHS, P_A, Q_A, KP, KQ, 0, 0));
        dplasma_dplrnt( parsec_ctx, 0, (parsec_tiled_matrix_dc_t *)&dcB, random_seed+1);

        PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
            two_dim_block_cyclic, (&dcX, matrix_RealDouble, matrix_Tile,
                                   rank_A, MB_A, NB_A, LDB, NRHS, 0, 0,
                                   *N, NRHS, P_A, Q_A, KP, KQ, 0, 0));
        dplasma_dlacpy( parsec_ctx, PlasmaUpperLower,
                        (parsec_tiled_matrix_dc_t *)&dcB, (parsec_tiled_matrix_dc_t *)&dcX );

        dplasma_dpotrs(parsec_ctx, uplo_parsec,
                       (parsec_tiled_matrix_dc_t *)&dcA_out,
                       (parsec_tiled_matrix_dc_t *)&dcX );

        ret |= check_daxmb( parsec_ctx, (rank_A == 0) ? loud : 0, uplo_parsec,
                            (parsec_tiled_matrix_dc_t *)&dcA0,
                            (parsec_tiled_matrix_dc_t *)&dcB,
                            (parsec_tiled_matrix_dc_t *)&dcX);

        /* Cleanup */
        parsec_data_free(dcA_out.mat); dcA_out.mat = NULL;
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA_out );
        parsec_data_free(dcA0.mat); dcA0.mat = NULL;
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA0 );
        parsec_data_free(dcB.mat); dcB.mat = NULL;
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcB );
        parsec_data_free(dcX.mat); dcX.mat = NULL;
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcX );
    }
#endif

    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)dcA);
}

GENERATE_F77_BINDINGS (PDPOTRF,
                       pdpotrf,
                       pdpotrf_,
                       pdpotrf__,
                       pdpotrf_w,
                       (char * UPLO, int * N, double * A, int * IA, int * JA, int * DESCA, int * info),
                       (UPLO, N, A, IA, JA, DESCA, info))
