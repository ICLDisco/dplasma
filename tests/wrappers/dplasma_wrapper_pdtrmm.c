#include "common.h"

 /*
 *  Purpose
 *  =======
 *
 *  PDTRMM  performs one of the matrix-matrix operations
 *
 *     sub( B ) := alpha * op( sub( A ) ) * sub( B ),
 *
 *     or
 *
 *     sub( B ) := alpha * sub( B ) * op( sub( A ) ),
 *
 *  where
 *
 *     sub( A ) denotes A(IA:IA+M-1,JA:JA+M-1)  if SIDE = 'L',
 *                      A(IA:IA+N-1,JA:JA+N-1)  if SIDE = 'R', and,
 *
 *     sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1).
 *
 *  Alpha  is a scalar,  sub( B )  is an m by n submatrix,  sub( A ) is a
 *  unit, or non-unit, upper or lower triangular submatrix and op( X ) is
 *  one of
 *
 *     op( X ) = X   or   op( X ) = X'.
 *
 *  Notes
 *  =====
 *
 *  A description  vector  is associated with each 2D block-cyclicly dis-
 *  tributed matrix.  This  vector  stores  the  information  required to
 *  establish the  mapping  between a  matrix entry and its corresponding
 *  process and memory location.
 *
 *  In  the  following  comments,   the character _  should  be  read  as
 *  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
 *  block cyclicly distributed matrix.  Its description vector is DESC_A:
 *
 *  NOTATION         STORED IN       EXPLANATION
 *  ---------------- --------------- ------------------------------------
 *  DTYPE_A (global) DESCA[ DTYPE_ ] The descriptor type.
 *  CTXT_A  (global) DESCA[ CTXT_  ] The BLACS context handle, indicating
 *                                   the NPROW x NPCOL BLACS process grid
 *                                   A  is  distributed over. The context
 *                                   itself  is  global,  but  the handle
 *                                   (the integer value) may vary.
 *  M_A     (global) DESCA[ M_     ] The  number of rows in the distribu-
 *                                   ted matrix A, M_A >= 0.
 *  N_A     (global) DESCA[ N_     ] The number of columns in the distri-
 *                                   buted matrix A, N_A >= 0.
 *  IMB_A   (global) DESCA[ IMB_   ] The number of rows of the upper left
 *                                   block of the matrix A, IMB_A > 0.
 *  INB_A   (global) DESCA[ INB_   ] The  number  of columns of the upper
 *                                   left   block   of   the  matrix   A,
 *                                   INB_A > 0.
 *  MB_A    (global) DESCA[ MB_    ] The blocking factor used to  distri-
 *                                   bute the last  M_A-IMB_A  rows of A,
 *                                   MB_A > 0.
 *  NB_A    (global) DESCA[ NB_    ] The blocking factor used to  distri-
 *                                   bute the last  N_A-INB_A  columns of
 *                                   A, NB_A > 0.
 *  RSRC_A  (global) DESCA[ RSRC_  ] The process row over which the first
 *                                   row of the matrix  A is distributed,
 *                                   NPROW > RSRC_A >= 0.
 *  CSRC_A  (global) DESCA[ CSRC_  ] The  process column  over  which the
 *                                   first column of  A  is  distributed.
 *                                   NPCOL > CSRC_A >= 0.
 *  LLD_A   (local)  DESCA[ LLD_   ] The  leading dimension  of the local
 *                                   array  storing  the  local blocks of
 *                                   the distributed matrix A,
 *                                   IF( Lc( 1, N_A ) > 0 )
 *                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
 *                                   ELSE
 *                                      LLD_A >= 1.
 *
 *  Let K be the number of  rows of a matrix A starting at the global in-
 *  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
 *  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
 *  receive if these K rows were distributed over NPROW processes.  If  K
 *  is the number of columns of a matrix  A  starting at the global index
 *  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
 *  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
 *  these K columns were distributed over NPCOL processes.
 *
 *  The values of Lr() and Lc() may be determined via a call to the func-
 *  tion PB_Cnumroc:
 *  Lr( IA, K ) = PB_Cnumroc( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
 *  Lc( JA, K ) = PB_Cnumroc( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
 *
 *  Arguments
 *  =========
 *
 *  SIDE    (global input) CHARACTER*1
 *          On entry,  SIDE  specifies whether  op( sub( A ) ) multiplies
 *          sub( B ) from the left or right as follows:
 *
 *          SIDE = 'L' or 'l'  sub( B ) := alpha*op( sub( A ) )*sub( B ),
 *
 *          SIDE = 'R' or 'r'  sub( B ) := alpha*sub( B )*op( sub( A ) ).
 *
 *  UPLO    (global input) CHARACTER*1
 *          On entry,  UPLO  specifies whether the submatrix  sub( A ) is
 *          an upper or lower triangular submatrix as follows:
 *
 *             UPLO = 'U' or 'u'   sub( A ) is an upper triangular
 *                                 submatrix,
 *
 *             UPLO = 'L' or 'l'   sub( A ) is a  lower triangular
 *                                 submatrix.
 *
 *  TRANSA  (global input) CHARACTER*1
 *          On entry,  TRANSA  specifies the form of op( sub( A ) ) to be
 *          used in the matrix multiplication as follows:
 *
 *             TRANSA = 'N' or 'n'   op( sub( A ) ) = sub( A ),
 *
 *             TRANSA = 'T' or 't'   op( sub( A ) ) = sub( A )',
 *
 *             TRANSA = 'C' or 'c'   op( sub( A ) ) = sub( A )'.
 *
 *  DIAG    (global input) CHARACTER*1
 *          On entry,  DIAG  specifies  whether or not  sub( A )  is unit
 *          triangular as follows:
 *
 *             DIAG = 'U' or 'u'  sub( A )  is  assumed to be unit trian-
 *                                gular,
 *
 *             DIAG = 'N' or 'n'  sub( A ) is not assumed to be unit tri-
 *                                angular.
 *
 *  M       (global input) INTEGER
 *          On entry,  M  specifies the number of rows of  the  submatrix
 *          sub( B ). M  must be at least zero.
 *
 *  N       (global input) INTEGER
 *          On entry, N  specifies the number of columns of the submatrix
 *          sub( B ). N  must be at least zero.
 *
 *  ALPHA   (global input) DOUBLE PRECISION
 *          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *          supplied  as  zero  then  the  local entries of  the array  B
 *          corresponding to the entries of the submatrix  sub( B )  need
 *          not be set on input.
 *
 *  A       (local input) DOUBLE PRECISION array
 *          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
 *          at least Lc( 1, JA+M-1 ) when  SIDE = 'L' or 'l'   and  is at
 *          least Lc( 1, JA+N-1 ) otherwise.  Before  entry,  this  array
 *          contains the local entries of the matrix A.
 *          Before entry with  UPLO = 'U' or 'u', this array contains the
 *          local entries corresponding to  the entries of the upper tri-
 *          angular submatrix  sub( A ), and the local entries correspon-
 *          ding to the entries of the strictly lower triangular part  of
 *          the submatrix  sub( A )  are not referenced.
 *          Before entry with  UPLO = 'L' or 'l', this array contains the
 *          local entries corresponding to  the entries of the lower tri-
 *          angular submatrix  sub( A ), and the local entries correspon-
 *          ding to the entries of the strictly upper triangular part  of
 *          the submatrix  sub( A )  are not referenced.
 *          Note  that  when DIAG = 'U' or 'u', the local entries corres-
 *          ponding to the  diagonal elements  of the submatrix  sub( A )
 *          are not referenced either, but are assumed to be unity.
 *
 *  IA      (global input) INTEGER
 *          On entry, IA  specifies A's global row index, which points to
 *          the beginning of the submatrix sub( A ).
 *
 *  JA      (global input) INTEGER
 *          On entry, JA  specifies A's global column index, which points
 *          to the beginning of the submatrix sub( A ).
 *
 *  DESCA   (global and local input) INTEGER array
 *          On entry, DESCA  is an integer array of dimension DLEN_. This
 *          is the array descriptor for the matrix A.
 *
 *  B       (local input/local output) DOUBLE PRECISION array
 *          On entry, B is an array of dimension (LLD_B, Kb), where Kb is
 *          at least Lc( 1, JB+N-1 ).  Before  entry, this array contains
 *          the local entries of the matrix  B.
 *          On exit, the local entries of this array corresponding to the
 *          to  the entries of the submatrix sub( B ) are  overwritten by
 *          the local entries of the m by n transformed submatrix.
 *
 *  IB      (global input) INTEGER
 *          On entry, IB  specifies B's global row index, which points to
 *          the beginning of the submatrix sub( B ).
 *
 *  JB      (global input) INTEGER
 *          On entry, JB  specifies B's global column index, which points
 *          to the beginning of the submatrix sub( B ).
 *
 *  DESCB   (global and local input) INTEGER array
 *          On entry, DESCB  is an integer array of dimension DLEN_. This
 *          is the array descriptor for the matrix B.
 *
 *  -- Written on April 1, 1998 by
 *     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
 *
 *  ---------------------------------------------------------------------
 */
#ifdef CHECK_RESULTS
static int check_solution( parsec_context_t *parsec, int loud,
                           dplasma_enum_t side, dplasma_enum_t uplo, dplasma_enum_t trans, dplasma_enum_t diag,
                           double alpha,
                           int Am, int An, int Aseed,
                           int M,  int N,  int Cseed,
                           two_dim_block_cyclic_t *dcCfinal );
#endif
void pdtrmm_w( char* SIDE, char* UPLO, char* TRANS, char* DIAG,
               int * M, int * N, double * ALPHA,
               double * A, int * IA, int * JA, int * DESCA,
               double * B, int * IB, int * JB, int * DESCB )
{

#ifdef COUNT_WRAPPED_CALLS
    count_PDTRMM++;
#endif
    if( (*M == 0) || (*N == 0)){
      /* NOP */
      return;
    }

    int i0=0;
    int KP = 1;
    int KQ = 1;

    dplasma_enum_t side  = OP_SIDE(*SIDE);
    dplasma_enum_t uplo  = OP_UPLO(*UPLO);
    dplasma_enum_t trans = OP_TRANS(*TRANS);
    dplasma_enum_t diag  = OP_DIAG(*DIAG);

    PASTE_SETUP(A);
    PASTE_SETUP(B);

#ifdef WRAPPER_VERBOSE_CALLS
    if(rank_A == 0){
      PARSEC_DEBUG_VERBOSE(0, parsec_debug_output,  "PDTRMM M%d N%d "
                              "IA%d JA%d A%p MBA%d NBA%d "
                              "IB%d JB%d B%p MBB%d NBB%d %c%c%c%c ",
                              *M, *N,
                              *IA, *JA, A, DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_],
                              *IB, *JB, B, DESCB[WRAPPER_MB1_], DESCB[WRAPPER_NB1_],
                              *SIDE, *UPLO, *TRANS, *DIAG);
    }
#endif

    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  " M%d N%d IA%d JA%d (ictxt)DESCA[WRAPPER_CTXT1_] %d, "
          "(gM)DESCA[WRAPPER_M1_] %d, (gN)DESCA[WRAPPER_N1_] %d, (MB)DESCA[WRAPPER_MB1_] %d, (NB)DESCA[WRAPPER_NB1_] %d, "
          "DESCA[WRAPPER_RSRC1_] %d, DESCA[WRAPPER_CSRC1_] %d, (LDD)DESCA[WRAPPER_LLD1_] %d "
          "SIDE %c, UPLO %c, TRANS %c, DIAG %c mloc_A %d nloc_A %d\n",
          *M, *N, *IA, *JA, DESCA[WRAPPER_CTXT1_],
          DESCA[WRAPPER_M1_], DESCA[WRAPPER_N1_], DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_],
          DESCA[WRAPPER_RSRC1_], DESCA[WRAPPER_CSRC1_], DESCA[WRAPPER_LLD1_],
          *SIDE, *UPLO, *TRANS, *DIAG, mloc_A, nloc_A);
    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  " M%d N%d IB%d JB%d (ictxt)DESCB[WRAPPER_CTXT1_] %d, "
          "(gM)DESCB[WRAPPER_M1_] %d, (gN)DESCB[WRAPPER_N1_] %d, (MB)DESCB[WRAPPER_MB1_] %d, (NB)DESCB[WRAPPER_NB1_] %d, "
          "DESCB[WRAPPER_RSRC1_] %d, DESCB[WRAPPER_CSRC1_] %d, (LDD)DESCB[WRAPPER_LLD1_] %d "
          "SIDE %c, UPLO %c, TRANS %c, DIAG %c mloc_B %d nloc_B %d\n",
          *M, *N, *IB, *JB, DESCB[WRAPPER_CTXT1_],
          DESCB[WRAPPER_M1_], DESCB[WRAPPER_N1_], DESCB[WRAPPER_MB1_], DESCB[WRAPPER_NB1_],
          DESCB[WRAPPER_RSRC1_], DESCB[WRAPPER_CSRC1_], DESCB[WRAPPER_LLD1_],
          *SIDE, *UPLO, *TRANS, *DIAG, mloc_B, nloc_B);

    assert(comm_index_A == comm_index_B);
    parsec_init_wrapped_call((void*)comm_A);

    int Am, An, Bm, Bn;
    int Am_loc, An_loc, Bm_loc, Bn_loc;
    int AGm, AGn, BGm, BGn;

    if ( side == dplasmaLeft ) {
        Am = *M; An = *M;
    } else {
        Am = *N; An = *N;
    }

    Bm = *M; Bn = *N;
    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,
        "A-%c %dx%d TRMM B %dx%d",
        *TRANS, Am, An, Bm, Bn);

    // Am_loc = mloc_A/*LDD_A*/; An_loc = nloc_A;
    // Bm_loc = mloc_B/*LDD_B*/; Bn_loc = nloc_B;
    Am_loc = LDD_A; An_loc = nloc_A;
    Bm_loc = LDD_B; Bn_loc = nloc_B;

    AGm = gM_A; AGn = gN_A;
    BGm = gM_B; BGn = gN_B;

    PASTE_CODE_INIT_LAPACK_MATRIX(dcA_lapack, two_dim_block_cyclic, A,
                                  (&dcA_lapack, matrix_RealDouble, matrix_Lapack,
                                   nodes_A, rank_A,
                                   MB_A, NB_A,
                                   AGm, AGn,
                                   cIA, cJA,
                                   Am, An,
                                   KP, KQ,
                                   iP_A, jQ_A,
                                   P_A,
                                   An_loc, Am_loc));

    PASTE_CODE_INIT_LAPACK_MATRIX(dcB_lapack, two_dim_block_cyclic, B,
                                  (&dcB_lapack, matrix_RealDouble, matrix_Lapack,
                                   nodes_B, rank_B,
                                   MB_B, NB_B,
                                   BGm, BGn,
                                   cIB, cJB,
                                   Bm, Bn,
                                   KP, KQ,
                                   iP_B, jQ_B,
                                   P_B,
                                   Bn_loc, Bm_loc));

#ifdef WRAPPER_VERBOSE
    //PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", *dcA, P_A, Q_A);
    PRINT(parsec_ctx, comm_B, PlasmaUpperLower, "dcB", *dcB, P_B, Q_B);
#endif

#ifdef MEASURE_INTERNAL_TIMES
    PASTE_CODE_FLOPS(FLOPS_DTRMM, (side, (DagDouble_t)*M, (DagDouble_t)*N))
#endif
    int P = P_A;
    int Q = Q_A;
    int gN = gN_A;
    int NB = NB_A;


    // int redisA = 1, redisB = 1;
    /* Redistribute if:
     * - Unaligned offsets with blocks
     * TODO check for tile compatibility and redistribute better
     */
    int redisA = 0, redisB = 0;
    if( (cIA % MB_A != 0) || ( cJA % NB_A != 0)) redisA = 1;
    if( (cIB % MB_B != 0) || ( cJB % NB_B != 0)) redisB = 1;

    if( (MB_A != MB_B)
        ||(NB_A != NB_B)) {
      /* redistribute all, there are internal calls with different block sizes */
        redisA = 1;
        redisB = 1;
        MB_A = MB_B = NB_A = NB_B = REDIS_BLOCKSZ;
    }

    PASTE_CODE_REDIS_INPUT(A, parsec_ctx, redisA, comm_A, BLOCKING_REDIS);
    PASTE_CODE_REDIS_INPUT(B, parsec_ctx, redisB, comm_B, DEFAULT_REDIS); /*We can't force BLOCKING and then set owners once datacopies created*/

    INI_WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec_ctx, dtrmm,
                              ( side, uplo,
                                trans, diag,
                                *ALPHA,
                               (parsec_tiled_matrix_dc_t *)dcA,
                               (parsec_tiled_matrix_dc_t *)dcB ),
                              dplasma_dtrmm_Destruct( PARSEC_dtrmm ),
                              rank_A, P_A, Q_A, NB_A, gN_A, comm_A);

#ifndef DO_BLOCKING_REDISTRIBUTION
#ifdef ALL_REDIS_NONBLOCKING
    parsec_data_collection_set_owner((parsec_data_collection_t *)dcA, PARSEC_dtrmm);
#endif
    parsec_data_collection_set_owner((parsec_data_collection_t *)dcB, PARSEC_dtrmm);
#endif
    PASTE_CODE_REDIS_OUTPUT_INI(B, parsec_ctx, redisB, comm_B, DEFAULT_REDIS);

    FINI_WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec_ctx, dtrmm,
                              ( side, uplo,
                                trans, diag,
                                *ALPHA,
                               (parsec_tiled_matrix_dc_t *)dcA,
                               (parsec_tiled_matrix_dc_t *)dcB ),
                              dplasma_dtrmm_Destruct( PARSEC_dtrmm ),
                              rank_A, P_A, Q_A, NB_A, gN_A, comm_A);

    PASTE_CODE_CLEANUP_REDIS(A, parsec_ctx, redisA, comm_A);
    PASTE_CODE_REDIS_OUTPUT_FINI(B, parsec_ctx, redisB, comm_B, DEFAULT_REDIS);
    PASTE_CODE_CLEANUP_REDIS(B, parsec_ctx, redisB, comm_B);

#ifdef WRAPPER_VERBOSE
    //PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", *dcA, P_A, Q_A);
    PRINT(parsec_ctx, comm_B, PlasmaUpperLower, "dcB", *dcB, P_B, Q_B);
#endif

#ifdef CHECK_RESULTS
    {
      int check=1;
      int loud=5;
      int rank;
      MPI_Comm_rank(comm_A, &rank);
      /* Only our checking */
      int Aseed = 3872;
      int Bseed = 2873;

      /* Check the factorization */
      PASTE_CODE_ALLOCATE_MATRIX(dcB_out, check,
          two_dim_block_cyclic, (&dcB_out, matrix_RealDouble, matrix_Tile,
                  nodes_B, rank_B, MB_B, NB_B, BGm, BGn, 0, 0,
                  Bm, Bn, KP, KQ, 0, 0, P_B));
      dcopy_lapack_tile(parsec_ctx, dcB, &dcB_out, mloc_B, nloc_B);

      int info_solution = check_solution(parsec_ctx, rank == 0 ? loud : 0,
                                         side, uplo, trans, diag,
                                         *ALPHA, Am, An, Aseed,
                                                 *M, *N, Bseed,
                                         &dcB_out);

      if ( rank == 0 ) {
          if (info_solution == 0) {
              printf(" ---- TESTING ZTRMM (%c, %c, %c, %c) ...... PASSED !\n",
                     *SIDE, *UPLO, *TRANS, *DIAG);
          }
          else {
              printf(" ---- TESTING ZTRMM (%c, %c, %c, %c) ... FAILED !\n",
                     *SIDE, *UPLO, *TRANS, *DIAG);
          }
          printf("***************************************************\n");
      }
      /* Cleanup */
      parsec_data_free(dcB_out.mat);
      parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcB_out);
    }
#endif

    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)dcA);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)dcB);

}

#ifdef CHECK_RESULTS
static int check_solution( parsec_context_t *parsec, int loud,
                           dplasma_enum_t side, dplasma_enum_t uplo, dplasma_enum_t trans, dplasma_enum_t diag,
                           double alpha,
                           int Am, int An, int Aseed,
                           int M,  int N,  int Cseed,
                           two_dim_block_cyclic_t *dcCfinal )
{
    int info_solution = 1;
    double Anorm, Cinitnorm, Cdplasmanorm, Clapacknorm, Rnorm;
    double eps, result;
    int MB = dcCfinal->super.mb;
    int NB = dcCfinal->super.nb;
    int LDA = Am;
    int LDC = M;
    int rank  = dcCfinal->super.super.myrank;

    eps = LAPACKE_dlamch_work('e');

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        two_dim_block_cyclic, (&dcA, matrix_RealDouble, matrix_Lapack,
                               1, rank, MB, NB, LDA, An, 0, 0,
                               Am, An, 1, 1, 0, 0, 1));
    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        two_dim_block_cyclic, (&dcC, matrix_RealDouble, matrix_Lapack,
                               1, rank, MB, NB, LDC, N, 0, 0,
                               M, N, 1, 1, 0, 0, 1));

    dplasma_dplrnt( parsec, 0, (parsec_tiled_matrix_dc_t *)&dcA, Aseed);
    dplasma_dplrnt( parsec, 0, (parsec_tiled_matrix_dc_t *)&dcC, Cseed );

    Anorm        = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)&dcA );
    Cinitnorm    = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)&dcC );
    Cdplasmanorm = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)dcCfinal );

    if ( rank == 0 ) {
        cblas_dtrmm(CblasColMajor,
                    (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                    (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                    M, N,
                    (alpha), dcA.mat, LDA,
                                        dcC.mat, LDC );
    }

    Clapacknorm = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)&dcC );

    dplasma_dgeadd( parsec, dplasmaNoTrans,
                    -1.0, (parsec_tiled_matrix_dc_t*)dcCfinal,
                     1.0, (parsec_tiled_matrix_dc_t*)&dcC );

    Rnorm = dplasma_dlange( parsec, dplasmaMaxNorm, (parsec_tiled_matrix_dc_t*)&dcC );

    result = Rnorm / (Clapacknorm * max(M,N) * eps);

    if ( rank == 0 ) {
        if ( loud > 2 ) {
            printf("  ||A||_inf = %e, ||C||_inf = %e\n"
                   "  ||lapack(a*A*C)||_inf = %e, ||dplasma(a*A*C)||_inf = %e, ||R||_m = %e, res = %e\n",
                   Anorm, Cinitnorm, Clapacknorm, Cdplasmanorm, Rnorm, result);
        }

        if (  isinf(Clapacknorm) || isinf(Cdplasmanorm) ||
              isnan(result) || isinf(result) || (result > 10.0) ) {
            info_solution = 1;
        }
        else {
            info_solution = 0;
        }
    }

#if defined(PARSEC_HAVE_MPI)
    MPI_Bcast(&info_solution, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    parsec_data_free(dcA.mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA);
    parsec_data_free(dcC.mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcC);

    return info_solution;
}
#endif

GENERATE_F77_BINDINGS (PDTRMM,
                       pdtrmm,
                       pdtrmm_,
                       pdtrmm__,
                       pdtrmm_w,
                       ( char* SIDE, char* UPLO, char* TRANS, char* DIAG, int * M, int * N, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB ),
                       ( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB ))
