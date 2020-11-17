#include "common.h"

/*
 *  Purpose
 *  =======
 *
 *  PDGEMM  performs one of the matrix-matrix operations
 *
 *     sub( C ) := alpha*op( sub( A ) )*op( sub( B ) ) + beta*sub( C ),
 *
 *  where
 *
 *     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),  and, op( X )  is one  of
 *     op( X ) = X   or   op( X ) = X'.
 *
 *  Thus, op( sub( A ) ) denotes A(IA:IA+M-1,JA:JA+K-1)  if TRANSA = 'N',
 *                               A(IA:IA+K-1,JA:JA+M-1)' if TRANSA = 'T',
 *                               A(IA:IA+K-1,JA:JA+M-1)' if TRANSA = 'C',
 *
 *  and,  op( sub( B ) ) denotes B(IB:IB+K-1,JB:JB+N-1)  if TRANSB = 'N',
 *                               B(IB:IB+N-1,JB:JB+K-1)' if TRANSB = 'T',
 *                               B(IB:IB+N-1,JB:JB+K-1)' if TRANSB = 'C',
 *
 *  Alpha and beta are scalars.  A, B and C are matrices;  op( sub( A ) )
 *  is an  m by k submatrix,  op( sub( B ) )  is an  k by n submatrix and
 *  sub( C ) is an m by n submatrix.
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
 *  TRANSB  (global input) CHARACTER*1
 *          On entry,  TRANSB  specifies the form of op( sub( B ) ) to be
 *          used in the matrix multiplication as follows:
 *
 *             TRANSB = 'N' or 'n'   op( sub( B ) ) = sub( B ),
 *
 *             TRANSB = 'T' or 't'   op( sub( B ) ) = sub( B )',
 *
 *             TRANSB = 'C' or 'c'   op( sub( B ) ) = sub( B )'.
 *
 *  M       (global input) INTEGER
 *          On entry,  M  specifies  the number of rows of the  submatrix
 *          op( sub( A ) ) and of the submatrix sub( C ). M  must  be  at
 *          least  zero.
 *
 *  N       (global input) INTEGER
 *          On entry, N specifies the number of columns of the  submatrix
 *          op( sub( B ) )  and  the  number of columns of the  submatrix
 *          sub( C ). N must be at least zero.
 *
 *  K       (global input) INTEGER
 *          On entry, K specifies the number of columns of the  submatrix
 *          op( sub( A ) )  and  the  number of rows   of  the  submatrix
 *          op( sub( B ) ). K must be at least  zero.
 *
 *  ALPHA   (global input) DOUBLE PRECISION
 *          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *          supplied  as zero then the local entries of the arrays  A and
 *          B corresponding to the entries of  the  submatrices  sub( A )
 *          and sub( B ) respectively need not be set on input.
 *
 *  A       (local input) DOUBLE PRECISION array
 *          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
 *          at least Lc( 1, JA+K-1 ) when  TRANSA = 'N' or 'n', and is at
 *          least  Lc( 1, JA+M-1 )  otherwise.  Before  entry, this array
 *          contains the local entries of the matrix A.
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
 *  B       (local input) DOUBLE PRECISION array
 *          On entry, B is an array of dimension (LLD_B, Kb), where Kb is
 *          at least Lc( 1, JB+N-1 ) when  TRANSB = 'N' or 'n', and is at
 *          least Lc( 1, JB+K-1 )  otherwise.  Before  entry,  this array
 *          contains the local entries of the matrix B.
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
 *  BETA    (global input) DOUBLE PRECISION
 *          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
 *          supplied  as  zero  then  the  local entries of  the array  C
 *          corresponding to  the  entries of the submatrix sub( C ) need
 *          not be set on input.
 *
 *  C       (local input/local output) DOUBLE PRECISION array
 *          On entry, C is an array of dimension (LLD_C, Kc), where Kc is
 *          at least Lc( 1, JC+N-1 ).  Before  entry, this array contains
 *          the local entries of the matrix  C.
 *          On exit, the entries of this array corresponding to the local
 *          entries of the  submatrix  sub( C )  are  overwritten  by the
 *          local entries of the m by n updated submatrix.
 *
 *  IC      (global input) INTEGER
 *          On entry, IC  specifies C's global row index, which points to
 *          the beginning of the submatrix sub( C ).
 *
 *  JC      (global input) INTEGER
 *          On entry, JC  specifies C's global column index, which points
 *          to the beginning of the submatrix sub( C ).
 *
 *  DESCC   (global and local input) INTEGER array
 *          On entry, DESCC  is an integer array of dimension DLEN_. This
 *          is the array descriptor for the matrix C.
 *
 *  -- Written on April 1, 1998 by
 *     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
 *
 *  ---------------------------------------------------------------------
 */

#ifdef CHECK_RESULTS
static int check_solution( parsec_context_t *parsec, int loud,
                           dplasma_enum_t transA, dplasma_enum_t transB,
                           double alpha, int Am, int An, int Aseed,
                                         int Bm, int Bn, int Bseed,
                           double beta,  int M,  int N,  int Cseed,
                           two_dim_block_cyclic_t *dcCfinal );
#endif


/* pdgemmwrap_ for fine grain control of what is wrapped */
// void pdgemmwrap_
void pdgemm_w( char* TRANSA, char* TRANSB,
               int * M, int * N, int * K,
               double * ALPHA,
               double * A, int * IA, int * JA, int * DESCA,
               double * B, int * IB, int * JB, int * DESCB,
               double * BETA,
               double * C, int * IC, int * JC, int * DESCC ){

#ifdef COUNT_WRAPPED_CALLS
    count_PDGEMM++;
#endif

    if( (*M == 0) || (*N == 0) || (*K == 0)){
      /* NOP */
      return;
    }

    int KP = 1;
    int KQ = 1;

    PASTE_SETUP(A);
    PASTE_SETUP(B);
    PASTE_SETUP(C);

#ifdef WRAPPER_VERBOSE_CALLS
    if(rank_A == 0){
      printf("V-PDGEMM M%d N%d K%d "
             "IA%d JA%d A%p MBA%d NBA%d %c "
             "IB%d JB%d B%p MBB%d NBB%d %c "
             "IC%d JC%d C%p MBC%d NBC%d \n",
             *M, *N, *K,
             *IA, *JA, A, DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_], *TRANSA,
             *IB, *JB, B, DESCB[WRAPPER_MB1_], DESCB[WRAPPER_NB1_], *TRANSB,
             *IC, *JC, C, DESCC[WRAPPER_MB1_], DESCC[WRAPPER_NB1_]);
    }
#endif

    dplasma_enum_t tA = OP_TRANS(*TRANSA);
    dplasma_enum_t tB = OP_TRANS(*TRANSB);

    /*trust dplasma precisionn generator*/
    if(tA == dplasmaConjTrans) tA=dplasmaTrans;
    if(tB == dplasmaConjTrans) tB=dplasmaTrans;

    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  "M%d N%d K%d IA%d JA%d (ictxt)DESCA[WRAPPER_CTXT1_] %d, "
          "(gM)DESCA[WRAPPER_M1_] %d, (gN)DESCA[WRAPPER_N1_] %d, (MB)DESCA[WRAPPER_MB1_] %d, (NB)DESCA[WRAPPER_NB1_] %d, "
          "DESCA[WRAPPER_RSRC1_] %d, DESCA[WRAPPER_CSRC1_] %d, (LLD)DESCA[WRAPPER_LLD1_] %d TRANS%s (%c)",
          *M, *N, *K, *IA, *JA, DESCA[WRAPPER_CTXT1_],
          DESCA[WRAPPER_M1_], DESCA[WRAPPER_N1_], DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_],
          DESCA[WRAPPER_RSRC1_], DESCA[WRAPPER_CSRC1_], DESCA[WRAPPER_LLD1_], (tA==PlasmaNoTrans)? "PlasmaNoTrans": "PlasmaTrans", *TRANSA);

    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  "M%d N%d K%d IB%d JB%d (ictxt)DESCB[WRAPPER_CTXT1_] %d, "
          "(gM)DESCB[WRAPPER_M1_] %d, (gN)DESCB[WRAPPER_N1_] %d, (MB)DESCB[WRAPPER_MB1_] %d, (NB)DESCB[WRAPPER_NB1_] %d, "
          "DESCB[WRAPPER_RSRC1_] %d, DESCB[WRAPPER_CSRC1_] %d, (LLD)DESCB[WRAPPER_LLD1_] %d TRANS%s (%c)",
          *M, *N, *K, *IB, *JB, DESCB[WRAPPER_CTXT1_],
          DESCB[WRAPPER_M1_], DESCB[WRAPPER_N1_], DESCB[WRAPPER_MB1_], DESCB[WRAPPER_NB1_],
          DESCB[WRAPPER_RSRC1_], DESCB[WRAPPER_CSRC1_], DESCB[WRAPPER_LLD1_], (tB==PlasmaNoTrans)? "PlasmaNoTrans": "PlasmaTrans", *TRANSB);

    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  "M%d N%d K%d IC%d JC%d (ictxt)DESCC[WRAPPER_CTXT1_] %d, "
          "(gM)DESCC[WRAPPER_M1_] %d, (gN)DESCC[WRAPPER_N1_] %d, (MB)DESCC[WRAPPER_MB1_] %d, (NB)DESCC[WRAPPER_NB1_] %d, "
          "DESCC[WRAPPER_RSRC1_] %d, DESCC[WRAPPER_CSRC1_] %d, (LLD)DESCC[WRAPPER_LLD1_] %d",
          *M, *N, *K, *IC, *JC, DESCC[WRAPPER_CTXT1_],
          DESCC[WRAPPER_M1_], DESCC[WRAPPER_N1_], DESCC[WRAPPER_MB1_], DESCC[WRAPPER_NB1_],
          DESCC[WRAPPER_RSRC1_], DESCC[WRAPPER_CSRC1_], DESCC[WRAPPER_LLD1_]);

    assert(comm_index_A == comm_index_B);
    assert(comm_index_A == comm_index_C);

    parsec_init_wrapped_call((void*)comm_A);

//    A: M x K
//    B: K x N
//    C: M x N

    int Am, An, Bm, Bn, Cm, Cn;
    // int Am_loc, An_loc, Bm_loc, Bn_loc, Cm_loc, Cn_loc;
    // int AGm, AGn, BGm, BGn, CGm, CGn;

    /*We are expecting the original dim (without transposing)*/
    if ( tA == PlasmaNoTrans ) {
        Am = *M; An = *K;
    } else {
        Am = *K; An = *M;
    }
    if ( tB == PlasmaNoTrans ) {
        Bm = *K; Bn = *N;
    } else {
        Bm = *N; Bn = *K;
    }
    Cm = *M; Cn = *N;

    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,
        "A-%c %dx%d B-%c %dx%d C %dx%d",
        *TRANSA, Am, An, *TRANSB, Bm, Bn, Cn, Cm);

    two_dim_block_cyclic_t dcA_lapack;
    two_dim_block_cyclic_lapack_init(&dcA_lapack, matrix_RealDouble, matrix_Lapack,
                                     rank_A,
                                     MB_A, NB_A,
                                     gM_A, gN_A,
                                     cIA, cJA,
                                     Am, An,
                                     P_A, Q_A,
                                     KP, KQ,
                                     iP_A, jQ_A,
                                     LLD_A, nloc_A);
    dcA_lapack.mat = A;
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcA_lapack, "dcA_lapack");

    two_dim_block_cyclic_t dcB_lapack;
    two_dim_block_cyclic_lapack_init(&dcB_lapack, matrix_RealDouble, matrix_Lapack,
                                     rank_B,
                                     MB_B, NB_B,
                                     gM_B, gN_B,
                                     cIB, cJB,
                                     Bm, Bn,
                                     P_B, Q_B,
                                     KP, KQ,
                                     iP_B, jQ_B,
                                     LLD_B, nloc_B);
    dcB_lapack.mat = B;
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcB_lapack, "dcB_lapack");

    two_dim_block_cyclic_t dcC_lapack;
    two_dim_block_cyclic_lapack_init(&dcC_lapack, matrix_RealDouble, matrix_Lapack,
                                     rank_C,
                                     MB_C, NB_C,
                                     gM_C, gN_C,
                                     cIC, cJC,
                                     Cm, Cn,
                                     P_C, Q_C,
                                     KP, KQ,
                                     iP_C, jQ_C,
                                     LLD_C, nloc_C);
    dcC_lapack.mat = C;
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcC_lapack, "dcC_lapack");

    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", (&dcA_lapack));
    PRINT(parsec_ctx, comm_B, PlasmaUpperLower, "dcB", (&dcB_lapack));
    PRINT(parsec_ctx, comm_C, PlasmaUpperLower, "dcC", (&dcC_lapack));

    /* Redistribute if:
     * - Unaligned offsets with blocks.
     * - Different block sizes.
     */
    int redisA = 0, redisB = 0, redisC = 0;
    /* Non aligned offsets? */
    redisA = ( (cIA % MB_A != 0) || ( cJA % NB_A != 0) );
    redisB = ( (cIB % MB_B != 0) || ( cJB % NB_B != 0) );
    redisC = ( (cIC % MB_C != 0) || ( cJC % NB_C != 0) );

    /* Different block sizes? */
    if( (MB_A != MB_B) || (MB_A != MB_C) || (NB_A != NB_B) || (NB_A != NB_C) ) {
        /* redistribute all, different internal block sizes */
        redisA = redisB = redisC = 1;
    }

    /* If redistributing one matrix, redistribute all
     * TODO optimization: check for tile compatibility and avoid redistributions?
     */
    if( redisA || redisB || redisC ) {
        /* redistribute all, different internal block sizes */
        redisA = redisB = redisC = 1;
    }

    two_dim_block_cyclic_t *dcA = redistribute_lapack_input(&dcA_lapack, redisA, comm_A, rank_A, "redisA");
    two_dim_block_cyclic_t *dcB = redistribute_lapack_input(&dcB_lapack, redisB, comm_B, rank_B, "redisB");
    two_dim_block_cyclic_t *dcC = redistribute_lapack_input(&dcC_lapack, redisC, comm_C, rank_C, "redisC");

#ifdef MEASURE_INTERNAL_TIMES
    PASTE_CODE_FLOPS(FLOPS_DGEMM, ((DagDouble_t)*M,(DagDouble_t)*N,(DagDouble_t)*K));
#endif

    WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec_ctx, dgemm,
                              (tA, tB, *ALPHA,
                               (parsec_tiled_matrix_dc_t *)dcA,
                               (parsec_tiled_matrix_dc_t *)dcB,
                               *BETA,
                               (parsec_tiled_matrix_dc_t *)dcC),
                              dplasma_dgemm_Destruct( PARSEC_dgemm ),
                              rank_A, P_A, Q_A, NB_A, gN_A, comm_A);

    dcA = redistribute_lapack_output_cleanup(&dcA_lapack, dcA, 0, comm_A, rank_A, "redisA");
    dcB = redistribute_lapack_output_cleanup(&dcB_lapack, dcB, 0, comm_B, rank_B, "redisB");
    dcC = redistribute_lapack_output_cleanup(&dcC_lapack, dcC, 1, comm_C, rank_C, "redisC");

    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", dcA);
    PRINT(parsec_ctx, comm_B, PlasmaUpperLower, "dcB", dcB);
    PRINT(parsec_ctx, comm_C, PlasmaUpperLower, "dcC", dcC);

#ifdef CHECK_RESULTS
    /* Only for our testing, checking using same seed */
    int check = 1;
    int loud = 5;
    /* Check the factorization */
    PASTE_CODE_ALLOCATE_MATRIX(dcC_out, check,
        two_dim_block_cyclic, (&dcC_out, matrix_RealDouble, matrix_Tile,
                rank_C, MB_C, NB_C, gM_C, gN_C, 0, 0,
                Cm, Cn, P_C, Q_C,  KP, KQ, 0, 0));
    dcopy_lapack_tile(parsec_ctx, dcC, &dcC_out, mloc_C, nloc_C);

    int Aseed = 3872;
    int Bseed = 4674;
    int Cseed = 2873;
    /* Check the solution */
    int info_solution = check_solution( parsec_ctx, (rank_A == 0) ? loud : 0,
                                        tA, tB,
                                        *ALPHA, Am, An, Aseed,
                                                Bm, Bn, Bseed,
                                        *BETA,  *M, *N, Cseed,
                                        &dcC_out);

    if ( rank_A == 0 ) {
        if (info_solution == 0) {
            printf(" ---- TESTING GEMM (%c, %c) ...... PASSED !\n",
                   *TRANSA, *TRANSB);
        }
        else {
            printf(" ---- TESTING GEMM (%c, %c) ... FAILED !\n",
                   *TRANSA, *TRANSB);
        }
        printf("***************************************************\n");
    }
    /* Cleanup */
    parsec_data_free(dcC_out.mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcC_out);
#endif

    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)dcA );
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)dcB );
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)dcC );

}

#ifdef CHECK_RESULTS
static int check_solution( parsec_context_t *parsec, int loud,
                           dplasma_enum_t transA, dplasma_enum_t transB,
                           double alpha, int Am, int An, int Aseed,
                                         int Bm, int Bn, int Bseed,
                           double beta,  int M,  int N,  int Cseed,
                           two_dim_block_cyclic_t *dcCfinal )
{
    int info_solution = 1;
    double Anorm, Bnorm, Cinitnorm, Cdplasmanorm, Clapacknorm, Rnorm;
    double eps, result;
    int K  = ( transA == dplasmaNoTrans ) ? An : Am ;
    int MB = dcCfinal->super.mb;
    int NB = dcCfinal->super.nb;
    int LDA = Am;
    int LDB = Bm;
    int LDC = M;
    int rank  = dcCfinal->super.super.myrank;

    eps = LAPACKE_dlamch_work('e');

    PASTE_CODE_ALLOCATE_MATRIX(dcA, 1,
        two_dim_block_cyclic, (&dcA, matrix_RealDouble, matrix_Lapack,
                               rank, MB, NB, LDA, An, 0, 0,
                               Am, An, 1, 1, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcB, 1,
        two_dim_block_cyclic, (&dcB, matrix_RealDouble, matrix_Lapack,
                               rank, MB, NB, LDB, Bn, 0, 0,
                               Bm, Bn, 1, 1, 1, 1, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcC, 1,
        two_dim_block_cyclic, (&dcC, matrix_RealDouble, matrix_Lapack,
                               rank, MB, NB, LDC, N, 0, 0,
                               M, N, 1, 1, 1, 1, 0, 0));

    dplasma_dplrnt( parsec, 0, (parsec_tiled_matrix_dc_t *)&dcA, Aseed );
    dplasma_dplrnt( parsec, 0, (parsec_tiled_matrix_dc_t *)&dcB, Bseed );
    dplasma_dplrnt( parsec, 0, (parsec_tiled_matrix_dc_t *)&dcC, Cseed );

    Anorm        = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)&dcA );
    Bnorm        = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)&dcB );
    Cinitnorm    = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)&dcC );
    Cdplasmanorm = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)dcCfinal );

    if ( rank == 0 ) {
        cblas_dgemm(CblasColMajor,
                    (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                    M, N, K,
                    (alpha), dcA.mat, LDA,
                                        dcB.mat, LDB,
                    (beta),  dcC.mat, LDC );
    }

    Clapacknorm = dplasma_dlange( parsec, dplasmaInfNorm, (parsec_tiled_matrix_dc_t*)&dcC );

    dplasma_dgeadd( parsec, dplasmaNoTrans, -1.0, (parsec_tiled_matrix_dc_t*)dcCfinal,
                                           1.0, (parsec_tiled_matrix_dc_t*)&dcC );

    Rnorm = dplasma_dlange( parsec, dplasmaMaxNorm, (parsec_tiled_matrix_dc_t*)&dcC);

    if ( rank == 0 ) {
        if ( loud > 2 ) {
            printf("  ||A||_inf = %e, ||B||_inf = %e, ||C||_inf = %e\n"
                   "  ||lapack(a*A*B+b*C)||_inf = %e, ||dplasma(a*A*B+b*C)||_inf = %e, ||R||_m = %e\n",
                   Anorm, Bnorm, Cinitnorm, Clapacknorm, Cdplasmanorm, Rnorm);
        }

        result = Rnorm / ((Anorm + Bnorm + Cinitnorm) * max(M,N) * eps);
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
    parsec_data_free(dcB.mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcB);
    parsec_data_free(dcC.mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcC);

    return info_solution;
}

#endif

GENERATE_F77_BINDINGS (PDGEMM,
                       pdgemm,
                       pdgemm_,
                       pdgemm__,
                       pdgemm_w,
                       ( char* TRANSA, char* TRANSB, int * M, int * N, int * K, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB, double * BETA, double * C, int * IC, int * JC, int * DESCC ),
                       ( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC ))
