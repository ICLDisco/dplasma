#include "common.h"
  /************************************************************************
    DPLASMA GETRF_1D doesn't support 2D block cyclic.
    We redistribute to 1D when lapack matrix A is 2DBC.
    Wrapper uses its own temporary IPIV in TILED format.
  ************************************************************************/

/*     SUBROUTINE PDGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
 *
 *  -- ScaLAPACK routine (version 1.7) --
 *     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
 *     and University of California, Berkeley.
 *     May 25, 2001
 *
 *     .. Scalar Arguments ..
       INTEGER            IA, INFO, JA, M, N
 *     ..
 *     .. Array Arguments ..
       INTEGER            DESCA( * ), IPIV( * )
       DOUBLE PRECISION   A( * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  PDGETRF computes an LU factorization of a general M-by-N distributed
 *  matrix sub( A ) = (IA:IA+M-1,JA:JA+N-1) using partial pivoting with
 *  row interchanges.
 *
 *  The factorization has the form sub( A ) = P * L * U, where P is a
 *  permutation matrix, L is lower triangular with unit diagonal ele-
 *  ments (lower trapezoidal if m > n), and U is upper triangular
 *  (upper trapezoidal if m < n). L and U are stored in sub( A ).
 *
 *  This is the right-looking Parallel Level 3 BLAS version of the
 *  algorithm.
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
 *  M       (global input) INTEGER
 *          The number of rows to be operated on, i.e. the number of rows
 *          of the distributed submatrix sub( A ). M >= 0.
 *
 *  N       (global input) INTEGER
 *          The number of columns to be operated on, i.e. the number of
 *          columns of the distributed submatrix sub( A ). N >= 0.
 *
 *  A       (local input/local output) DOUBLE PRECISION pointer into the
 *          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
 *          On entry, this array contains the local pieces of the M-by-N
 *          distributed matrix sub( A ) to be factored. On exit, this
 *          array contains the local pieces of the factors L and U from
 *          the factorization sub( A ) = P*L*U; the unit diagonal ele-
 *          ments of L are not stored.
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
 *  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
 *          This array contains the pivoting information.
 *          IPIV(i) -> The global row local row i was swapped with.
 *          This array is tied to the distributed matrix A.
 *
 *  INFO    (global output) INTEGER
 *          = 0:  successful exit
 *          < 0:  If the i-th argument is an array and the j-entry had
 *                an illegal value, then INFO = -(i*100+j), if the i-th
 *                argument is a scalar and had an illegal value, then
 *                INFO = -i.
 *          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
 *                The factorization has been completed, but the factor U
 *                is exactly singular, and division by zero will occur if
 *                it is used to solve a system of equations.
 *
 *  =====================================================================
 */

#ifdef CHECK_RESULTS
static int check_solution( parsec_context_t *parsec, int loud,
                          parsec_tiled_matrix_t *dcA,
                          parsec_tiled_matrix_t *dcB,
                          parsec_tiled_matrix_t *dcX );

static int check_inverse( parsec_context_t *parsec, int loud,
                         parsec_tiled_matrix_t *dcA,
                         parsec_tiled_matrix_t *dcInvA,
                         parsec_tiled_matrix_t *dcI );
#endif

void pdgetrf_w(int * M,
              int * N,
              double * A,
              int * IA,
              int * JA,
              int * DESCA,
              int * IPIV,
              int * info){


#ifdef COUNT_WRAPPED_CALLS
    count_PDGETRF_1D++;
#endif
    *info=0;
    if( (*M == 0) || (*N == 0)){
      /* NOP */
      return;
    }

    /* Current DPLASMA GETRF_1D doesn't support 2D block cyclic.
     * We redistribute to 1D when lapack matrix A is 2DBC.
     * Wrapper uses its own temporary IPIV in TILED format.
     */

    int KP = 1;
    int KQ = 1;
    PASTE_SETUP(A);

#ifdef WRAPPER_VERBOSE_CALLS
    if(rank_A == 0){
      printf("V-PDGETRF_1D M%d N%d "
             "IA%d JA%d A%p MBA%d NBA%d \n",
             *M, *N,
             *IA, *JA, A, DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_]);
    }
#endif

    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  "M%d N%d IA%d JA%d (ictxt)DESCA[WRAPPER_CTXT1_] %d, "
          "(gM)DESCA[WRAPPER_M1_] %d, (gN)DESCA[WRAPPER_N1_] %d, (MB)DESCA[WRAPPER_MB1_] %d, (NB)DESCA[WRAPPER_NB1_] %d, "
          "DESCA[WRAPPER_RSRC1_] %d, DESCA[WRAPPER_CSRC1_] %d, (LLD)DESCA[WRAPPER_LLD1_] %d\n",
          *M, *N, *IA, *JA, DESCA[WRAPPER_CTXT1_],
          DESCA[WRAPPER_M1_], DESCA[WRAPPER_N1_], DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_],
          DESCA[WRAPPER_RSRC1_], DESCA[WRAPPER_CSRC1_], DESCA[WRAPPER_LLD1_]);

    parsec_init_wrapped_call(comm_A);

    parsec_matrix_block_cyclic_t dcA_lapack;
    parsec_matrix_block_cyclic_lapack_init(&dcA_lapack, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_LAPACK,
                                     rank_A,
                                     MB_A, NB_A,
                                     gM_A, gN_A,
                                     cIA, cJA,
                                     *M, *N,
                                     P_A, Q_A,
                                     KP, KQ,
                                     iP_A, jQ_A,
                                     LLD_A, nloc_A);
    dcA_lapack.mat = A;
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcA_lapack, "dcA_lapack");

    /* SCALAPACK:
     *  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
     *          This array contains the pivoting information.
     *          IPIV(i) -> The global row local row i was swapped with.
     *          This array is tied to the distributed matrix A.
     *  However, in the original ScaLAPACK GETRF every process has the full
     *  IPIV containing global row indexes, therefore, the wrapper needs to rebuild it.
     *
     * DPLASMA:
     * @param[out] IPIV
     *          Descriptor of the IPIV matrix. Should be of size 1-by-min(M,N).
     *          On exit, contains the pivot indices; for 1 <= i <= min(M,N), row i
     *          of the matrix was interchanged with row IPIV(i).
     */

    int redisA = 0;
    if( (cIA % MB_A != 0) || ( cJA % NB_A != 0)) redisA = 1;

    int redisP = 1;
    int redisQ = dcA_lapack.grid.rows*dcA_lapack.grid.cols;
    parsec_matrix_block_cyclic_t *dcA = redistribute_lapack_input_1D(&dcA_lapack, redisA, comm_A, rank_A, "redisA", redisP, redisQ);
    /* Matrix A is only redistributed to lapack if redisA or P_A != 1*/

    int redisMB = dcA->super.mb;
    (void)redisMB;
    int redisNB = dcA->super.nb;
    int gN_IPIV = dplasma_imin(*M, *N);

    PASTE_CODE_ALLOCATE_MATRIX(dcIPIV_tmp, 1,
                               parsec_matrix_block_cyclic,
                               (&dcIPIV_tmp, PARSEC_MATRIX_INTEGER, PARSEC_MATRIX_TILE,
                                rank_A,
                                1, redisNB,
                                1, gN_IPIV,
                                0, 0,
                                1, gN_IPIV,
                                redisP, redisQ,
                                KP, KQ,
                                0, 0));

    parsec_matrix_block_cyclic_t *dcIPIV = &dcIPIV_tmp;

    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", dcA);
    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcIPIV", dcIPIV);

#ifdef CHECK_RESULTS
    int check = 1;
    int check_inv = 1;

    if ( *M != *N && check ) {
        fprintf(stderr, "Check is impossible if M != N\n");
        check = 0;
    }

    int cLDA = *M;
    int NRHS  = KP;
    int LDB   = *M;

    PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
                               parsec_matrix_block_cyclic, (&dcA0, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank_A, redisMB, redisNB, cLDA, *N, 0, 0,
                                                      *M, *N, redisP, redisQ, KP, KQ, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcA_out, check,
                               parsec_matrix_block_cyclic, (&dcA_out, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank_A, redisMB, redisNB, cLDA, *N, 0, 0,
                                                      *M, *N, redisP, redisQ, KP, KQ, 0, 0));
    /* Random B check */
    PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
                               parsec_matrix_block_cyclic, (&dcB, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank_A, redisMB, redisNB, LDB, NRHS, 0, 0,
                                                      *M, NRHS, redisP, redisQ, KP, KQ, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
                               parsec_matrix_block_cyclic, (&dcX, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank_A, redisMB, redisNB, LDB, NRHS, 0, 0,
                                                      *M, NRHS, redisP, redisQ, KP, KQ, 0, 0));
    /* Inverse check */
    PASTE_CODE_ALLOCATE_MATRIX(dcInvA, check_inv,
                               parsec_matrix_block_cyclic, (&dcInvA, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank_A, redisMB, redisNB, cLDA, *N, 0, 0,
                                                      *M, *N, redisP, redisQ, KP, KQ, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcI, check_inv,
                               parsec_matrix_block_cyclic, (&dcI, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                                      rank_A, redisMB, redisNB, cLDA, *N, 0, 0,
                                                      *M, *N, redisP, redisQ, KP, KQ, 0, 0));


    if ( check ) {
        if(dcA == &dcA_lapack){
            dcopy_lapack_tile(parsec_ctx, dcA, &dcA0, mloc_A, nloc_A);
        }else{
            dplasma_dlacpy( parsec_ctx, PlasmaUpperLower,
                            (parsec_tiled_matrix_t *)dcA,
                            (parsec_tiled_matrix_t *)&dcA0 );
        }

        dplasma_dplrnt( parsec_ctx, 0, (parsec_tiled_matrix_t *)&dcB, 2354);
        dplasma_dlacpy( parsec_ctx, PlasmaUpperLower,
                        (parsec_tiled_matrix_t *)&dcB,
                        (parsec_tiled_matrix_t *)&dcX );
    }
    if ( check_inv ) {
        dplasma_dlaset( parsec_ctx, PlasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&dcI);
        dplasma_dlaset( parsec_ctx, PlasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&dcInvA);
    }

#endif


#ifdef MEASURE_INTERNAL_TIMES
    PASTE_CODE_FLOPS(FLOPS_DGETRF, ((DagDouble_t)*M,(DagDouble_t)*N));
#endif
    WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec_ctx, dgetrf_1d,
                          ((parsec_tiled_matrix_t*)dcA,
                           (parsec_tiled_matrix_t*)dcIPIV, info),
                          dplasma_dgetrf_1d_Destruct( PARSEC_dgetrf_1d ),
                          rank_A, redisP, redisQ, NB_A, gN_A, comm_A);

#ifdef CHECK_RESULTS
    /* check before redistributing back to lapack format */
    int loud=5;
    if ( check ) {
        if(dcA == &dcA_lapack){
            dcopy_lapack_tile(parsec_ctx, dcA, &dcA_out, mloc_A, nloc_A);
        }else{
            dplasma_dlacpy( parsec_ctx, PlasmaUpperLower,
                            (parsec_tiled_matrix_t *)dcA,
                            (parsec_tiled_matrix_t *)&dcA_out );
        }
        /*
         * First check with a right hand side
         */
        dplasma_dgetrs(parsec_ctx, PlasmaNoTrans,
                       (parsec_tiled_matrix_t *)&dcA_out,
                       (parsec_tiled_matrix_t *)dcIPIV,
                       (parsec_tiled_matrix_t *)&dcX );

        /* Check the solution */
        check_solution( parsec_ctx, (rank_A == 0) ? loud : 0,
                               (parsec_tiled_matrix_t *)&dcA0,
                               (parsec_tiled_matrix_t *)&dcB,
                               (parsec_tiled_matrix_t *)&dcX);

        /*
         * Second check with inverse
         */
        if ( check_inv ) {
            dplasma_dgetrs(parsec_ctx, PlasmaNoTrans,
                           (parsec_tiled_matrix_t *)&dcA_out,
                           (parsec_tiled_matrix_t *)dcIPIV,
                           (parsec_tiled_matrix_t *)&dcInvA );

            /* Check the solution */
            check_inverse(parsec_ctx, (rank_A == 0) ? loud : 0,
                                 (parsec_tiled_matrix_t *)&dcA0,
                                 (parsec_tiled_matrix_t *)&dcInvA,
                                 (parsec_tiled_matrix_t *)&dcI);
        }
    }

    if ( check ) {
        parsec_data_free(dcA0.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA0);
        parsec_data_free(dcA_out.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcA_out);
        parsec_data_free(dcB.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcB);
        parsec_data_free(dcX.mat);
        parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcX);
        if ( check_inv ) {
            parsec_data_free(dcInvA.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcInvA);
            parsec_data_free(dcI.mat);
            parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&dcI);
        }
    }

    if ( *M != *N ) {
        fprintf(stderr, "Check is impossible if M != N\n");
    }
#endif

    dcA = redistribute_lapack_output_cleanup(&dcA_lapack, dcA, 1, comm_A, rank_A, "redisA");

    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", dcA);

    /* SCALAPACK IPIV:
     * looking at checking on pdegtrf, the IPIV is passed to
     *   pdgetrs_( "N", &n, &s, Alu, &i1, &i1, descA, ippiv, X, &i1, &i1, descB, &info );
     * which invokes
     *   CALL pdlapiv( 'Forward', 'Row', 'Col', n, nrhs, b, ib, jb,
     *  $                descb, ipiv, ia, 1, descip, idum1 )
     * on pdlapiv documentation:
     *  DIREC   (global input) CHARACTER*1
     *          Specifies in which order the permutation is applied:
     *            = 'F' (Forward) Applies pivots Forward from top of matrix.
     *                  Computes P*sub( A ).
     *            = 'B' (Backward) Applies pivots Backward from bottom of
     *                  matrix. Computes inv( P )*sub( A ).
     *
     *  ROWCOL  (global input) CHARACTER*1
     *          Specifies if the rows or columns are to be permuted:
     *             = 'R' Rows will be permuted,
     *             = 'C' Columns will be permuted.
     *
     *  PIVROC  (global input) CHARACTER*1
     *          Specifies whether IPIV is distributed over a process row
     *          or column:
     *          = 'R' IPIV distributed over a process row
     *          = 'C' IPIV distributed over a process column
     *
     *  IPIV    (local input) INTEGER array, dimension (LIPIV) where LIPIV is
     *          when ROWCOL='R' or 'r':
     *             >= LOCr( IA+M-1 ) + MB_A      if PIVROC='C' or 'c',
     *             >= LOCc( M + MOD(JP-1,NB_P) ) if PIVROC='R' or 'r', and,
     *          when ROWCOL='C' or 'c':
     *             >= LOCr( N + MOD(IP-1,MB_P) ) if PIVROC='C' or 'c',
     *             >= LOCc( JA+N-1 ) + NB_A      if PIVROC='R' or 'r'.
     *          This array contains the pivoting information. IPIV(i) is the
     *          global row (column), local row (column) i was swapped with.
     *          When ROWCOL='R' or 'r' and PIVROC='C' or 'c', or ROWCOL='C'
     *          or 'c' and PIVROC='R' or 'r', the last piece of this array of
     *          size MB_A (resp. NB_A) is used as workspace. In those cases,
     *          this array is tied to the distributed matrix A.
     */

    /* Translate to global rows naming expected by scalapack */
    for(int n = 0; n < dcIPIV->super.lln; n++){
        int local_tile_row = (n / dcIPIV->super.nb);
        int global_tile_row = local_tile_row * (dcIPIV->grid.cols) + dcIPIV->grid.crank;
        int global_row_offset = global_tile_row * dcIPIV->super.nb;
        int global_row = global_row_offset + n % dcIPIV->super.nb;
        if(((int*)dcIPIV->mat)[ n ] == 0) { /* no swap */
            ((int*)dcIPIV->mat)[ n ] = global_row;
        }else{
            ((int*)dcIPIV->mat)[ n ] = ((int*)dcIPIV->mat)[ n ] + global_row_offset;
        }
    }

    /* For the output ScaLAPACK IPIV every process has the IPIV for its local rows.
     * therefore, it's not a vector representing row or a column, it's a matrix.
     * distributed 2dbc.
     * On dplasma it's a vector 1xM distributed among 1xP*Q processes.
     * On ScaLAPACK it's a matrix distributed amon PxQ processes where each
     * process has a vector 1xmloc.
     * Therefore, we can't redistribute dplasma IPIV to ScaLAPACK IPIV.
     * We reconstructed gathering first a rank 0.
     */

    /* Reconstruct IPIV local rows as requiered by the spec.
     * We ran 1D 1xP*Q getrf_1D.
     * DPLASMA IPIV is spread across the 1xQ*P column processes.
     * --> Gather global IPIV on rank 0.
     */
    int *recv_count = (int*) malloc(sizeof(int)*nodes_A);
    int *displs = (int*) malloc(sizeof(int)*nodes_A);
    recv_count[rank_A] = dcIPIV->super.nb_local_tiles;
    MPI_Allgather(&(recv_count[rank_A]), 1, MPI_INT, recv_count, 1, MPI_INT, comm_A);

    MPI_Datatype dt_ipiv, tmp_dt_ipiv;
    MPI_Type_vector(1,
                    redisNB,
                    redisQ*redisNB,
                    MPI_INT,
                    &tmp_dt_ipiv);
    MPI_Aint ub, lb, extent;
    MPI_Type_get_extent(tmp_dt_ipiv, &lb, &extent);
    ub = (lb + extent)*(redisQ);
    MPI_Type_create_resized(tmp_dt_ipiv, 0, ub, &dt_ipiv);
    MPI_Type_commit(&dt_ipiv);


    int r;
    displs[0] = 0;
    for(r = 1; r < nodes_A; r++){
        displs[r] = displs[r-1] + redisNB;
    }

    /* send/recv redis tiles info */
    int ipiv_size_redis = (rank_A == 0)?
                          dcIPIV->super.n * dcIPIV->super.nb: /* gather all IPIV: dcIPIV->super.nb=redisNB & TILED matrix: sending FULL tile */
                          gN_IPIV; /* rest of ranks only care about the relevant part of the global IPIV*/
    int *tmp_ipiv = (int*)malloc(sizeof(int)*(ipiv_size_redis));
    if(rank_A == 0){
        MPI_Request req;
        MPI_Isend(dcIPIV->mat, recv_count[rank_A]*redisNB, MPI_INT, 0,
                 rank_A, comm_A, &req);
        for(r = 0; r < nodes_A; r++){
            MPI_Recv(&(tmp_ipiv[r*redisNB]), recv_count[r], dt_ipiv, r,
                     r, comm_A, MPI_STATUS_IGNORE);
        }
        MPI_Wait(&req, MPI_STATUS_IGNORE);
    }else{
        MPI_Send(dcIPIV->mat, recv_count[rank_A]*redisNB, MPI_INT, 0,
                 rank_A, comm_A);
    }
    MPI_Type_free(&tmp_dt_ipiv);
    MPI_Type_free(&dt_ipiv);

    /* Broadcast only the relevant part of the global IPIV */
    MPI_Bcast(tmp_ipiv, gN_IPIV, MPI_INT, 0, comm_A);

    for(int m = 0; m < mloc_A; m++){ /* for each local row of the original lapack matrix*/
        int local_tile_row = (m / dcA_lapack.super.mb);
        int global_tile_row = local_tile_row * (dcA_lapack.grid.rows) + dcA_lapack.grid.rrank;
        int global_row = global_tile_row * dcA_lapack.super.mb + m % dcA_lapack.super.mb;
        IPIV[m] = tmp_ipiv[global_row];
    }

    free(tmp_ipiv);

    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)dcA);
    parsec_data_free(dcIPIV->mat);
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)dcIPIV);
}

 /*-------------------------------------------------------------------*/

#ifdef CHECK_RESULTS
static int check_solution( parsec_context_t *parsec, int loud,
                            parsec_tiled_matrix_t *dcA,
                            parsec_tiled_matrix_t *dcB,
                            parsec_tiled_matrix_t *dcX )
 {
     int info_solution;
     double Rnorm = 0.0;
     double Anorm = 0.0;
     double Bnorm = 0.0;
     double Xnorm, result;
     int m = dcB->m;
     double eps = LAPACKE_dlamch_work('e');

     Anorm = dplasma_dlange(parsec, PlasmaInfNorm, dcA);
     Bnorm = dplasma_dlange(parsec, PlasmaInfNorm, dcB);
     Xnorm = dplasma_dlange(parsec, PlasmaInfNorm, dcX);

     /* Compute b - A*x */
     dplasma_dgemm( parsec, PlasmaNoTrans, PlasmaNoTrans, -1.0, dcA, dcX, 1.0, dcB);

     Rnorm = dplasma_dlange(parsec, PlasmaInfNorm, dcB);

     result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * m * eps ) ;

     if ( loud > 2 ) {
         printf("============\n");
         printf("Checking the Residual of the solution \n");
         if ( loud > 3 )
             printf( "-- ||A||_oo = %e, ||X||_oo = %e, ||B||_oo= %e, ||A X - B||_oo = %e\n",
                     Anorm, Xnorm, Bnorm, Rnorm );
         printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);
     }

     if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
         if( loud ) printf("-- Solution is suspicious ! \n");
         info_solution = 1;
     }
     else{
         if( loud ) printf("-- Solution is CORRECT (RES)! \n");
         info_solution = 0;
     }

     return info_solution;
 }

 static int check_inverse( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_t *dcA,
                           parsec_tiled_matrix_t *dcInvA,
                           parsec_tiled_matrix_t *dcI )
 {
     int info_solution;
     double Anorm    = 0.0;
     double InvAnorm = 0.0;
     double Rnorm, result;
     int m = dcA->m;
     double eps = LAPACKE_dlamch_work('e');

     Anorm    = dplasma_dlange(parsec, PlasmaInfNorm, dcA   );
     InvAnorm = dplasma_dlange(parsec, PlasmaInfNorm, dcInvA);

     /* Compute I - A*A^{-1} */
     dplasma_dgemm( parsec, PlasmaNoTrans, PlasmaNoTrans, -1.0, dcA, dcInvA, 1.0, dcI);

     Rnorm = dplasma_dlange(parsec, PlasmaInfNorm, dcI);

     result = Rnorm / ( ( Anorm * InvAnorm ) * m * eps ) ;

     if ( loud > 2 ) {
         printf("============\n");
         printf("Checking the Residual of the solution \n");
         if ( loud > 3 )
             printf( "-- ||A||_oo = %e, ||A^{-1}||_oo = %e, ||A A^{-1} - I||_oo = %e\n",
                     Anorm, InvAnorm, Rnorm );
         printf("-- ||AA^{-1}-I||_oo/((||A||_oo||A^{-1}||_oo).N.eps) = %e \n", result);
     }

     if (  isnan(Rnorm) || isinf(Rnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
         if( loud ) printf("-- Solution is suspicious ! \n");
         info_solution = 1;
     }
     else{
         if( loud ) printf("-- Solution is CORRECT (INV)! \n");
         info_solution = 0;
     }

     return info_solution;
 }
#endif

GENERATE_F77_BINDINGS (PDGETRF,
                       pdgetrf,
                       pdgetrf_,
                       pdgetrf__,
                       pdgetrf_w,
                       (int * M, int * N, double * A, int * IA, int * JA, int * DESCA, int * IPIV, int * info),
                       (M, N, A, IA, JA, DESCA, IPIV, info))
