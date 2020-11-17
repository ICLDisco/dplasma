#include "common.h"

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

static int check_solution( parsec_context_t *parsec, int loud,
                          parsec_tiled_matrix_dc_t *dcA,
                          parsec_tiled_matrix_dc_t *dcB,
                          parsec_tiled_matrix_dc_t *dcX );

static int check_inverse( parsec_context_t *parsec, int loud,
                         parsec_tiled_matrix_dc_t *dcA,
                         parsec_tiled_matrix_dc_t *dcInvA,
                         parsec_tiled_matrix_dc_t *dcI );

void pdgetrf_w(int * M,
              int * N,
              double * A,
              int * IA,
              int * JA,
              int * DESCA,
              int * IPIV,
              int * info){


#ifdef COUNT_WRAPPED_CALLS
    count_PDGETRF_NOPIV++;
#endif
    *info=0;
    if( (*M == 0) || (*N == 0)){
      /* NOP */
      return;
    }

    /* Only for diagonal dominant matrix. IPIV is faked at the end. */

    int KP = 1;
    int KQ = 1;

    PASTE_SETUP(A);

    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  "M%d N%d IA%d JA%d (ictxt)DESCA[WRAPPER_CTXT1_] %d, "
          "(gM)DESCA[WRAPPER_M1_] %d, (gN)DESCA[WRAPPER_N1_] %d, (MB)DESCA[WRAPPER_MB1_] %d, (NB)DESCA[WRAPPER_NB1_] %d, "
          "DESCA[WRAPPER_RSRC1_] %d, DESCA[WRAPPER_CSRC1_] %d, (LDD)DESCA[WRAPPER_LLD1_] %d\n",
          *M, *N, *IA, *JA, DESCA[WRAPPER_CTXT1_],
          DESCA[WRAPPER_M1_], DESCA[WRAPPER_N1_], DESCA[WRAPPER_MB1_], DESCA[WRAPPER_NB1_],
          DESCA[WRAPPER_RSRC1_], DESCA[WRAPPER_CSRC1_], DESCA[WRAPPER_LLD1_]);

    parsec_init_wrapped_call((void*)comm_A);

    PASTE_CODE_INIT_LAPACK_MATRIX(dcA, two_dim_block_cyclic, A,
                                  (&dcA, matrix_RealDouble, matrix_Lapack,
                                   nodes_A, rank_A,
                                   MB_A, NB_A,
                                   gM_A, gN_A,
                                   cIA, cJA,
                                   *M, *N,
                                   KP, KQ,
                                   iP_A, jQ_A,
                                   P_A,
                                   nloc_A, LDD_A));

#ifdef CHECK_RESULTS
    int check = 1;
    int check_inv = 1;

    if ( *M != *N && check ) {
        fprintf(stderr, "Check is impossible if M != N\n");
        check = 0;
    }

    int cLDA = *M;// max(M, iparam[IPARAM_LDA])
    //cLDA = max(*M, LDA); no in case submatrix

    int NRHS  = KP;
    //iparam[IPARAM_LDB] = -'m';
    int LDB   = -'m';//max(N, iparam[IPARAM_LDB]);
    LDB = *M;//max(*M, LDB);no in case submatrix

    PASTE_CODE_ALLOCATE_MATRIX(dcA_out, check,
                               two_dim_block_cyclic, (&dcA_out, matrix_RealDouble, matrix_Tile,
                                                      nodes_A, rank_A, MB_A, NB_A, cLDA, *N, 0, 0,
                                                      *M, *N, KP, KQ, 0, 0, P_A));
    PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
                               two_dim_block_cyclic, (&dcA0, matrix_RealDouble, matrix_Tile,
                                                      nodes_A, rank_A, MB_A, NB_A, cLDA, *N, 0, 0,
                                                      *M, *N, KP, KQ, 0, 0, P_A));

    /* Random B check */
    PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
                               two_dim_block_cyclic, (&dcB, matrix_RealDouble, matrix_Tile,
                                                      nodes_A, rank_A, MB_A, NB_A, LDB, NRHS, 0, 0,
                                                      *M, NRHS, KP, KQ, 0, 0, P_A));
    PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
                               two_dim_block_cyclic, (&dcX, matrix_RealDouble, matrix_Tile,
                                                      nodes_A, rank_A, MB_A, NB_A, LDB, NRHS, 0, 0,
                                                      *M, NRHS, KP, KQ, 0, 0, P_A));
    /* Inverse check */
    PASTE_CODE_ALLOCATE_MATRIX(dcInvA, check_inv,
                               two_dim_block_cyclic, (&dcInvA, matrix_RealDouble, matrix_Tile,
                                                      nodes_A, rank_A, MB_A, NB_A, cLDA, *N, 0, 0,
                                                      *M, *N, KP, KQ, 0, 0, P_A));
    PASTE_CODE_ALLOCATE_MATRIX(dcI, check_inv,
                               two_dim_block_cyclic, (&dcI, matrix_RealDouble, matrix_Tile,
                                                      nodes_A, rank_A, MB_A, NB_A, cLDA, *N, 0, 0,
                                                      *M, *N, KP, KQ, 0, 0, P_A));


    if ( check ) {
        dplasma_dlacpy( parsec_ctx, PlasmaUpperLower,
                        (parsec_tiled_matrix_dc_t *)&dcA,
                        (parsec_tiled_matrix_dc_t *)&dcA0 );
        dplasma_dplrnt( parsec_ctx, 0, (parsec_tiled_matrix_dc_t *)&dcB, 2354);
        dplasma_dlacpy( parsec_ctx, PlasmaUpperLower,
                        (parsec_tiled_matrix_dc_t *)&dcB,
                        (parsec_tiled_matrix_dc_t *)&dcX );
    }
    if ( check_inv ) {
        dplasma_dlaset( parsec_ctx, PlasmaUpperLower, 0., 1., (parsec_tiled_matrix_dc_t *)&dcI);
        dplasma_dlaset( parsec_ctx, PlasmaUpperLower, 0., 1., (parsec_tiled_matrix_dc_t *)&dcInvA);
    }

#endif


#ifdef WRAPPER_VERBOSE
    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", dcA, P_A, Q_A);
#endif

#ifdef MEASURE_INTERNAL_TIMES
    PASTE_CODE_FLOPS(FLOPS_DGETRF, ((DagDouble_t)*M,(DagDouble_t)*N));
#endif
    WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec_ctx, dgetrf_nopiv,
                              ((parsec_tiled_matrix_dc_t*)&dcA, info),
                              dplasma_dgetrf_nopiv_Destruct( PARSEC_dgetrf_nopiv ),
                              rank_A, P_A, Q_A, NB_A, gN_A, comm_A);


#ifdef CHECK_RESULTS
    int loud=5;
    if ( check ) {
        dcopy_lapack_tile(parsec_ctx, &dcA, &dcA_out, mloc_A, nloc_A);

        dplasma_dtrsm( parsec_ctx, PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit,    1.0,
                       (parsec_tiled_matrix_dc_t*)&dcA_out,
                       (parsec_tiled_matrix_dc_t*)&dcX);
        dplasma_dtrsm( parsec_ctx, PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit, 1.0,
                       (parsec_tiled_matrix_dc_t*)&dcA_out,
                       (parsec_tiled_matrix_dc_t*)&dcX);

        /* Check the solution */
        check_solution( parsec_ctx, (rank_A == 0) ? loud : 0,
                               (parsec_tiled_matrix_dc_t *)&dcA0,
                               (parsec_tiled_matrix_dc_t *)&dcB,
                               (parsec_tiled_matrix_dc_t *)&dcX);

        /*
         * Second check with inverse
         */
        if ( check_inv ) {
            dplasma_dtrsm( parsec_ctx, PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit,    1.0,
                           (parsec_tiled_matrix_dc_t*)&dcA_out,
                           (parsec_tiled_matrix_dc_t*)&dcInvA);
            dplasma_dtrsm( parsec_ctx, PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit, 1.0,
                           (parsec_tiled_matrix_dc_t*)&dcA_out,
                           (parsec_tiled_matrix_dc_t*)&dcInvA);

            /* Check the solution */
            check_inverse(parsec_ctx, (rank_A == 0) ? loud : 0,
                                 (parsec_tiled_matrix_dc_t *)&dcA0,
                                 (parsec_tiled_matrix_dc_t *)&dcInvA,
                                 (parsec_tiled_matrix_dc_t *)&dcI);
        }
    }

    if ( check ) {
        parsec_data_free(dcA_out.mat);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA_out);

        parsec_data_free(dcA0.mat);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA0);
        parsec_data_free(dcB.mat);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcB);
        parsec_data_free(dcX.mat);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcX);
        if ( check_inv ) {
            parsec_data_free(dcInvA.mat);
            parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcInvA);
            parsec_data_free(dcI.mat);
            parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcI);
        }
    }
    if ( *M != *N ) {
        fprintf(stderr, "Check is impossible if M != N\n");
    }

#endif

    /* Faking IPIV */
    int lrw=1;
    int cblock = myrow_A;
    int i_ipiv;
    for(i_ipiv = 0; i_ipiv<mloc_A; i_ipiv++){
        ((int*)IPIV)[i_ipiv] = lrw + cblock*MB_A;
        lrw++;
        if(lrw>MB_A){
            lrw=1;
            cblock+=P_A;
        }
    }

#ifdef WRAPPER_VERBOSE
    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", dcA, P_A, Q_A);
    /*
     *  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
     *          This array contains the pivoting information.
     *          IPIV(i) -> The global row local row i was swapped with.
     *          This array is tied to the distributed matrix A.
     *  However, in the original ScaLAPACK GETRF every process has the full
     *  IPIV containing global row indexes, therefore, we  build it.
     */
    int gN_IPIV = dplasma_imin(*M, *N);
    int nloc_IPIV = nloc_A;//mloc + MB;//LETS TRY //+ MB;
    int LD_IPIV = 1;/*assuming this is not a submatrix and that it is a row*/
    two_dim_block_cyclic_t  dcIPIV;
    two_dim_block_cyclic_lapack_init(&dcIPIV, matrix_Integer, matrix_Lapack,
                               nodes_A, rank_A,
                               1, NB_A,
                               1, gN_IPIV,
                               cIA, cJA,
                               1, gN_IPIV,
                               KP, KQ,
                               iP_A, jQ_A,
                               P_A,
                               nloc_IPIV, LD_IPIV);
    dcIPIV.mat = (void*)IPIV;
    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcIPIV", dcIPIV, P_A, Q_A);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcIPIV);
#endif

    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA);
}

 /*-------------------------------------------------------------------*/

static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_dc_t *dcA,
                           parsec_tiled_matrix_dc_t *dcB,
                           parsec_tiled_matrix_dc_t *dcX )
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
        if( loud ) printf("-- Solution is CORRECT ! \n");
        info_solution = 0;
    }

    return info_solution;
}

static int check_inverse( parsec_context_t *parsec, int loud,
                          parsec_tiled_matrix_dc_t *dcA,
                          parsec_tiled_matrix_dc_t *dcInvA,
                          parsec_tiled_matrix_dc_t *dcI )
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
        if( loud ) printf("-- Solution is CORRECT ! \n");
        info_solution = 0;
    }

    return info_solution;
}

GENERATE_F77_BINDINGS (PDGETRF,
                       pdgetrf,
                       pdgetrf_,
                       pdgetrf__,
                       pdgetrf_w,
                       (int * M, int * N, double * A, int * IA, int * JA, int * DESCA, int * IPIV, int * info),
                       (M, N, A, IA, JA, DESCA, IPIV, info))
