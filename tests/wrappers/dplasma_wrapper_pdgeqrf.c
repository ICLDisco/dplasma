#include "common.h"
  /************************************************************************
    Version not compliant with scalapack.
    LAPACK 3.7.0 introduces the LATSQR routine implementing the
    sequential TSQR (Tall Skinny QR) factorization algorithm which
    corresponds to the DPLASMA implementation of the QR factorization in the GEQRF
    Wrapper would need to use CORHR reconstruct the results expected by scalapack (not release).
  ************************************************************************/

/*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
*  INTEGER            IA, INFO, JA, LWORK, M, N
*     ..
*     .. Array Arguments ..
*  INTEGER            DESCA( * )
*  DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGEQRF computes a QR factorization of a real distributed M-by-N
*  matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) = Q * R.
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
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW *MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL *NB_A
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
*          On entry, the local pieces of the M-by-N distributed matrix
*          sub( A ) which is to be factored.  On exit, the elements on
*          and above the diagonal of sub( A ) contain the min(M,N) by N
*          upper trapezoidal matrix R (R is upper triangular if M >= N);
*          the elements below the diagonal, with the array TAU,
*          represent the orthogonal matrix Q as a product of elementary
*          reflectors (see Further Details).
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
*  TAU     (local output) DOUBLE PRECISION array, dimension
*          LOCc(JA+MIN(M,N)-1). This array contains the scalar factors
*          TAU of the elementary reflectors. TAU is tied to the
*          distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= NB_A * ( Mp0 + Nq0 + NB_A ), where
*
*          IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ),
*          Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ),
*
*          and NUMROC, INDXG2P are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = *100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(ja) H(ja+1) . . . H(ja+k-1), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(j) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with v(1:i-1) = 0
*  and v(i) = 1; v(i+1:m) is stored on exit in A(ia+i:ia+m-1,ja+i-1),
*  and tau in TAU(ja+i-1).
*
*  =====================================================================
*
*/

#ifdef CHECK_RESULTS
static int check_orthogonality(parsec_context_t *parsec, int loud,
                               parsec_tiled_matrix_dc_t *Q);
static int check_factorization(parsec_context_t *parsec, int loud,
                               parsec_tiled_matrix_dc_t *Aorig,
                               parsec_tiled_matrix_dc_t *A,
                               parsec_tiled_matrix_dc_t *Q);
static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_dc_t *dcA,
                           parsec_tiled_matrix_dc_t *dcB,
                           parsec_tiled_matrix_dc_t *dcX );
#endif

void pdgeqrf_w(int * M,
              int * N,
              double * A,
              int * IA,
              int * JA,
              int * DESCA,
              double * TAU,
              double * WORK,
              int * LWORK,
              int * info){


#ifdef COUNT_WRAPPED_CALLS
    count_PDGEQRF++;
#endif
    *info=0;
    if( (*M == 0) || (*N == 0)){
      /* NOP */
      return;
    }

    int KP = 1;
    int KQ = 1;

    PASTE_SETUP(A);

#ifdef WRAPPER_VERBOSE_CALLS
    if(rank_A == 0){
      printf("V-PDGEQRF M%d N%d "
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

    //TODO Not doing all the check done in scalapack
    WORK[0]= (double)(NB_A * ( mloc_A + nloc_A + NB_A ));
    if(*LWORK==-1){
        *info = 0;
        return;
    }

    parsec_init_wrapped_call((void*)comm_A);

    two_dim_block_cyclic_t dcA_lapack;
    two_dim_block_cyclic_lapack_init(&dcA_lapack, matrix_RealDouble, matrix_Lapack,
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

    int IB = 1;
    char *var_IB= getenv("PARSEC_WRAPPER_IB");
    if(var_IB!=NULL){
        IB = atoi(var_IB);
    }

 /* @param[out] T
  *          Descriptor of the matrix T distributed exactly as the A matrix. T.mb
  *          defines the IB parameter of tile QR algorithm. This matrix must be
  *          of size A.mt * T.mb - by - A.nt * T.nb, with T.nb == A.nb.
  *          On exit, contains auxiliary information required to compute the Q
  *          matrix, and/or solve the problem.
  */
    int MB_T = IB;
    int NB_T = NB_A;
    int MT = dcA_lapack.super.mt;
    int M_T = MB_T * MT;
    int N_T = NB_T * dcA_lapack.super.nt;
    int nloc_T = NB_T * dcA_lapack.super.lnt;
    int LLD_T = MB_T * dcA_lapack.super.lmt;

    two_dim_block_cyclic_t dcT_lapack;
    two_dim_block_cyclic_lapack_init(&dcT_lapack, matrix_RealDouble, matrix_Lapack,
                                     rank_A,
                                     MB_T, NB_T,
                                     M_T, N_T,
                                     cIA, cJA,
                                     M_T, N_T,
                                     P_A, Q_A,
                                     KP, KQ,
                                     iP_A, jQ_A,
                                     LLD_T, nloc_T);
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcT_lapack, "dcT_lapack");

    dcT_lapack.mat = malloc( (size_t)LLD_T * (size_t)nloc_T *
            (size_t)parsec_datadist_getsizeoftype(dcT_lapack.super.mtype));

    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", (&dcA_lapack));
    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcT", (&dcT_lapack));

#ifdef CHECK_RESULTS
    dplasma_dlaset( parsec_ctx, PlasmaUpperLower,
            0., 0., (parsec_tiled_matrix_dc_t *)&dcT_lapack);
#endif

#ifdef CHECK_RESULTS
    int check = 1;
    int cLDA = *M;
    int NRHS  = KP;
    int LDB   = *M;//max(*M, LDB);no in case submatrix
    PASTE_CODE_ALLOCATE_MATRIX(dcA0, check,
        two_dim_block_cyclic, (&dcA0, matrix_RealDouble, matrix_Tile,
                               rank_A, MB_A, NB_A, cLDA, *N, 0, 0,
                               *M, *N, P_A, Q_A, KP, KQ, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcA_out, check,
        two_dim_block_cyclic, (&dcA_out, matrix_RealDouble, matrix_Tile,
                               rank_A, MB_A, NB_A, cLDA, *N, 0, 0,
                               *M, *N, P_A, Q_A, KP, KQ, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcT_out, check,
        two_dim_block_cyclic, (&dcT_out, matrix_RealDouble, matrix_Tile,
                               rank_A, MB_T, NB_T, M_T, N_T, 0, 0,
                               M_T, N_T, P_A, Q_A, KP, KQ, 0, 0));

    PASTE_CODE_ALLOCATE_MATRIX(dcQ, check,
        two_dim_block_cyclic, (&dcQ, matrix_RealDouble, matrix_Tile,
                               rank_A, MB_A, NB_A, cLDA, *N, 0, 0,
                               *M, *N, P_A, Q_A, KP, KQ, 0, 0));

    /* Check the solution */
    PASTE_CODE_ALLOCATE_MATRIX(dcB, check,
        two_dim_block_cyclic, (&dcB, matrix_RealDouble, matrix_Tile,
                               rank_A, MB_A, NB_A, LDB, NRHS, 0, 0,
                               *M, NRHS, P_A, Q_A, KP, KQ, 0, 0));
    PASTE_CODE_ALLOCATE_MATRIX(dcX, check,
        two_dim_block_cyclic, (&dcX, matrix_RealDouble, matrix_Tile,
                               rank_A, MB_A, NB_A, LDB, NRHS, 0, 0,
                               *M, NRHS, P_A, Q_A, KP, KQ, 0, 0));

    if( check ){
        dcopy_lapack_tile(parsec_ctx, &dcA_lapack, &dcA0, mloc_A, nloc_A);
    }
#endif

    int redisA = 0, redisT = 0;
    two_dim_block_cyclic_t *dcA = redistribute_lapack_input(&dcA_lapack, redisA, comm_A, rank_A, "redisA");
    two_dim_block_cyclic_t *dcT = redistribute_lapack_input(&dcT_lapack, redisT, comm_A, rank_A, "redisT");

#ifdef MEASURE_INTERNAL_TIMES
    PASTE_CODE_FLOPS(FLOPS_DGEQRF, ((DagDouble_t)*M, (DagDouble_t)*N));
#endif

    WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(parsec_ctx, dgeqrf,
                              ((parsec_tiled_matrix_dc_t*)dcA,
                               (parsec_tiled_matrix_dc_t*)dcT),
                              dplasma_dgeqrf_Destruct( PARSEC_dgeqrf ),
                              rank_A, P_A, Q_A, NB_A, gN_A, comm_A);

    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcA", dcA);
    PRINT(parsec_ctx, comm_A, PlasmaUpperLower, "dcT", dcT);

#ifdef CHECK_RESULTS
    if( check ) {
        dcopy_lapack_tile(parsec_ctx, dcA, &dcA_out, mloc_A, nloc_A);
        dcopy_lapack_tile(parsec_ctx, dcT, &dcT_out, dcT->super.llm, dcT->super.lln);
        int loud=5;
        int ret;
        if (*M >= *N) {
            if(loud > 2) printf("+++ Generate the Q ...");
            dplasma_dorgqr( parsec_ctx,
                            (parsec_tiled_matrix_dc_t *)&dcA_out,
                            (parsec_tiled_matrix_dc_t *)&dcT_out,
                            (parsec_tiled_matrix_dc_t *)&dcQ);

            if(loud > 2) printf("Done\n");

            if(loud > 2) printf("+++ Solve the system ...");
            dplasma_dplrnt( parsec_ctx, 0, (parsec_tiled_matrix_dc_t *)&dcX, 2354);
            dplasma_dlacpy( parsec_ctx, PlasmaUpperLower,
                            (parsec_tiled_matrix_dc_t *)&dcX,
                            (parsec_tiled_matrix_dc_t *)&dcB );
            dplasma_dgeqrs( parsec_ctx,
                            (parsec_tiled_matrix_dc_t *)&dcA_out,
                            (parsec_tiled_matrix_dc_t *)&dcT_out,
                            (parsec_tiled_matrix_dc_t *)&dcX );
            if(loud > 2) printf("Done\n");

            /* Check the orthogonality, factorization and the solution */
            ret |= check_orthogonality( parsec_ctx, (rank_A == 0) ? loud : 0,
                                        (parsec_tiled_matrix_dc_t *)&dcQ);
            ret |= check_factorization( parsec_ctx, (rank_A == 0) ? loud : 0,
                                        (parsec_tiled_matrix_dc_t *)&dcA0,
                                        (parsec_tiled_matrix_dc_t *)&dcA_out,
                                        (parsec_tiled_matrix_dc_t *)&dcQ );
            ret |= check_solution( parsec_ctx, (rank_A == 0) ? loud : 0,
                                   (parsec_tiled_matrix_dc_t *)&dcA0,
                                   (parsec_tiled_matrix_dc_t *)&dcB,
                                   (parsec_tiled_matrix_dc_t *)&dcX );

        } else {
            printf("Check cannot be performed when N > M\n");
        }

        parsec_data_free(dcA0.mat);
        parsec_data_free(dcA_out.mat);
        parsec_data_free(dcT_out.mat);
        parsec_data_free(dcQ.mat);
        parsec_data_free(dcB.mat);
        parsec_data_free(dcX.mat);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA0);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA_out);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcT_out);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcQ);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcB);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcX);
    }

#endif

    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)dcA);
    parsec_data_free(dcT->mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)dcT);
}


#ifdef CHECK_RESULTS
/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int check_orthogonality(parsec_context_t *parsec, int loud, parsec_tiled_matrix_dc_t *Q)
{
    two_dim_block_cyclic_t *twodQ = (two_dim_block_cyclic_t *)Q;
    double normQ = 999999.0;
    double result;
    double eps = LAPACKE_dlamch_work('e');
    int info_ortho;
    int M = Q->m;
    int N = Q->n;
    int minMN = min(M, N);

    PASTE_CODE_ALLOCATE_MATRIX(Id, 1,
        two_dim_block_cyclic, (&Id, matrix_RealDouble, matrix_Tile,
                               twodQ->grid.rank,
                               Q->mb, Q->nb, minMN, minMN, 0, 0,
                               minMN, minMN,
                               twodQ->grid.rows, twodQ->grid.cols, twodQ->grid.krows, twodQ->grid.kcols, twodQ->grid.ip, twodQ->grid.jq));

    dplasma_dlaset( parsec, PlasmaUpperLower, 0., 1., (parsec_tiled_matrix_dc_t *)&Id);

    /* Perform Id - Q'Q */
    if ( M >= N ) {
        dplasma_dsyrk( parsec, PlasmaUpper, PlasmaTrans,
                       1.0, Q, -1.0, (parsec_tiled_matrix_dc_t*)&Id );
    } else {
        dplasma_dsyrk( parsec, PlasmaUpper, PlasmaNoTrans,
                       1.0, Q, -1.0, (parsec_tiled_matrix_dc_t*)&Id );
    }

    normQ = dplasma_dlansy(parsec, PlasmaInfNorm, PlasmaUpper, (parsec_tiled_matrix_dc_t*)&Id);

    result = normQ / (minMN * eps);
    if ( loud ) {
        printf("============\n");
        printf("Checking the orthogonality of Q \n");
        printf("||Id-Q'*Q||_oo / (N*eps) = %e \n", result);
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        if ( loud ) printf("-- Orthogonality is suspicious ! \n");
        info_ortho=1;
    }
    else {
        if ( loud ) printf("-- Orthogonality is CORRECT ! \n");
        info_ortho=0;
    }

    parsec_data_free(Id.mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&Id);
    return info_ortho;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int
check_factorization(parsec_context_t *parsec, int loud,
                    parsec_tiled_matrix_dc_t *Aorig,
                    parsec_tiled_matrix_dc_t *A,
                    parsec_tiled_matrix_dc_t *Q)
{
    parsec_tiled_matrix_dc_t *subA;
    two_dim_block_cyclic_t *twodA = (two_dim_block_cyclic_t *)A;
    double Anorm, Rnorm;
    double result;
    double eps = LAPACKE_dlamch_work('e');
    int info_factorization;
    int M = A->m;
    int N = A->n;
    int minMN = min(M, N);

    PASTE_CODE_ALLOCATE_MATRIX(Residual, 1,
        two_dim_block_cyclic, (&Residual, matrix_RealDouble, matrix_Tile,
                               twodA->grid.rank,
                               A->mb, A->nb, M, N, 0, 0,
                               M, N,
                               twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols,
                               twodA->grid.ip, twodA->grid.jq));

    PASTE_CODE_ALLOCATE_MATRIX(R, 1,
        two_dim_block_cyclic, (&R, matrix_RealDouble, matrix_Tile,
                               twodA->grid.rank,
                               A->mb, A->nb, N, N, 0, 0,
                               N, N,
                               twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols,
                               twodA->grid.ip, twodA->grid.jq));

    /* Copy the original A in Residual */
    dplasma_dlacpy( parsec, PlasmaUpperLower, Aorig, (parsec_tiled_matrix_dc_t *)&Residual );

    /* Extract the R */
    dplasma_dlaset( parsec, PlasmaUpperLower, 0., 0., (parsec_tiled_matrix_dc_t *)&R);

    subA = tiled_matrix_submatrix( A, 0, 0, N, N );
    dplasma_dlacpy( parsec, PlasmaUpper, subA, (parsec_tiled_matrix_dc_t *)&R );
    free(subA);

    /* Perform Residual = Aorig - Q*R */
    dplasma_dgemm( parsec, PlasmaNoTrans, PlasmaNoTrans,
                   -1.0, Q, (parsec_tiled_matrix_dc_t *)&R,
                    1.0, (parsec_tiled_matrix_dc_t *)&Residual);

    /* Free R */
    parsec_data_free(R.mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&R);

    Rnorm = dplasma_dlange(parsec, PlasmaInfNorm, (parsec_tiled_matrix_dc_t*)&Residual);
    Anorm = dplasma_dlange(parsec, PlasmaInfNorm, Aorig);

    result = Rnorm / ( Anorm * minMN * eps);

    if ( loud ) {
        printf("============\n");
        printf("Checking the QR Factorization \n");
        printf("-- ||A-QR||_oo/(||A||_oo.N.eps) = %e \n", result );
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        if ( loud ) printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else {
        if ( loud ) printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    parsec_data_free(Residual.mat);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&Residual);
    return info_factorization;
}

static int check_solution( parsec_context_t *parsec, int loud,
                           parsec_tiled_matrix_dc_t *dcA,
                           parsec_tiled_matrix_dc_t *dcB,
                           parsec_tiled_matrix_dc_t *dcX )
{
    parsec_tiled_matrix_dc_t *subX;
    int info_solution;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double Bnorm = 0.0;
    double Xnorm, result;
    double eps = LAPACKE_dlamch_work('e');

    subX = tiled_matrix_submatrix( dcX, 0, 0, dcA->n, dcX->n );

    Anorm = dplasma_dlange(parsec, PlasmaInfNorm, dcA);
    Bnorm = dplasma_dlange(parsec, PlasmaInfNorm, dcB);
    Xnorm = dplasma_dlange(parsec, PlasmaInfNorm, subX);

    /* Compute A*x-b */
    dplasma_dgemm( parsec, PlasmaNoTrans, PlasmaNoTrans, 1.0, dcA, subX, -1.0, dcB);

    /* Compute A' * ( A*x - b ) */
    dplasma_dgemm( parsec, PlasmaTrans, PlasmaNoTrans,
                   1.0, dcA, dcB, 0., subX );

    Rnorm = dplasma_dlange( parsec, PlasmaInfNorm, subX );
    free(subX);

    result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * dcA->n * eps ) ;

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

#endif

GENERATE_F77_BINDINGS (PDGEQRF,
                       pdgeqrf,
                       pdgeqrf_,
                       pdgeqrf__,
                       pdgeqrf_w,
                       (int * M, int * N, double * A, int * IA, int * JA, int * DESCA, double * TAU, double * WORK, int * LWORK, int * info),
                       (M, N, A, IA, JA, DESCA, TAU, WORK, LWORK, info))
