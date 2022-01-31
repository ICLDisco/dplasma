#ifndef SCALAPACK_WRAPPER
#define SCALAPACK_WRAPPER

#include <stdio.h>
#include <stdlib.h>
#include "../../tests/common.h"
#include "../../tests/common_timing.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/redistribute/redistribute_internal.h"

/* Performed DPLASMA tests checking at the end of wrapper. */
// #define CHECK_RESULTS

/* Verbose modes for the wrapped calls. */
/* Print call with parameters. */
// #define WRAPPER_VERBOSE_CALLS
/* Print matrices info/data. */
// #define WRAPPER_VERBOSE_FULL

/* Measure times of the wrapped calls. */
#define MEASURE_INTERNAL_TIMES

/* Count the number of wrapped calls by type and number of redistributions performed. */
#define COUNT_WRAPPED_CALLS

/* Generate wrappers with signature w<function_call>: allow wrapping only
 * application level calls if they have been modified to w<function_call>.
 */
// #define APPLEVEL

/* Activate redistribution to enable operation (e.g. submatrixes with offsets not aligned).
 * Mandatory for TRSM. */
#define ACTIVATE_REDIS

/* Activate redistribution to optimized matrices with tiles inferior to REDIS_BLOCKSZ_MINxREDIS_BLOCKSZ_MIN*/
// #define ACTIVATE_REDIS_SIZE

/* Force always redistribution for every input matrix */
// #define FORCE_REDIS

/* Environ variable to control redistribution:
 * - REDIS_BLOCKSZ: block size of the redistributed matrixes. Default value: 512.
 * - REDIS_BLOCKSZ_MIN: if defined ACTIVATE_REDIS_SIZE: redistribution when matrix
     block sizes < REDIS_BLOCKSZ_MIN. Default value: 128.
 */
#if defined(ACTIVATE_REDIS_SIZE) || defined(FORCE_REDIS)
#define ACTIVATE_REDIS
#endif
/* Setting up redistribution:
 * - Lapack matrix correctly defined (real ia, ja).
 * - Tiled redistributed matrix without offsets.
 * - Redistribution specifying offset in the incomplete tile.
 * Parsec doesn't operate with offsets not aligned. Wrapper relies
 * on redistribution. Data_of will return the beginning
 * of the full tile (not the begining of the valid data of the submatrix),
 * and then redistribution applies offset within the incomplete tile.
 */

extern int REDIS_BLOCKSZ; /*Block size of the redistributed matrixes*/
extern int REDIS_BLOCKSZ_MIN; /*Max block size to redistribute if defined ACTIVATE_REDIS_SIZE */

#ifdef APPLEVEL
/* generate wrappers with signature w<function_call>*/
#ifdef SCALAPACK_SUP_UNDERSCORE
#define GENERATE_F77_BINDINGS(upper_case, lower_case, single_underscore, double_underscore,\
                              wrapper_function, signature, params)\
            void w##lower_case signature { wrapper_function params; }
#else
#define GENERATE_F77_BINDINGS(upper_case, lower_case, single_underscore, double_underscore,\
                              wrapper_function, signature, params)\
            void w##single_underscore signature { wrapper_function params; }
#endif
#else
#ifdef SCALAPACK_SUP_UNDERSCORE
#define GENERATE_F77_BINDINGS(upper_case, lower_case, single_underscore, double_underscore,\
                              wrapper_function, signature, params)\
            void lower_case signature { wrapper_function params; }
#else
#define GENERATE_F77_BINDINGS(upper_case, lower_case, single_underscore, double_underscore,\
                              wrapper_function, signature, params)\
            void single_underscore signature { wrapper_function params; }
#endif
#endif

extern parsec_context_t* parsec_ctx;

#ifdef COUNT_WRAPPED_CALLS
extern int count_REDIS_IN;
extern int count_REDIS_OUT;
extern int count_PDGEMM;
extern int count_PDLATSQR;
extern int count_PDGETRF_1D;
extern int count_PDGETRF_NOPIV;
extern int count_PDPOTRF;
extern int count_PDTRMM;
extern int count_PDTRSM;
#define DO_COUNT_REDIS_IN() count_REDIS_IN++
#define DO_COUNT_REDIS_OUT() count_REDIS_OUT++
#else
#define DO_COUNT_REDIS_IN()
#define DO_COUNT_REDIS_OUT()
#endif

#ifdef MEASURE_INTERNAL_TIMES
extern double redis_time;
#endif

void parsec_init_wrapped_call(MPI_Comm comm);
void parsec_wrapper_devices_release_memory_(void);
void parsec_wrapper_devices_reset_load_(void);

#ifdef MEASURE_INTERNAL_TIMES
  #define REDIS_TIME_INI(COMM)  SYNC_TIME_START_COMM(COMM)
  #define REDIS_TIME_FINI(COMM, RANK, MATRIX) \
          SYNC_TIME_STOP_COMM(COMM);\
          if(RANK == 0) {\
              printf("[****] TIME(s) %12.5f : REDIS %s \n", sync_time_elapsed, MATRIX);\
              redis_time+=sync_time_elapsed;\
          }
#else
  #define REDIS_TIME_INI(COMM)
  #define REDIS_TIME_FINI(COMM, RANK, MATRIX)
#endif //#ifdef MEASURE_INTERNAL_TIMES

parsec_matrix_block_cyclic_t *
redistribute_lapack_input_internal(parsec_matrix_block_cyclic_t * dc_lapack,
                                   int do_redis, MPI_Comm comm, int rank, char* name,
                                   int P, int Q);

parsec_matrix_block_cyclic_t *
redistribute_lapack_input(parsec_matrix_block_cyclic_t * dc_lapack,
                          int do_redis, MPI_Comm comm, int rank, char* name);

parsec_matrix_block_cyclic_t *
redistribute_lapack_input_1D(parsec_matrix_block_cyclic_t * dc_lapack,
                          int do_redis, MPI_Comm comm, int rank, char* name, int redisP, int redisQ);

parsec_matrix_block_cyclic_t *
redistribute_lapack_output_cleanup(parsec_matrix_block_cyclic_t * dc_lapack, parsec_matrix_block_cyclic_t * dc_redis,
                           int do_redis, MPI_Comm comm, int rank, char* name);



#ifdef MEASURE_INTERNAL_TIMES
    #define WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(PARSEC, KERNEL, PARAMS, DESTRUCT, rank, P, Q, NB, N, COMM)\
        SYNC_TIME_START_COMM(COMM);\
        parsec_taskpool_t* PARSEC_##KERNEL = dplasma_##KERNEL##_New PARAMS;\
        PARSEC_CHECK_ERROR(parsec_context_add_taskpool(PARSEC, PARSEC_##KERNEL), "parsec_context_add_taskpool");\
        SYNC_TIME_STOP_COMM(COMM);\
        double stime_A = sync_time_elapsed;\
        SYNC_TIME_START_COMM(COMM);\
        PARSEC_CHECK_ERROR(parsec_context_start(PARSEC), "parsec_context_start");\
        TIME_START();\
        PARSEC_CHECK_ERROR(parsec_context_wait(PARSEC), "parsec_context_wait");\
        SYNC_TIME_STOP_COMM(COMM);\
        double stime_B = sync_time_elapsed;\
        SYNC_TIME_START_COMM(COMM);\
        DESTRUCT;\
        SYNC_TIME_STOP_COMM(COMM);\
        double stime_C = sync_time_elapsed;\
        if(rank==0){\
            printf("[****] TIME(s) %12.5f : " #KERNEL "\tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops"\
                      " - ENQ&PROG&DEST %12.5f : %14f gflops"\
                      " - ENQ %12.5f - DEST %12.5f | REDIS %12.5f | REDIS_NB= %4d \n",\
                              stime_B, P, Q, NB, N,\
                              gflops=(flops/1e9)/stime_B,\
                              (stime_A+stime_B+stime_C),\
                              (flops/1e9)/(stime_A+stime_B+stime_C),\
                              stime_A,stime_C, redis_time, \
                              REDIS_BLOCKSZ);\
        }\
        redis_time=0.0;
#else
    #define WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(PARSEC, KERNEL, PARAMS, DESTRUCT, rank, P, Q, NB, N, COMM)\
        parsec_taskpool_t* PARSEC_##KERNEL = dplasma_##KERNEL##_New PARAMS;\
        PARSEC_CHECK_ERROR(parsec_context_add_taskpool(PARSEC, PARSEC_##KERNEL), "parsec_context_add_taskpool");\
        PARSEC_CHECK_ERROR(parsec_context_start(PARSEC), "parsec_context_start");\
        PARSEC_CHECK_ERROR(parsec_context_wait(PARSEC), "parsec_context_wait");\
        DESTRUCT;
#endif



#define WRAPPER_DUMP_TILE_MATRIX(name, dcA, parsec_ctx, uplo) do{           \
    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output, "dc %p %s Grid PxQ %dx%d, type %d entire lmxln %dx%d "\
      "|| mxn %dx%d submatrix at %d,%d || tiles mbxnb %dx%d (%d elems) -> lmtxlnt %dx%d total tiles; "\
      "submatrix mtxnt %dx%d; local  submat slmxsln %dx%d  mat llmxlln %dx%d (%d tiles) uplo %s\n",\
    dcA, name, dcA->grid.rows, dcA->grid.cols, dcA->super.dtype, \
    dcA->super.lm, dcA->super.ln,\
    dcA->super.m, dcA->super.n, dcA->super.i, dcA->super.j,\
    dcA->super.mb, dcA->super.nb,dcA->super.bsiz,\
    dcA->super.lmt, dcA->super.lnt, dcA->super.mt, dcA->super.nt, \
    dcA->super.slm, dcA->super.sln,dcA->super.llm, dcA->super.lln,\
    dcA->super.nb_local_tiles, uplo);\
}while(0)

#ifdef WRAPPER_VERBOSE_FULL
    #define PRINT(parsec_ctx, comm, uplo_parsec, mname, dc) do {\
        MPI_Barrier(comm);\
        WRAPPER_DUMP_TILE_MATRIX(mname, dc, parsec_ctx, "");\
        /*dplasma_dprint( parsec_ctx, uplo_parsec, (parsec_tiled_matrix_t *)&dc );*/\
        MPI_Barrier(comm);\
    }while(0)
#else
    #define PRINT(parsec_ctx, comm, uplo_parsec, mname, dc)
#endif


void dcopy_lapack_tile(parsec_context_t *parsec,
                       parsec_matrix_block_cyclic_t *dcIn,
                       parsec_matrix_block_cyclic_t *dcOut,
                       int mloc, int nloc);

void dcopy_lapack_lapack(parsec_context_t *parsec,
                         parsec_matrix_block_cyclic_t *dcIn,
                         parsec_matrix_block_cyclic_t *dcOut,
                         int mloc, int nloc);

#define OP_TRANS(TRANS) ( ((TRANS=='N')||(TRANS=='n')) ? dplasmaNoTrans : \
                          ((TRANS=='T')||(TRANS=='t')) ? dplasmaTrans : \
                          ((TRANS=='C')||(TRANS=='c')) ? dplasmaConjTrans : -1 )


#define OP_SIDE(SIDE) ( ((SIDE=='L')||(SIDE=='l')) ? dplasmaLeft : \
                        ((SIDE=='R')||(SIDE=='r')) ? dplasmaRight : -1 )

#define OP_UPLO(UPLO) ( ((UPLO=='U')||(UPLO=='u')) ? dplasmaUpper : \
                        ((UPLO=='L')||(UPLO=='l')) ? dplasmaLower : -1 )

#define OP_DIAG(DIAG) ( ((DIAG=='U')||(DIAG=='u')) ? dplasmaUnit : \
                        ((DIAG=='N')||(DIAG=='n')) ? dplasmaNonUnit : -1 )


#define PASTE_SETUP(MATRIX)\
 int ictxt_##MATRIX = DESC##MATRIX[WRAPPER_CTXT1_];\
 int gM_##MATRIX = DESC##MATRIX[WRAPPER_M1_];\
 int gN_##MATRIX = DESC##MATRIX[WRAPPER_N1_];\
 int MB_##MATRIX = DESC##MATRIX[WRAPPER_MB1_];\
 int NB_##MATRIX = DESC##MATRIX[WRAPPER_NB1_];\
 int LLD_##MATRIX = DESC##MATRIX[WRAPPER_LLD1_];\
 \
 int cI##MATRIX = (*I##MATRIX)-1;\
 int cJ##MATRIX = (*J##MATRIX)-1;\
 \
 int myrow_##MATRIX, mycol_##MATRIX;\
 int P_##MATRIX, Q_##MATRIX;\
 int iP_##MATRIX = DESC##MATRIX[WRAPPER_RSRC1_];\
 int jQ_##MATRIX = DESC##MATRIX[WRAPPER_CSRC1_];\
 Cblacs_gridinfo( ictxt_##MATRIX, &P_##MATRIX, &Q_##MATRIX, &myrow_##MATRIX, &mycol_##MATRIX );\
 int mloc_##MATRIX = numroc_( &gM_##MATRIX, &MB_##MATRIX, &myrow_##MATRIX, &DESC##MATRIX[WRAPPER_RSRC1_], &P_##MATRIX );\
 int nloc_##MATRIX = numroc_( &gN_##MATRIX, &NB_##MATRIX, &mycol_##MATRIX, &DESC##MATRIX[WRAPPER_CSRC1_], &Q_##MATRIX );\
 PARSEC_DEBUG_VERBOSE(3, parsec_debug_output,  "" #MATRIX "Cblacs_gridinfo gMxgN %dx%d nbxmb %dx%d PxQ %dx%d rowxcol %dx%d  mloc %d nloc %d\n", \
    gM_##MATRIX, gN_##MATRIX, MB_##MATRIX, NB_##MATRIX, P_##MATRIX, Q_##MATRIX, myrow_##MATRIX, mycol_##MATRIX, mloc_##MATRIX, nloc_##MATRIX);\
 \
 int comm_index_##MATRIX;\
 Cblacs_get(ictxt_##MATRIX, 10, &comm_index_##MATRIX);\
 MPI_Comm comm_##MATRIX = Cblacs2sys_handle(comm_index_##MATRIX);\
 int rank_##MATRIX,nodes_##MATRIX;\
 MPI_Comm_rank(comm_##MATRIX, &rank_##MATRIX);\
 MPI_Comm_size(comm_##MATRIX, &nodes_##MATRIX);\
 (void)LLD_##MATRIX;(void)cI##MATRIX;(void)cJ##MATRIX;\
 (void)iP_##MATRIX;(void)jQ_##MATRIX;(void)mloc_##MATRIX;(void)nloc_##MATRIX;


// ScaLAPACK routines and constants

//#define SGET_SYSCONTXT    0
//#define SGET_MSGIDS       1
//#define SGET_DEBUGLVL     2
//#define SGET_BLACSCONTXT 10
#define WRAPPER_SGET_BLACSCONTXT 10
//#define SGET_NR_BS       11
//#define SGET_NB_BS       12
//#define SGET_NR_CO       13
//#define SGET_NB_CO       14
//#define SGET_TOPSREPEAT  15
//#define SGET_TOPSCOHRNT  16


/* ScalaPACK: indexes of array descriptor of a distributed matrix
 */
#define    WRAPPER_DTYPE1_             0                   /* Descriptor Type */
#define    WRAPPER_CTXT1_              1                     /* BLACS context */
#define    WRAPPER_M1_                 2             /* Global Number of Rows */
#define    WRAPPER_N1_                 3          /* Global Number of Columns */
#define    WRAPPER_MB1_                4                 /* Row Blocking Size */
#define    WRAPPER_NB1_                5              /* Column Blocking Size */
#define    WRAPPER_RSRC1_              6            /* Starting Processor Row */
#define    WRAPPER_CSRC1_              7         /* Starting Processor Column */
#define    WRAPPER_LLD1_               8           /* Local Leading Dimension */
#define    WRAPPER_DLEN1_              9                 /* Descriptor Length */

#ifdef SCALAPACK_SUP_UNDERSCORE
#define numroc_ numroc
#endif

extern MPI_Comm Cblacs2sys_handle(int BlacsCtxt);
extern void     Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);
extern void     Cblacs_get(int ConTxt, int what, int *val);
extern int      numroc_(int *n, int *nb, int *iproc, int *srcproc, int *nprocs);
#endif
