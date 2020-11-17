#ifndef SCALAPACK_WRAPPER
#define SCALAPACK_WRAPPER

#include <stdio.h>
#include <stdlib.h>
#include "../common.h"
#include "../flops.h"
#include "../common_timing.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/redistribute/redistribute_internal.h"

// #define CHECK_RESULTS
// #define WRAPPER_VERBOSE
#define MEASURE_INTERNAL_TIMES
// #define WRAPPER_VERBOSE_CALLS
// #define COUNT_WRAPPED_CALLS

// #define APPLEVEL

#define ACTIVATE_REDIS /* Enable redistribution for submatrixes with offsets not aligned */
#define ACTIVATE_REDIS_SIZE /* Enable redistribution to optimized matrixes with tiles < REDIS_BLOCKSZ_MINxREDIS_BLOCKSZ_MIN*/
#define FORZE_REDIS_SIZE /* Force always redistribution */
/* env vars to set up:
 * - REDIS_BLOCKSZ: block size of the redistributed matrixes
 * - REDIS_BLOCKSZ_MIN: if undefined FORZE_REDIS_SIZE: redistribution when matrix block sizes < REDIS_BLOCKSZ_MIN
 */
extern int REDIS_BLOCKSZ; /*Block size of the redistributed matrixes*/
extern int REDIS_BLOCKSZ_MIN; /*Max block size to redistribute if undefined FORZE_REDIS_SIZE */

#define DO_BLOCKING_REDISTRIBUTION /* redistributions are blocking */
//#undef DO_BLOCKING_REDISTRIBUTION /* redistributions are not blocking */
#ifndef DO_BLOCKING_REDISTRIBUTION
#define ALL_REDIS_NONBLOCKING /* force all redistributions to be not blocking */
#endif

#ifdef APPLEVEL
#define GENERATE_F77_BINDINGS(upper_case,\
                              lower_case,\
                              single_underscore,\
                              double_underscore,\
                              wrapper_function,\
                              signature,\
                              params)\
            void w##single_underscore signature { wrapper_function params; }
#else
#define GENERATE_F77_BINDINGS(upper_case,\
                              lower_case,\
                              single_underscore,\
                              double_underscore,\
                              wrapper_function,\
                              signature,\
                              params)\
            void single_underscore signature { wrapper_function params; }
#endif

extern parsec_context_t* parsec_ctx;

#ifdef COUNT_WRAPPED_CALLS
extern int count_REDIS_IN;
extern int count_REDIS_OUT;
extern int count_PDGEMM;
extern int count_PDGEQRF;
extern int count_PDGETRF_1D;
extern int count_PDGETRF_NOPIV;
extern int count_PDPOTRF;
extern int count_PDTRMM;
extern int count_PDTRSM;
#endif

#ifdef MEASURE_INTERNAL_TIMES
extern double redis_time;
#endif

void parsec_init_wrapped_call(MPI_Comm comm);


#define PASTE_CODE_INIT_LAPACK_MATRIX(DC, TYPE, PTR, INIT_PARAMS)        \
    TYPE##_t DC;                                                         \
    TYPE##_lapack_init INIT_PARAMS;                                      \
    DC.mat = PTR;                                                        \
    parsec_data_collection_set_key((parsec_data_collection_t*)&DC, #DC);


#define PASTE_CODE_INIT_TILED_REDIS_MATRIX(DC, TYPE, INIT_PARAMS)        \
    TYPE##_init INIT_PARAMS;                                             \
    DC.mat = parsec_data_allocate(                                       \
                (size_t)DC.super.nb_local_tiles *                        \
                (size_t)DC.super.bsiz *                                  \
                (size_t)parsec_datadist_getsizeoftype(DC.super.mtype));  \
    parsec_data_collection_set_key((parsec_data_collection_t*)&DC, #DC);



#define WRAPPER_DUMP_TILE_MATRIX(name, dcA, P, Q, parsec_ctx, uplo) do{           \
    PARSEC_DEBUG_VERBOSE(3, parsec_debug_output, "dc %p %s Grid PxQ %dx%d, type %d entire lmxln %dx%d "\
      "|| mxn %dx%d submatrix at %d,%d || tiles mbxnb %dx%d (%d elems) -> lmtxlnt %dx%d total tiles; "\
      "submatrix mtxnt %dx%d; local  submat slmxsln %dx%d  mat llmxlln %dx%d (%d tiles) uplo %s\n",\
    &dcA, name, P, Q, dcA.super.dtype, \
    dcA.super.lm, dcA.super.ln,\
    dcA.super.m, dcA.super.n, dcA.super.i, dcA.super.j,\
    dcA.super.mb, dcA.super.nb,dcA.super.bsiz,\
    dcA.super.lmt, dcA.super.lnt, dcA.super.mt, dcA.super.nt, \
    dcA.super.slm, dcA.super.sln,dcA.super.llm, dcA.super.lln,\
    dcA.super.nb_local_tiles, uplo);\
}while(0)


#ifdef ACTIVATE_REDIS
/* Setting up redistribution:
 * - Lapack matrix correctly defined (real ia, ja).
 * - Tiled redistributed matrix without offset.
 * - Redistribution specifying offset in the incomplete tile.
 * Parsec doesn't operate with offsets not aligned. Wrapper relies
 * on Qinglei redistribution. Data_of will return the beginning
 * of the full tile (not the begining of the valid data of the submatrix),
 * and then redistribution applies offset within the incomplete tile.
 */

#ifdef ACTIVATE_REDIS_SIZE
  #ifdef FORZE_REDIS_SIZE
  #define CHECK_SIZE_REDIS(MATRIX, PARSEC_CTX, REDIS, COMM)\
      REDIS = 1;
  #else
  #define CHECK_SIZE_REDIS(MATRIX, PARSEC_CTX, REDIS, COMM)\
      if (MB_##MATRIX < REDIS_BLOCKSZ_MIN) {\
            REDIS = 1;\
      }
  #endif //#ifdef ACTIVATE_REDIS_SIZE
#else
  #define CHECK_SIZE_REDIS(MATRIX, PARSEC_CTX, REDIS, COMM)
#endif //#ifdef ACTIVATE_REDIS_SIZE

#ifdef COUNT_WRAPPED_CALLS
  #define DO_COUNT_REDIS_IN() count_REDIS_IN++
  #define DO_COUNT_REDIS_OUT() count_REDIS_OUT++
#else
  #define DO_COUNT_REDIS_IN()
  #define DO_COUNT_REDIS_OUT()
#endif //#ifdef COUNT_WRAPPED_CALLS

#ifdef ACTIVATE_REDIS_SIZE
  #define SET_SIZE_REDIS(MATRIX)\
          int cmb = REDIS_BLOCKSZ;\
          int cnb = REDIS_BLOCKSZ;
#else
  #define SET_SIZE_REDIS(MATRIX)\
          int cmb = MB_##MATRIX;\
          int cnb = NB_##MATRIX;
#endif //#ifdef ACTIVATE_REDIS_SIZE

#ifdef MEASURE_INTERNAL_TIMES
  #define REDIS_TIME_INI(COMM)  SYNC_TIME_START_COMM(COMM)
  #define REDIS_TIME_FINI(COMM, MATRIX) \
          SYNC_TIME_STOP_COMM(COMM);\
          if(rank_##MATRIX==0) {\
              printf("[****] TIME(s) %12.5f : REDIS " #MATRIX "\n", sync_time_elapsed);\
              redis_time+=sync_time_elapsed;\
          }
#else
  #define REDIS_TIME_INI(COMM)
  #define REDIS_TIME_FINI(COMM, MATRIX)
#endif //#ifdef MEASURE_INTERNAL_TIMES


#define DEFAULT_REDIS 0
#define BLOCKING_REDIS 1
#define INPUT_REDIS 0
#define OUTPUT_REDIS 1

#define DO_BLOCKING_REDIS(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TARGET)\
        if(REDIS_TARGET == INPUT_REDIS){\
            parsec_redistribute(PARSEC_CTX,\
                              (parsec_tiled_matrix_dc_t *)&dc##MATRIX##_lapack,\
                              (parsec_tiled_matrix_dc_t *)&dc##MATRIX##_redis,\
                              dc##MATRIX##_lapack.super.m, dc##MATRIX##_lapack.super.n,\
                              cI##MATRIX % dc##MATRIX##_lapack.super.mb, cJ##MATRIX % dc##MATRIX##_lapack.super.nb,\
                              0, 0);\
        }else{\
            parsec_redistribute(PARSEC_CTX,\
                              (parsec_tiled_matrix_dc_t *)&dc##MATRIX##_redis,\
                              (parsec_tiled_matrix_dc_t *)&dc##MATRIX##_lapack,\
                              dc##MATRIX##_lapack.super.m, dc##MATRIX##_lapack.super.n,\
                              0, 0,\
                              cI##MATRIX % dc##MATRIX##_lapack.super.mb, cJ##MATRIX % dc##MATRIX##_lapack.super.nb);\
        }

#define DO_ASYNC_REDIS_INI(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TARGET)\
        if(REDIS_TARGET == INPUT_REDIS){\
            TASKPOOL_NAME##MATRIX = parsec_redistribute_New(\
                                    (parsec_tiled_matrix_dc_t *)&dc##MATRIX##_lapack,\
                                    (parsec_tiled_matrix_dc_t *)&dc##MATRIX##_redis,\
                                    dc##MATRIX##_lapack.super.m, dc##MATRIX##_lapack.super.n,\
                                    cI##MATRIX % dc##MATRIX##_lapack.super.mb, cJ##MATRIX % dc##MATRIX##_lapack.super.nb,\
                                    0, 0);\
        }else{\
            TASKPOOL_NAME##MATRIX = parsec_redistribute_New(\
                                    (parsec_tiled_matrix_dc_t *)&dc##MATRIX##_redis,\
                                    (parsec_tiled_matrix_dc_t *)&dc##MATRIX##_lapack,\
                                    dc##MATRIX##_lapack.super.m, dc##MATRIX##_lapack.super.n,\
                                    0, 0,\
                                    cI##MATRIX % dc##MATRIX##_lapack.super.mb, cJ##MATRIX % dc##MATRIX##_lapack.super.nb);\
        }\
        PARSEC_CHECK_ERROR(parsec_context_add_taskpool(PARSEC_CTX, TASKPOOL_NAME##MATRIX), "parsec_context_add_taskpool");\
        parsec_data_collection_set_owner((parsec_data_collection_t *)&dc##MATRIX##_redis, TASKPOOL_NAME##MATRIX);


#ifdef DO_BLOCKING_REDISTRIBUTION
  #define DO_REDIS(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TYPE, REDIS_TARGET) DO_BLOCKING_REDIS(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TARGET)
#else
  #ifdef ALL_REDIS_NONBLOCKING
  #define DO_REDIS(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TYPE, REDIS_TARGET) DO_ASYNC_REDIS_INI(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TARGET)
  #else
  #define DO_REDIS(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TYPE, REDIS_TARGET)\
          if(REDIS_TYPE == BLOCKING_REDIS) {DO_BLOCKING_REDIS(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TARGET);}\
          else {DO_ASYNC_REDIS_INI(PARSEC_CTX, MATRIX, TASKPOOL_NAME, REDIS_TARGET);}
  #endif /* ALL_REDIS_NONBLOCKING */
#endif /* DO_BLOCKING_REDISTRIBUTION */

#define PASTE_CODE_REDIS_INPUT(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)\
    two_dim_block_cyclic_t dc##MATRIX##_redis;\
    two_dim_block_cyclic_t * dc##MATRIX;\
    parsec_taskpool_t* tp_redis_in_##MATRIX = NULL;\
    parsec_taskpool_t* tp_redis_out_##MATRIX = NULL;\
    CHECK_SIZE_REDIS(MATRIX, PARSEC_CTX, REDIS, COMM);\
    if(REDIS){\
        DO_COUNT_REDIS_IN();\
        SET_SIZE_REDIS(MATRIX);\
        PASTE_CODE_INIT_TILED_REDIS_MATRIX(dc##MATRIX##_redis, two_dim_block_cyclic,\
                              (&dc##MATRIX##_redis, matrix_RealDouble, matrix_Tile,\
                               nodes_##MATRIX, rank_##MATRIX,\
                               cmb, cnb,\
                               MATRIX##m, MATRIX##n,\
                               0, 0,\
                               MATRIX##m, MATRIX##n,\
                               KP, KQ,\
                               iP_##MATRIX, jQ_##MATRIX,\
                               P_##MATRIX));\
        REDIS_TIME_INI(COMM);\
        DO_REDIS(PARSEC_CTX, MATRIX, tp_redis_in_, REDIS_TYPE, INPUT_REDIS);\
        REDIS_TIME_FINI(COMM, MATRIX);\
        dc##MATRIX = &dc##MATRIX##_redis;\
    }else{\
        WRAPPER_DUMP_TILE_MATRIX("lapack", dc##MATRIX##_lapack, P_##MATRIX, Q_##MATRIX, PARSEC_CTX, "");\
        dc##MATRIX = &dc##MATRIX##_lapack;\
    }

#define DO_REDIS_OUTPUT(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)\
    if(REDIS){\
        DO_COUNT_REDIS_OUT();\
        REDIS_TIME_INI(COMM);\
        /* not waiting for taskpool fails */\
        DO_REDIS(PARSEC_CTX, MATRIX, tp_redis_out_, REDIS_TYPE, OUTPUT_REDIS);\
        REDIS_TIME_FINI(COMM, MATRIX);\
    }


#ifdef DO_BLOCKING_REDISTRIBUTION
  #define PASTE_CODE_REDIS_OUTPUT_INI(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)
  #define PASTE_CODE_REDIS_OUTPUT_FINI(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE) DO_REDIS_OUTPUT(MATRIX, PARSEC_CTX, REDIS, COMM, BLOCKING_REDIS)
#else
  #define PASTE_CODE_REDIS_OUTPUT_INI(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)  DO_REDIS_OUTPUT(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)
  #define PASTE_CODE_REDIS_OUTPUT_FINI(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)
#endif

#define PASTE_CODE_CLEANUP_REDIS(MATRIX, PARSEC_CTX, REDIS, COMM)\
    if(REDIS){\
        if(tp_redis_in_##MATRIX != NULL)  parsec_redistribute_Destruct(tp_redis_in_##MATRIX);\
        if(tp_redis_out_##MATRIX != NULL) parsec_redistribute_Destruct(tp_redis_out_##MATRIX);\
        parsec_data_free(dc##MATRIX##_redis.mat);\
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dc##MATRIX##_redis);\
        dc##MATRIX = &dc##MATRIX##_lapack;\
    }

#else
  #define PASTE_CODE_REDIS_INPUT(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)\
          two_dim_block_cyclic_t * dc##MATRIX = &dc##MATRIX##_lapack;
  #define PASTE_CODE_REDIS_OUTPUT_INI(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)
  #define PASTE_CODE_REDIS_OUTPUT_FINI(MATRIX, PARSEC_CTX, REDIS, COMM, REDIS_TYPE)
  #define PASTE_CODE_CLEANUP_REDIS(MATRIX, PARSEC_CTX, REDIS, COMM)
#endif /*ACTIVATE_REDIS*/



#ifdef MEASURE_INTERNAL_TIMES
    #define INI_WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(PARSEC, KERNEL, PARAMS, DESTRUCT, rank, P, Q, NB, N, COMM)\
        SYNC_TIME_START_COMM(COMM);\
        parsec_taskpool_t* PARSEC_##KERNEL = dplasma_##KERNEL##_New PARAMS;\
        PARSEC_CHECK_ERROR(parsec_context_add_taskpool(PARSEC, PARSEC_##KERNEL), "parsec_context_add_taskpool");\
        SYNC_TIME_STOP_COMM(COMM);\
        double stime_A = sync_time_elapsed;\


    #define FINI_WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(PARSEC, KERNEL, PARAMS, DESTRUCT, rank, P, Q, NB, N, COMM)\
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
                      " - ENQ %12.5f - DEST %12.5f | REDIS %12.5f \n",\
                              stime_B, P, Q, NB, N,\
                              gflops=(flops/1e9)/stime_B,\
                              (stime_A+stime_B+stime_C),\
                              (flops/1e9)/(stime_A+stime_B+stime_C),\
                              stime_A,stime_C, redis_time);\
        }\
        redis_time=0.0;
#else

    #define INI_WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(PARSEC, KERNEL, PARAMS, DESTRUCT, rank, P, Q, NB, N, COMM)\
        parsec_taskpool_t* PARSEC_##KERNEL = dplasma_##KERNEL##_New PARAMS;\
        PARSEC_CHECK_ERROR(parsec_context_add_taskpool(PARSEC, PARSEC_##KERNEL), "parsec_context_add_taskpool");


    #define FINI_WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(PARSEC, KERNEL, PARAMS, DESTRUCT, rank, P, Q, NB, N, COMM)\
        PARSEC_CHECK_ERROR(parsec_context_start(PARSEC), "parsec_context_start");\
        PARSEC_CHECK_ERROR(parsec_context_wait(PARSEC), "parsec_context_wait");\
        DESTRUCT;


    #define WRAPPER_PASTE_CODE_ENQUEUE_PROGRESS_DESTRUCT_KERNEL(PARSEC, KERNEL, PARAMS, DESTRUCT, rank, P, Q, NB, N, COMM)\
        parsec_taskpool_t* PARSEC_##KERNEL = dplasma_##KERNEL##_New PARAMS;\
        PARSEC_CHECK_ERROR(parsec_context_add_taskpool(PARSEC, PARSEC_##KERNEL), "parsec_context_add_taskpool");\
        PARSEC_CHECK_ERROR(parsec_context_start(PARSEC), "parsec_context_start");\
        PARSEC_CHECK_ERROR(parsec_context_wait(PARSEC), "parsec_context_wait");\
        DESTRUCT;
#endif

#ifdef WRAPPER_VERBOSE
    #define PRINT(parsec_ctx, comm, uplo_parsec, mname, dc, P, Q)do{\
        MPI_Barrier(comm);\
        WRAPPER_DUMP_TILE_MATRIX(mname, dc, P, Q, parsec_ctx, "");\
        /*dplasma_dprint( parsec_ctx, uplo_parsec, (parsec_tiled_matrix_dc_t *)&dc );*/\
        MPI_Barrier(comm);\
    }while(0)
#else
    #define PRINT
#endif

static void dcopy_lapack_tile(parsec_context_t *parsec,
                              two_dim_block_cyclic_t *dcIn,
                              two_dim_block_cyclic_t *dcOut,
                              int mloc, int nloc){

    int LD_lapack = dcIn->super.llm;
    int LD_tile   = dcOut->super.mb;

    /*# of local row tiles A: mloc != llm because llm is stored not submatrix */
    int lrowtiles = (mloc % dcIn->super.mb == 0)? mloc/dcIn->super.mb: (mloc/dcIn->super.mb) + 1;
    int lcoltiles = (nloc % dcIn->super.nb == 0)? nloc/dcIn->super.nb: (nloc/dcIn->super.nb) + 1;
    int m, n;
    for(m=0; m< dcIn->nb_elem_r; m++){
        for(n=0; n< dcIn->nb_elem_c; n++){
            int tempmm = ( m == (lrowtiles-1) ) ? mloc - m*dcIn->super.mb : dcIn->super.mb;
            int tempnn = ( n == (lcoltiles-1) ) ? nloc - n*dcIn->super.nb : dcIn->super.nb;
            int start_block_lapack = m*dcIn->super.mb   + n*dcIn->super.nb*LD_lapack;
            int start_block_tile   = m*dcIn->super.bsiz + n*dcIn->super.bsiz*dcIn->nb_elem_r;
            int ii, jj;
            for(jj=0; jj<tempnn; jj++) {
                for(ii=0; ii<tempmm; ii++) {
                    int ind_lapack = start_block_lapack + jj*LD_lapack + ii;
                    int ind_tile   = start_block_tile   + jj*LD_tile   + ii;
                    // printf("COPY ind lapack %d (%d + %d) to ind_tile %d (%d + %d)\n",
                    //     ind_lapack, start_block_lapack, jj*LD_lapack + ii,
                    //     ind_tile, start_block_tile, jj*LD_tile   + ii);
                    if(dcIn->super.mtype == matrix_RealDouble){
                        ((double*)dcOut->mat)[ind_tile] = ((double*)dcIn->mat)[ind_lapack];
                    }else if(dcIn->super.mtype == matrix_Integer){
                        ((int*)dcOut->mat)[ind_tile] = ((int*)dcIn->mat)[ind_lapack];
                    }else{
                        assert(1==0);
                    }
                }
            }
        }
    }
}

static void dcopy_lapack_lapack(parsec_context_t *parsec,
                              two_dim_block_cyclic_t *dcIn,
                              two_dim_block_cyclic_t *dcOut,
                              int mloc, int nloc){

    int LD_lapack = dcIn->super.llm;
    int LD_tile = dcOut->super.mb;

    /*# of local row tiles A: mloc != llm because llm is stored not submatrix */
    int lrowtiles = (mloc % dcIn->super.mb == 0)? mloc/dcIn->super.mb: (mloc/dcIn->super.mb) + 1;
    int lcoltiles = (nloc % dcIn->super.nb == 0)? nloc/dcIn->super.nb: (nloc/dcIn->super.nb) + 1;
    int m, n;
    for(m=0; m< dcIn->nb_elem_r; m++){
        for(n=0; n< dcIn->nb_elem_c; n++){
            int tempmm = ( m == (lrowtiles-1) ) ? mloc - m*dcIn->super.mb : dcIn->super.mb;
            int tempnn = ( n == (lcoltiles-1) ) ? nloc - n*dcIn->super.nb : dcIn->super.nb;
            int ii, jj;
            for(ii=0; ii<tempmm; ii++) {
                for(jj=0; jj<tempnn; jj++) {
                    int ind_lapack = m*dcIn->super.mb + n*dcIn->super.nb*LD_lapack
                            + jj*LD_lapack + ii;
                    if(dcIn->super.mtype == matrix_RealDouble){
                        ((double*)dcOut->mat)[ind_lapack] = ((double*)dcIn->mat)[ind_lapack];
                    }else if(dcIn->super.mtype == matrix_Integer){
                        ((int*)dcOut->mat)[ind_lapack] = ((int*)dcIn->mat)[ind_lapack];
                    }else{
                        assert(1==0);
                    }
                }
            }
        }
    }
}

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
 int LDD_##MATRIX = DESC##MATRIX[WRAPPER_LLD1_];\
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

extern MPI_Comm Cblacs2sys_handle(int BlacsCtxt);
extern void     Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);
extern void     Cblacs_get(int ConTxt, int what, int *val);
extern int      numroc_(int *n, int *nb, int *iproc, int *srcproc, int *nprocs);
#endif
