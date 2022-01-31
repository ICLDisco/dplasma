#include "common.h"

parsec_context_t* parsec_ctx = NULL;

int REDIS_BLOCKSZ = 512; /*Block size of the redistributed matrixes*/
int REDIS_BLOCKSZ_MIN = 128; /*Max block size to redistribute if defined ACTIVATE_REDIS_SIZE */

#ifdef MEASURE_INTERNAL_TIMES
double time_elapsed = 0.0;
double sync_time_elapsed = 0.0;
double redis_time = 0.0;
#endif

#ifdef COUNT_WRAPPED_CALLS
int count_REDIS_IN = 0;
int count_REDIS_OUT = 0;
int count_PDGEMM = 0;
int count_PDLATSQR = 0;
int count_PDGETRF_1D = 0;
int count_PDGETRF_NOPIV = 0;
int count_PDPOTRF = 0;
int count_PDTRMM = 0;
int count_PDTRSM = 0;
#endif

parsec_matrix_block_cyclic_t *
redistribute_lapack_input_internal(parsec_matrix_block_cyclic_t * dc_lapack,
                                   int do_redis, MPI_Comm comm, int rank, char* name,
                                   int P, int Q)
{
    parsec_matrix_block_cyclic_t * dc_out = dc_lapack;
#ifdef ACTIVATE_REDIS

    #ifdef FORCE_REDIS
      do_redis = 1;
    #else
      #ifdef ACTIVATE_REDIS_SIZE
          if (dc_lapack->super.mb < REDIS_BLOCKSZ_MIN) {
              do_redis = 1;
          }
      #endif
      if ( (dc_lapack->grid.rows != P) || (dc_lapack->grid.cols != Q) ){
          do_redis = 1;
      }
    #endif

    if(do_redis)
    {
        DO_COUNT_REDIS_IN();

        int cmb = REDIS_BLOCKSZ;
        int cnb = REDIS_BLOCKSZ;
        dc_out = (parsec_matrix_block_cyclic_t*)malloc(sizeof(parsec_matrix_block_cyclic_t));
        parsec_matrix_block_cyclic_init(dc_out, dc_lapack->super.mtype, PARSEC_MATRIX_TILE,
                                  dc_lapack->super.super.myrank,
                                  cmb, cnb,
                                  dc_lapack->super.m,    dc_lapack->super.n,
                                  0, 0,
                                  dc_lapack->super.m,    dc_lapack->super.n,
                                  P, Q, /*dc_lapack->grid.rows,  dc_lapack->grid.cols,*/
                                  1, 1, /*dc_lapack->grid.krows, dc_lapack->grid.kcols,*/
                                  0, 0  /*dc_lapack->grid.ip,    dc_lapack->grid.jq*/
                                  );

        dc_out->mat = parsec_data_allocate(
                    (size_t)dc_out->super.nb_local_tiles *
                    (size_t)dc_out->super.bsiz *
                    (size_t)parsec_datadist_getsizeoftype(dc_out->super.mtype));
        if(NULL != name) {
            parsec_data_collection_set_key((parsec_data_collection_t*)&dc_out, name);
        }

        REDIS_TIME_INI(comm);

        parsec_redistribute(parsec_ctx,
                            (parsec_tiled_matrix_t *)dc_lapack,
                            (parsec_tiled_matrix_t *)dc_out,
                            dc_lapack->super.m, dc_lapack->super.n,
                            dc_lapack->super.i % dc_lapack->super.mb, dc_lapack->super.j % dc_lapack->super.nb,
                            0, 0);

        REDIS_TIME_FINI(comm, rank, name);
    }
#endif
    return dc_out;
}

parsec_matrix_block_cyclic_t *
redistribute_lapack_input(parsec_matrix_block_cyclic_t * dc_lapack,
                          int do_redis, MPI_Comm comm, int rank, char* name)
{
    return redistribute_lapack_input_internal(dc_lapack,
                                              do_redis, comm, rank, name,
                                              dc_lapack->grid.rows, dc_lapack->grid.cols);
}


parsec_matrix_block_cyclic_t *
redistribute_lapack_input_1D(parsec_matrix_block_cyclic_t * dc_lapack,
                          int do_redis, MPI_Comm comm, int rank, char* name, int redisP, int redisQ)
{
    return redistribute_lapack_input_internal(dc_lapack,
                                              do_redis, comm, rank, name,
                                              redisP, redisQ);
}

parsec_matrix_block_cyclic_t *
redistribute_lapack_output_cleanup(parsec_matrix_block_cyclic_t * dc_lapack, parsec_matrix_block_cyclic_t * dc_redis,
                           int do_redis, MPI_Comm comm, int rank, char* name)
{
#ifdef ACTIVATE_REDIS
    if(dc_lapack != dc_redis){
      if(do_redis){
          /* redistribute */
          DO_COUNT_REDIS_OUT();
          REDIS_TIME_INI(comm);
          parsec_redistribute(parsec_ctx,
                              (parsec_tiled_matrix_t *)dc_redis,
                              (parsec_tiled_matrix_t *)dc_lapack,
                              dc_lapack->super.m, dc_lapack->super.n,
                              0, 0,
                              dc_lapack->super.i % dc_lapack->super.mb, dc_lapack->super.j % dc_lapack->super.nb);
          REDIS_TIME_FINI(comm, rank, name);
        }
        /* Clean up */
        parsec_data_free(dc_redis->mat);
        parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)dc_redis);
    }
#endif
    return dc_lapack;
}


void
dcopy_lapack_tile(parsec_context_t *parsec,
                  parsec_matrix_block_cyclic_t *dcIn,
                  parsec_matrix_block_cyclic_t *dcOut,
                  int mloc, int nloc){
    (void)parsec;
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
                    if(dcIn->super.mtype == PARSEC_MATRIX_DOUBLE){
                        ((double*)dcOut->mat)[ind_tile] = ((double*)dcIn->mat)[ind_lapack];
                    }else if(dcIn->super.mtype == PARSEC_MATRIX_INTEGER){
                        ((int*)dcOut->mat)[ind_tile] = ((int*)dcIn->mat)[ind_lapack];
                    }else{
                        assert(1==0);
                    }
                }
            }
        }
    }
}

void
dcopy_lapack_lapack(parsec_context_t *parsec,
                    parsec_matrix_block_cyclic_t *dcIn,
                    parsec_matrix_block_cyclic_t *dcOut,
                    int mloc, int nloc){

    int LD_lapack = dcIn->super.llm;
    (void)parsec;

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
                    if(dcIn->super.mtype == PARSEC_MATRIX_DOUBLE){
                        ((double*)dcOut->mat)[ind_lapack] = ((double*)dcIn->mat)[ind_lapack];
                    }else if(dcIn->super.mtype == PARSEC_MATRIX_INTEGER){
                        ((int*)dcOut->mat)[ind_lapack] = ((int*)dcIn->mat)[ind_lapack];
                    }else{
                        assert(1==0);
                    }
                }
            }
        }
    }
}
