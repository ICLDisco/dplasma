#ifndef DPLASMA_LAPACK_ADT_H_HAS_BEEN_INCLUDED
#define DPLASMA_LAPACK_ADT_H_HAS_BEEN_INCLUDED

/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */

#include "dplasma.h"
#include "parsec/arena.h"
#include "parsec/class/parsec_hash_table.h"

/**************************/
/* Utils for LAPACK-TILED */
/**************************/

#define DFULL_T  0
#define DC00     1
#define DCM0     2
#define DC0N     3
#define DCMN     4
#define DTOP     5
#define DBOTTOM  6
#define DLEFT    7
#define DRIGHT   8
#define LOC_SZ   9

#define LAPACK     0
#define TILED      1
#define MAX_LAYOUT 2

static inline int dplasma_tile_location(const parsec_tiled_matrix_t *dM, int m, int n)
{

  int ret =
      ((m == 0)                      && (n == 0))                      ? DC00 :\
      (            (m == (dM->mt-1)) &&             (n == (dM->nt-1))) ? DCMN :\
      ((m == 0)                      &&             (n == (dM->nt-1))) ? DC0N :\
      (            (m == (dM->mt-1)) && (n == 0))                      ? DCM0 :\
      ((m == 0)                      && (n != 0) && (n != (dM->nt-1))) ? DTOP :\
      (            (m == (dM->mt-1)) && (n != 0) && (n != (dM->nt-1))) ? DBOTTOM :\
      ((m != 0) && (m != (dM->mt-1)) && (n == 0))                      ? DLEFT :\
      ((m != 0) && (m != (dM->mt-1)) &&             (n == (dM->nt-1))) ? DRIGHT :\
      ((m != 0) && (m != (dM->mt-1)) && (n != 0) && (n != (dM->nt-1))) ? DFULL_T : -1;

  PARSEC_DEBUG_VERBOSE(27, parsec_debug_output, "LOC %p(%d,%d) %d ", dM, m, n, ret);
  return ret;
}

/*************************************************************/
/* Management of LAPACK-TILED datatypes info and translation */
/*************************************************************/
#define LAPACK_ADT_KEY_STR_SZ 100

extern parsec_hash_table_t *dplasma_datatypes_lapack_helper;

typedef struct lapack_info_s {
    int lda;
    int rows;
    int cols;
    int loc;
    int shape;
    int layout;
} lapack_info_t;

struct dplasma_datatype_lapack_helper_s {
    parsec_hash_table_item_t ht_item;
    lapack_info_t info;
    parsec_arena_datatype_t adt;
};
typedef struct dplasma_datatype_lapack_helper_s dplasma_datatype_lapack_helper_t;

int dplasma_set_datatype_info( const dplasma_data_collection_t *dc, parsec_arena_datatype_t adt,
                               lapack_info_t info);

int dplasma_get_info_from_datatype( const dplasma_data_collection_t *dc, parsec_datatype_t dtt,
                                    lapack_info_t *info, parsec_arena_datatype_t **adt);
int dplasma_get_datatype_from_info( const dplasma_data_collection_t *dc,
                                    lapack_info_t *info, parsec_arena_datatype_t **adt);

int dplasma_cleanup_datatype_info( const dplasma_data_collection_t *dc,
                                   lapack_info_t info, parsec_arena_datatype_t *adt);

/***************************************************************************/
/* Management of dplasma_data_collection_t to wrap parsec_datacollections */
/***************************************************************************/

typedef struct dplasma_data_collection_s {
    parsec_data_collection_t super;
    parsec_tiled_matrix_t *dc_original;
} dplasma_data_collection_t;

dplasma_data_collection_t *dplasma_wrap_data_collection(parsec_tiled_matrix_t *A);
void dplasma_unwrap_data_collection(dplasma_data_collection_t *ddc);

#endif  /* DPLASMA_LAPACK_ADT_H_HAS_BEEN_INCLUDED */

