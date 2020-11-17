#ifndef DPLASMA_LAPACK_DATATYPE_H_HAS_BEEN_INCLUDED
#define DPLASMA_LAPACK_DATATYPE_H_HAS_BEEN_INCLUDED

/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */
#include "dplasma.h"

#include "types.h"
#include "utils/dplasma_lapack_adtt.h"

static inline
void SET_UP_ARENA_DTT(parsec_tiled_matrix_dc_t * desc,
                      parsec_datatype_t parsec_type,
                      int uplo, int diag,
                      parsec_arena_datatype_t *adt,
                      unsigned int rows, unsigned int cols, unsigned int ld,
                      int location, int shape, int layout, char* name)
{
  (void)name;
  DPLASMA_GET_OR_CONSTRUCT(adt, parsec_type,
                           PARSEC_ARENA_ALIGNMENT_SSE,
                           uplo, diag,
                           rows, cols, ld,
                           -1/* resized = -1 for all dplasma types */);
  dplasma_set_datatype_info(desc, *adt, ld, rows, cols, location, shape, layout);
  char dtt_name[MPI_MAX_OBJECT_NAME]="NULL";
  int len;
  if((adt->opaque_dtt != NULL) && (adt->opaque_dtt != PARSEC_DATATYPE_NULL))MPI_Type_get_name(adt->opaque_dtt, dtt_name, &len);
  PARSEC_DEBUG_VERBOSE(27, parsec_debug_output, "%d SETUP ADT dtt %p %30s [%s]", getpid(), adt->opaque_dtt, dtt_name, name);
}


/* NOTE: we are operating over the following condition:
 * For a given datacollection, if we reuse the datatype for one shape on different
 * locations, then other shapes will also be reusing a datatype for those locations.
 *
 * We are enabling lookup by two keys:
 *    shape + loc        --> parsec_datatype_t*
 *          This is fine, the application request shape & locations
 *    parsec_datatype_t* --> lda, loc, shape, adt
 *          This is set once by *dtt for the first location and shape it is used
 *
 * Therefore when translating a dtt to a different shape: we will first get the current
 * shape & location, and then obtain the target shape on that location.
 *
 * For the usage we are going to give here, this is sufficient.
 * !!!
  * Otherwise do everything by location, no passing DTT to JDF_C_CODE
  + would that be enough?
 *
 */
static inline
void SET_UP_ARENA_DATATYPES(dplasma_data_collection_t * ddc, parsec_datatype_t parsec_type,
                            int uplo, int diag,
                            parsec_arena_datatype_t *default_adt, int *shape)
{
    parsec_tiled_matrix_dc_t * desc = ddc->dc_original;

    int top    = GET_ROWS(desc, 0);
    int bottom = GET_ROWS(desc, desc->m);
    int left   = GET_COLS(desc, 0);
    int right  = GET_COLS(desc, desc->n);
    int ld     = (desc->storage == matrix_Lapack ) ? desc->llm : desc->mb;

    int full_rows = GET_ROWS(desc, 1); //desc->mb > desc->m ? desc->m : desc->mb;
    int full_cols = GET_COLS(desc, 1); //desc->nb > desc->n ? desc->n : desc->nb;
  PARSEC_DEBUG_VERBOSE(23, parsec_debug_output,
    "SET_UP_ARENA_DATATYPES top %d bottom %d left %d right %d ld %d full %d x %d",
    top, bottom, left, right, ld, full_rows, full_cols);

    if( desc->storage != matrix_Lapack ){ /*TILED DESC*/
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DC00,    (*shape), LAPACK, "ADT_C00_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DC0N,    (*shape), LAPACK, "ADT_C0N_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DCM0,    (*shape), LAPACK, "ADT_CM0_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DCMN,    (*shape), LAPACK, "ADT_CMN_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DTOP,    (*shape), LAPACK, "ADT_TOP_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DBOTTOM, (*shape), LAPACK, "ADT_BOTTOM_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DLEFT,   (*shape), LAPACK, "ADT_LEFT_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DRIGHT,  (*shape), LAPACK, "ADT_RIGHT_ARENA");

        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DFULL_T, (*shape), TILED, "ADT_TILE_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DC00,    (*shape), TILED, "ADT_TILE_C00_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DC0N,    (*shape), TILED, "ADT_TILE_C0N_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DCM0,    (*shape), TILED, "ADT_TILE_CM0_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DCMN,    (*shape), TILED, "ADT_TILE_CMN_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DTOP,    (*shape), TILED, "ADT_TILE_TOP_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DBOTTOM, (*shape), TILED, "ADT_TILE_BOTTOM_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DLEFT,   (*shape), TILED, "ADT_TILE_LEFT_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, desc->mb, desc->nb, ld,       DRIGHT,  (*shape), TILED, "ADT_TILE_RIGHT_ARENA");

    } else { /*LAPACK DESC*/
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, top,       left,      ld,       DC00,    (*shape), LAPACK, "ADT_C00_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, top,       right,     ld,       DC0N,    (*shape), LAPACK, "ADT_C0N_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, bottom,    left,      ld,       DCM0,    (*shape), LAPACK, "ADT_CM0_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, bottom,    right,     ld,       DCMN,    (*shape), LAPACK, "ADT_CMN_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, top,       full_cols, ld,       DTOP,    (*shape), LAPACK, "ADT_TOP_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, bottom,    full_cols, ld,       DBOTTOM, (*shape), LAPACK, "ADT_BOTTOM_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, full_rows, left,      ld,       DLEFT,   (*shape), LAPACK, "ADT_LEFT_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, full_rows, right,     ld,       DRIGHT,  (*shape), LAPACK, "ADT_RIGHT_ARENA");

        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, full_rows, full_cols, full_rows, DFULL_T, (*shape), TILED, "ADT_TILE_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, top,       left,      top,       DC00,    (*shape), TILED, "ADT_TILE_C00_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, top,       right,     top,       DC0N,    (*shape), TILED, "ADT_TILE_C0N_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, bottom,    left,      bottom,    DCM0,    (*shape), TILED, "ADT_TILE_CM0_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, bottom,    right,     bottom,    DCMN,    (*shape), TILED, "ADT_TILE_CMN_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, top,       full_cols, top,       DTOP,    (*shape), TILED, "ADT_TILE_TOP_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, bottom,    full_cols, bottom,    DBOTTOM, (*shape), TILED, "ADT_TILE_BOTTOM_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, full_rows, left,      full_rows, DLEFT,   (*shape), TILED, "ADT_TILE_LEFT_ARENA");
        SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, full_rows, right,     full_rows, DRIGHT,  (*shape), TILED, "ADT_TILE_RIGHT_ARENA");
    }
    /* keep the default type on default_adt */
    SET_UP_ARENA_DTT(desc, parsec_type, uplo, diag, default_adt, full_rows, full_cols, ld, DFULL_T, (*shape), LAPACK, "ADT_DEFAULT_ARENA");
    (*shape)++;
}

static inline
void CLEAN_UP_ARENA_DATATYPES(const dplasma_data_collection_t * ddc, int max_shape)
{
  (void)ddc;
  /* attempt all different shapes for the ddc */
  parsec_arena_datatype_t adt;
  for(int shape = 0; shape < max_shape; shape++){
    for(int layout = 0; layout < MAX_LAYOUT; layout++){
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DFULL_T, shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DC00,    shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DC0N,    shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DCM0,    shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DCMN,    shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DTOP,    shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DBOTTOM, shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DLEFT,   shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
      if( dplasma_cleanup_datatype_info(ddc->dc_original, DRIGHT,  shape, layout, &adt) == 0) dplasma_matrix_del2arena( &adt);
    }
  }
}

#endif  /* DPLASMA_LAPACK_DATATYPE_H_HAS_BEEN_INCLUDED */

