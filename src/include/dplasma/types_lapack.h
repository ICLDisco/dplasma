#ifndef DPLASMA_LAPACK_DATATYPE_H_HAS_BEEN_INCLUDED
#define DPLASMA_LAPACK_DATATYPE_H_HAS_BEEN_INCLUDED

/*
 * Copyright (c) 2010-2023 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */
#include "dplasma.h"

#include "types.h"
#include "utils/dplasma_lapack_adtt.h"

/* Support for TILED/LAPACK matrix with non homogeneous datatypes across tiles.
 * NOTE: we are operating with the following condition:
 * For a given data collection, if we reuse the datatype for one shape on different
 * locations, then other shapes will also be reusing a datatype for those locations.
 * The same happens between layouts.
 *
 * E.g.:
 *    if LAPACK ( DFULL_T | DC00 | DTOP) --> dtt1
 *    then TILE ( DFULL_T | DC00 | DTOP) --> dtt2  & dtt2 is a correct equivalence to dtt1
 *
 * As adtt_info stores key(dtt) -> first info that was set up, when translating that dtt
 * to the equivalent on a different shape/layout, we reuse the info location.
 *
 * This is sufficient for LAPACK/TILED support.
 *
 */

/* Obtain number of elements rows in tile in matrix tile row m. */
static inline int dplasma_get_rows(const parsec_tiled_matrix_t *dM, int m)
{
  int nrows = dM->mb;
  if(m == 0)     nrows =   ((dM->i % dM->mb) !=0          ? (dM->mb - (dM->i % dM->mb)) : dM->mb);
  if(m == dM->m) nrows =  (((dM->i + dM->m) % dM->mb) !=0 ? ((dM->i + dM->m) % dM->mb)  : dM->mb);
  return (nrows > dM->m) ? dM->m : nrows;
}

/* Obtain number of elements cols in tile in matrix tile col n. */
static inline int dplasma_get_cols(const parsec_tiled_matrix_t *dM, int n)
{
  int ncols = dM->nb;
  if(n == 0)     ncols =  ((dM->j % dM->nb) !=0          ? (dM->nb - (dM->j % dM->nb)) : dM->nb);
  if(n == dM->n) ncols = (((dM->j + dM->n) % dM->nb) !=0 ? ((dM->j + dM->n) % dM->nb)  : dM->nb);
  return (ncols > dM->n) ? dM->n : ncols;
}

/**
 * Setup info for the datatype if we can reuse it, obtain arena and datatype and setup info.
 */
static inline
void dplasma_setup_adtt_loc(dplasma_data_collection_t * ddc,
                            parsec_datatype_t parsec_type,
                            int uplo, int diag,
                            int rows, int cols, int ld,
                            int loc, int shape, int layout,
                            parsec_arena_datatype_t *adt_default,
                            int full_rows, int full_cols)
{
    lapack_info_t info;
    parsec_arena_datatype_t adt;

    info.lda = ld;
    info.rows = rows;
    info.cols = cols;
    info.loc = loc;
    info.shape = shape;
    info.layout = layout;

    if ( (adt_default != NULL) && (full_rows == rows) && (full_cols == cols) ) {
        /* Reuse external datatype */
        /* Reuse arena through dplasma reuse, otherwise we don't know how to release
         * during destruct.
         */
        PARSEC_OBJ_RETAIN(adt_default->arena);
        dplasma_set_datatype_info(ddc, *adt_default, &info);
    } else {
        /* Reuse dplasma arena & datatype */
        dplasma_get_or_construct_adt(&adt, parsec_type,
                                     PARSEC_ARENA_ALIGNMENT_SSE,
                                     uplo, diag,
                                     rows, cols, ld,
                                     -1/* resized = -1 for all dplasma types */);
        dplasma_set_datatype_info(ddc, adt, &info);
    }
}

/**
 * Setup parsec_arena_datatype_t for all locations for TILED/LAPACK matrix for a shape.
 * Shape represents an extraction (e.g. rectangle, lower triangle, upper triangle)
 * for a GIVEN DATACOLLECTION:
 * - on lapack different matrix can have different shapes (e.g. != LDA)
 * - always different extractions of data imply different shapes (UPPER, LOWER, etc).
 * Workaround to deal with external datatype provided on the data collection.
 * Assuming datatype from the data collection corresponds with the extraction of a full
 * tile (no upper, lower) and forcing its reuse on that case by checking if uplo & diag.
 * The other datatypes created by this call use the dplasma reuse of datatypes.
 */
static inline
void dplasma_setup_adtt_all_loc(dplasma_data_collection_t * ddc, parsec_datatype_t parsec_type,
                                int uplo, int diag, int *shape)
{
    parsec_tiled_matrix_t * desc = ddc->dc_original;
    parsec_arena_datatype_t *adt_default = NULL;
    parsec_arena_datatype_t adt;

    int top    = dplasma_get_rows(desc, 0);
    int bottom = dplasma_get_rows(desc, desc->m);
    int left   = dplasma_get_cols(desc, 0);
    int right  = dplasma_get_cols(desc, desc->n);
    int ld     = (desc->storage == PARSEC_MATRIX_LAPACK ) ? desc->llm : desc->mb;

    int full_rows = desc->mb;
    int full_cols = desc->nb;

    if( uplo == PARSEC_MATRIX_FULL ){
        /* Generating datatypes for full tiles. No other extraction of data.
         * Reuse data collection dtt.
         */
        ptrdiff_t lb = 0;
        ptrdiff_t extent;
        adt.opaque_dtt = ((parsec_data_collection_t *)desc)->default_dtt;
        parsec_type_extent(adt.opaque_dtt, &lb, &extent);
        dplasma_get_or_construct_arena(&adt.arena, extent, PARSEC_ARENA_ALIGNMENT_SSE);
        adt_default = &adt;
    }


    if( desc->storage != PARSEC_MATRIX_LAPACK ){ /*TILED DESC*/
        /* on tile format, data collection dtt is the entire MBxNB even if that memory is not logically on the matrix.
         * force usage of the data collection dtt.
         */
        /* Try to reuse adt_default when not scalapack storage */
        for(int loc = 0; loc < LOC_SZ; loc++){
            dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, desc->mb, desc->nb, ld, loc, (*shape), TILED,  adt_default, full_rows, full_cols);
            dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, desc->mb, desc->nb, ld, loc, (*shape), LAPACK, adt_default, full_rows, full_cols);
        }
    } else { /*LAPACK DESC*/
        /* on lapack format, dtt is the may never be MBxNB (if not enough rows).
         */

        int full_rows = dplasma_get_rows(desc, 1);
        int full_cols = dplasma_get_cols(desc, 1);

        /* Try to reuse adt_default for LAPACK layout */
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, full_rows, full_cols, ld,       DFULL_T,  (*shape), LAPACK, adt_default, full_rows, full_cols);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, top,       left,      ld,       DC00,     (*shape), LAPACK, adt_default, full_rows, full_cols);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, top,       right,     ld,       DC0N,     (*shape), LAPACK, adt_default, full_rows, full_cols);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, bottom,    left,      ld,       DCM0,     (*shape), LAPACK, adt_default, full_rows, full_cols);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, bottom,    right,     ld,       DCMN,     (*shape), LAPACK, adt_default, full_rows, full_cols);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, top,       full_cols, ld,       DTOP,     (*shape), LAPACK, adt_default, full_rows, full_cols);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, bottom,    full_cols, ld,       DBOTTOM,  (*shape), LAPACK, adt_default, full_rows, full_cols);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, full_rows, left,      ld,       DLEFT,    (*shape), LAPACK, adt_default, full_rows, full_cols);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, full_rows, right,     ld,       DRIGHT,   (*shape), LAPACK, adt_default, full_rows, full_cols);

        /* Always create adt for TILED layout on scalapack storage */
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, full_rows, full_cols, full_rows, DFULL_T, (*shape), TILED, NULL, -1, -1);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, top,       left,      top,       DC00,    (*shape), TILED, NULL, -1, -1);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, top,       right,     top,       DC0N,    (*shape), TILED, NULL, -1, -1);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, bottom,    left,      bottom,    DCM0,    (*shape), TILED, NULL, -1, -1);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, bottom,    right,     bottom,    DCMN,    (*shape), TILED, NULL, -1, -1);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, top,       full_cols, top,       DTOP,    (*shape), TILED, NULL, -1, -1);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, bottom,    full_cols, bottom,    DBOTTOM, (*shape), TILED, NULL, -1, -1);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, full_rows, left,      full_rows, DLEFT,   (*shape), TILED, NULL, -1, -1);
        dplasma_setup_adtt_loc(ddc, parsec_type, uplo, diag, full_rows, right,     full_rows, DRIGHT,  (*shape), TILED, NULL, -1, -1);
    }

    PARSEC_DEBUG_VERBOSE(22, parsec_debug_output,
    "dplasma_setup_adtt_all_loc top %d bottom %d left %d right %d ld %d full %d x %d",
    top, bottom, left, right, ld, full_rows, full_cols);

    (*shape)++;
}


/**
 * Clean up parsec_arena_datatype_t for all locations for TILED/LAPACK matrix.
 * generates base tile types depending on uplo/diag values
 */
static inline
void dplasma_clean_adtt_all_loc(const dplasma_data_collection_t * ddc, int max_shape)
{
    /* attempt all different shapes for the ddc */
    parsec_arena_datatype_t adt;
    lapack_info_t info;
    for(int shape = 0; shape < max_shape; shape++){
        for(int layout = 0; layout < MAX_LAYOUT; layout++){
            for(int loc = 0; loc < LOC_SZ; loc++){
                info.loc = loc;
                info.shape = shape;
                info.layout = layout;
                if ( dplasma_cleanup_datatype_info(ddc, &info, &adt) == 0) {
                    /* Retained when reusing it, not set on JDF array, not release by taskpool destructor */
                    PARSEC_OBJ_RELEASE(adt.arena);
                    /* A datatype being registered with multiple info we should
                     * only release the datatype once all info have been
                     * removed. We can track this using the refcount on the
                     * arena itself.
                     */
                    if( NULL == adt.arena )
                        dplasma_matrix_del2arena(&adt);
                }
            }
        }
    }

}

#endif  /* DPLASMA_LAPACK_DATATYPE_H_HAS_BEEN_INCLUDED */
