#ifndef DPLASMA_DATATYPE_H_HAS_BEEN_INCLUDED
#define DPLASMA_DATATYPE_H_HAS_BEEN_INCLUDED

/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */
#include "dplasma.h"

#include "utils/dplasma_arena_datatype.h"

/* Set up parsec_arena_datatype_t with the given specifications.
 * Both arena and datatype are reuse if one of the given specifications exits.
 */
static inline int
dplasma_get_or_construct_adt( parsec_arena_datatype_t *adt, parsec_datatype_t oldtype, size_t alignment,
                              int uplo, int diag,
                              unsigned int mb, unsigned int nb, unsigned int ld,
                              int resized)
{
    ptrdiff_t extent;
    int rc;
    rc = dplasma_get_or_construct_datatype( &(adt->opaque_dtt), oldtype,
                                            uplo, diag, mb, nb, ld, resized, &extent);
    assert(rc == 0);
    if(0 != rc) return rc;
    rc = dplasma_get_or_construct_arena( &(adt->arena), extent, alignment);
    assert(rc == 0);
    return rc;
}

static inline int
dplasma_add2arena_rectangle( parsec_arena_datatype_t *adt, size_t elem_size, size_t alignment,
                             parsec_datatype_t oldtype,
                             unsigned int tile_mb, unsigned int tile_nb, int resized )
{
    (void)elem_size;
    return dplasma_get_or_construct_adt( adt, oldtype, alignment,
                                         matrix_UpperLower, 1, tile_mb, tile_nb, tile_mb, resized);
}

static inline int
dplasma_add2arena_tile( parsec_arena_datatype_t *adt, size_t elem_size, size_t alignment,
                        parsec_datatype_t oldtype,
                        unsigned int tile_mb )
{
    (void)elem_size;
    return dplasma_get_or_construct_adt( adt, oldtype, alignment,
                                         matrix_UpperLower, 1, tile_mb, tile_mb, tile_mb, -1);
}

static inline int
dplasma_add2arena_upper( parsec_arena_datatype_t *adt, size_t elem_size, size_t alignment,
                         parsec_datatype_t oldtype,
                         unsigned int tile_mb, int diag )
{
    (void)elem_size;
    return dplasma_get_or_construct_adt( adt, oldtype, alignment,
                                         matrix_Upper, diag, tile_mb, tile_mb, tile_mb, -1);
}

static inline int
dplasma_add2arena_lower( parsec_arena_datatype_t *adt, size_t elem_size, size_t alignment,
                         parsec_datatype_t oldtype,
                         unsigned int tile_mb, int diag )
{
    (void)elem_size;
    return dplasma_get_or_construct_adt( adt, oldtype, alignment,
                                         matrix_Lower, diag, tile_mb, tile_mb, tile_mb, -1);
}

static inline int
dplasma_matrix_del2arena( parsec_arena_datatype_t *adt )
{
    return dplasma_cleanup_datatype( &adt->opaque_dtt );
}

#endif  /* DPLASMA_DATATYPE_H_HAS_BEEN_INCLUDED */

