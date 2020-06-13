#ifndef DPLASMA_DATATYPE_H_HAS_BEEN_INCLUDED
#define DPLASMA_DATATYPE_H_HAS_BEEN_INCLUDED

/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */
#include "dplasma.h"

#include "utils/dplasma_arena_datatype.h"

static inline int
dplasma_add2arena_rectangle( parsec_arena_datatype_t *adt, size_t elem_size, size_t alignment,
                             parsec_datatype_t oldtype,
                             unsigned int tile_mb, unsigned int tile_nb, int resized )
{
    (void)elem_size;
    ptrdiff_t extent; 
    int rc;
    rc = dplasma_get_or_construct_datatype( &adt->opaque_dtt, oldtype,
                                            matrix_UpperLower, 1, tile_mb, tile_nb, tile_mb,
                                            resized, &extent);
    if(0 != rc) return rc;
    rc = dplasma_get_or_construct_arena(&adt->arena, extent, alignment);
    return rc;
}

static inline int
dplasma_add2arena_tile( parsec_arena_datatype_t *adt, size_t elem_size, size_t alignment,
                        parsec_datatype_t oldtype, unsigned int tile_mb )
{
    (void)elem_size;
    ptrdiff_t extent; 
    int rc;
    rc = dplasma_get_or_construct_datatype( &adt->opaque_dtt, oldtype,
                                            matrix_UpperLower, 1, tile_mb, tile_mb, tile_mb,
                                            -1, &extent);
    if(0 != rc) return rc;
    rc = dplasma_get_or_construct_arena(&adt->arena, extent, alignment);
    return rc;
}

static inline int
dplasma_add2arena_upper( parsec_arena_datatype_t *adt, size_t elem_size, size_t alignment,
                         parsec_datatype_t oldtype, unsigned int tile_mb, int diag )
{
    (void)elem_size;
    ptrdiff_t extent; 
    int rc;
    rc = dplasma_get_or_construct_datatype( &adt->opaque_dtt, oldtype,
                                            matrix_Upper, diag, tile_mb, tile_mb, tile_mb,
                                            -1, &extent);
    if(0 != rc) return rc;
    rc = dplasma_get_or_construct_arena(&adt->arena, extent, alignment);
    return rc;
}

static inline int
dplasma_add2arena_lower( parsec_arena_datatype_t *adt, size_t elem_size, size_t alignment,
                         parsec_datatype_t oldtype, unsigned int tile_mb, int diag )
{
    (void)elem_size;
    ptrdiff_t extent; 
    int rc;
    rc = dplasma_get_or_construct_datatype( &adt->opaque_dtt, oldtype,
                                            matrix_Lower, diag, tile_mb, tile_mb, tile_mb,
                                            -1, &extent);
    if(0 != rc) return rc;
    rc = dplasma_get_or_construct_arena(&adt->arena, extent, alignment);
    return rc;
}

static inline int
dplasma_matrix_del2arena( parsec_arena_datatype_t *adt )
{
    return dplasma_cleanup_datatype( &adt->opaque_dtt );
}

#endif  /* DPLASMA_DATATYPE_H_HAS_BEEN_INCLUDED */
