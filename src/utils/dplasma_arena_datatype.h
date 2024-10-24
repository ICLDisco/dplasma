#ifndef DPLASMA_ARENA_DATATYPE_H_HAS_BEEN_INCLUDED
#define DPLASMA_ARENA_DATATYPE_H_HAS_BEEN_INCLUDED

/*
 * Copyright (c) 2010-2023 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */

#include "dplasma.h"
#include "parsec/arena.h"
#include "parsec/class/parsec_hash_table.h"

#undef REUSE_ARENA_DATATYPE
#define DATATYPE_KEY_STR_SZ 100

extern parsec_hash_table_t *dplasma_arenas;

struct dplasma_arena_s {
    parsec_hash_table_item_t ht_item;
    parsec_arena_t * arena;
};
typedef struct dplasma_arena_s dplasma_arena_t;

extern parsec_hash_table_t *dplasma_datatypes;

struct dplasma_datatype_s {
    parsec_hash_table_item_t ht_item;
    parsec_datatype_t dtt;
    ptrdiff_t extent;
};
typedef struct dplasma_datatype_s dplasma_datatype_t;

/* Function to obtain an arena with the given specifications.
 * If the hashmap doesn't contain one equivalent, a new arena is created.
 * The returned arena is retain to avoid it to be freed during the destructor
 * of the taskpool.
 * Arenas will be freed during parsec_fini (parsec external at fini callback) or
 * when invoking arenas_cleanup.
 * The user is not allowed to free directly the arenas.
 */
int dplasma_get_or_construct_arena(parsec_arena_t** parena,
                                   size_t elem_size,
                                   size_t alignment);

/* Function to clean up the arenas in the list that are tracked in the hashmap.
 * Arenas are released and will be destructed when any taskpool using them has been
 * destruct.
 * arenas: list of arenas to be cleaned up. If NULL all arenas are removed.
 * count: number of arenas in the list.
 */
void dplasma_cleanup_arenas(parsec_arena_t** arenas,
                            int count);


/* Function to obtain an datatatype with the given specifications.
 * If the hashmap doesn't contain one equivalent, a new datatype is created.
 * Datatypes will be freed during parsec_fini (parsec external at fini callback) or
 * when invoking datatypes_cleanup.
 * The user is not allowed to free directly the datatypes.
 */
int dplasma_get_or_construct_datatype(parsec_datatype_t *newtype, parsec_datatype_t oldtype,
                                      int uplo, int diag,
                                      unsigned int m, unsigned int n, unsigned int ld,
                                      int resized,
                                      ptrdiff_t * extent);

/* Function to clean up ONE datatype.
 * When reusing arenas and datatypes the call has no effect. The key of a datatype can only
 * be constructed with the arguments used during creation. Thus, datatypes will be cleaned during
 * callback during parsec context at fini.
 * When reusing is not activated, the datatype will be freed.
 */
int dplasma_cleanup_datatype( parsec_datatype_t *adt );

#endif  /* DPLASMA_ARENA_DATATYPE_H_HAS_BEEN_INCLUDED */
