/*
 * Copyright (c) 2020-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * $COPYRIGHT
 *
 */
#include "dplasma_arena_datatype.h"
#include <string.h>

#ifdef REUSE_ARENA_DATATYPE
static parsec_key_fn_t arena_key_fns = {
    .key_equal = parsec_hash_table_generic_64bits_key_equal,
    .key_print = parsec_hash_table_generic_64bits_key_print,
    .key_hash =  parsec_hash_table_generic_64bits_key_hash
};

static int dtt_key_equal(parsec_key_t keyA, parsec_key_t keyB, void *data)
{
    (void)data;
    return !strncmp((char*)keyA, (char*)keyB, DATATYPE_KEY_STR_SZ);
}

static char * dtt_key_print(char *buffer, size_t buffer_size,
                 parsec_key_t k, void *data)
{
  (void)data;
  char *key = (char*)k;
  snprintf(buffer, buffer_size, "<%s>", key);
  return buffer;
}

static uint64_t dtt_key_hash(parsec_key_t key, void *data)
{
    (void)data;
    char *str = (char*)key;
    uint64_t h = 1125899906842597ULL; // Large prime

    for (int i = 0; 0 != str[i]; i++) {
        h = 31*h + str[i];
    }
    return h;
}

static parsec_key_fn_t dtt_key_fns = {
    .key_equal = dtt_key_equal,
    .key_print = dtt_key_print,
    .key_hash = dtt_key_hash
};
#endif /* REUSE_ARENA_DATATYPE */

parsec_hash_table_t *dplasma_arenas = NULL;
parsec_hash_table_t *dplasma_datatypes = NULL;

static void dplasma_arenas_release(void *item, void*cb_data)
{
    (void)cb_data;
    dplasma_arena_t * arena_entry = (dplasma_arena_t *)item;
    assert(dplasma_arenas != NULL);
    parsec_hash_table_nolock_remove(dplasma_arenas, arena_entry->ht_item.key);
    PARSEC_OBJ_RELEASE(arena_entry->arena);
    free(arena_entry);
}

#ifdef REUSE_ARENA_DATATYPE
static void dplasma_datatypes_release(void *item, void*cb_data)
{
    (void)cb_data;
    dplasma_datatype_t * dtt_entry = (dplasma_datatype_t *)item;
    assert(dplasma_datatypes != NULL);
    parsec_hash_table_nolock_remove(dplasma_datatypes, dtt_entry->ht_item.key);
    parsec_type_free(&dtt_entry->dtt);
    free(((char*)dtt_entry->ht_item.key));
    free(dtt_entry);
}

static void dplasma_arenas_fini()
{
    if(dplasma_arenas != NULL){
        parsec_hash_table_for_all(dplasma_arenas, dplasma_arenas_release, NULL);
        parsec_hash_table_fini(dplasma_arenas);
        PARSEC_OBJ_RELEASE(dplasma_arenas);
    }
}

static void dplasma_datatypes_fini()
{
    if(dplasma_datatypes != NULL){
        parsec_hash_table_for_all(dplasma_datatypes, dplasma_datatypes_release, NULL);
        parsec_hash_table_fini(dplasma_datatypes);
        PARSEC_OBJ_RELEASE(dplasma_datatypes);
    }
}
#endif /* REUSE_ARENA_DATATYPE */

int dplasma_get_or_construct_arena(parsec_arena_t** parena,
                                   size_t elem_size,
                                   size_t alignment)
{
    int rc;
#ifdef REUSE_ARENA_DATATYPE
    /* Revisit if we are going to use different alignments, however,
     * dplasma uses always PARSEC_ARENA_ALIGNMENT_SSE with the
     * exception zhbrdt_wrapper which uses 16 (but still = PARSEC_ARENA_ALIGNMENT_SSE)
     */
    dplasma_arena_t * arena_entry = NULL;
    parsec_key_t k = (parsec_key_t)elem_size; //arena_key(elem_size, alignment);

    if( NULL == dplasma_arenas ) {
        parsec_hash_table_t *new_ht;
        new_ht = PARSEC_OBJ_NEW(parsec_hash_table_t);
        parsec_hash_table_init(new_ht,
                               offsetof(dplasma_arena_t, ht_item),
                               16, arena_key_fns, NULL);
        if(parsec_atomic_cas_ptr(&dplasma_arenas, NULL, new_ht)) {
            parsec_context_at_fini(dplasma_arenas_fini, NULL );
        } else {
            parsec_hash_table_fini(new_ht);
            PARSEC_OBJ_RELEASE(new_ht);
        }
    }

    parsec_hash_table_lock_bucket(dplasma_arenas, k);
    arena_entry = (dplasma_arena_t *)parsec_hash_table_nolock_find(dplasma_arenas, k);
    if( NULL != arena_entry ) {
        *parena = arena_entry->arena;
        PARSEC_OBJ_RETAIN(*parena);
        parsec_hash_table_unlock_bucket(dplasma_arenas, k);
        return 0;
    }
    /* New arena entry */
    *parena = PARSEC_OBJ_NEW(parsec_arena_t);
    rc = parsec_arena_construct(*parena, elem_size, alignment);

    if( 0 == rc ) {
        arena_entry = (dplasma_arena_t *)malloc(sizeof(dplasma_arena_t));
        arena_entry->ht_item.key = k;
        arena_entry->arena = *parena;
        parsec_hash_table_nolock_insert(dplasma_arenas, &arena_entry->ht_item);
        PARSEC_OBJ_RETAIN(*parena);
        parsec_hash_table_unlock_bucket(dplasma_arenas, k);
        return 0;
    }
    parsec_hash_table_unlock_bucket(dplasma_arenas, k);
    PARSEC_OBJ_RELEASE(*parena);
#else
    /* New arena */
    *parena = PARSEC_OBJ_NEW(parsec_arena_t);
    rc = parsec_arena_construct(*parena, elem_size, alignment);
#endif
    return rc;

}

void dplasma_cleanup_arenas(parsec_arena_t** arenas,
                           int count)
{
    int i;
    parsec_key_t k;
    dplasma_arena_t * arena_entry;
    if(dplasma_arenas == NULL)
        return;

    if(arenas == NULL){
        parsec_hash_table_for_all(dplasma_arenas, dplasma_arenas_release, NULL);
        return;
    }

    for(i = 0; i < count; i++){
        k = (parsec_key_t) arenas[i]->elem_size;
        arena_entry = (dplasma_arena_t *)parsec_hash_table_nolock_find(dplasma_arenas, k);
        if(arena_entry != NULL){
            dplasma_arenas_release( (void *)arena_entry, NULL);
        }
    }
}

int dplasma_get_or_construct_datatype(parsec_datatype_t *newtype, parsec_datatype_t oldtype,
                                      int uplo, int diag,
                                      unsigned int m, unsigned int n, unsigned int ld,
                                      int resized,
                                      ptrdiff_t * extent)
{
    int rc;
#ifdef REUSE_ARENA_DATATYPE
    dplasma_datatype_t * dtt_entry = NULL;
    char static_desc[DATATYPE_KEY_STR_SZ];
    snprintf( static_desc, DATATYPE_KEY_STR_SZ, "%"PRIxPTR"|%d|%d|%d|%d|%d|%d", (intptr_t)oldtype, uplo, diag, m, n, ld, resized);
    parsec_key_t k = (uint64_t)static_desc;

    if( NULL == dplasma_datatypes ) {
        parsec_hash_table_t *new_ht;
        new_ht = PARSEC_OBJ_NEW(parsec_hash_table_t);
        parsec_hash_table_init(new_ht,
                               offsetof(dplasma_datatype_t, ht_item),
                               16, dtt_key_fns, NULL);
        if(parsec_atomic_cas_ptr(&dplasma_datatypes, NULL, new_ht)) {
            parsec_context_at_fini(dplasma_datatypes_fini, NULL);
        } else {
            parsec_hash_table_fini(new_ht);
            PARSEC_OBJ_RELEASE(new_ht);
        }
    }

    parsec_hash_table_lock_bucket(dplasma_datatypes, k);
    dtt_entry = (dplasma_datatype_t *)parsec_hash_table_nolock_find(dplasma_datatypes, k);
    if( NULL != dtt_entry ) {
      *newtype = dtt_entry->dtt;
      *extent = dtt_entry->extent;
      parsec_hash_table_unlock_bucket(dplasma_datatypes, k);
      return 0;
    }
    /* New dtt entry */
    rc = parsec_matrix_define_datatype( newtype, oldtype,
                                        uplo, diag, m, n, ld,
                                        resized, extent);

    if( 0 == rc ) {
        dtt_entry = (dplasma_datatype_t *)malloc(sizeof(dplasma_datatype_t));
        dtt_entry->ht_item.key = (uint64_t)strdup(static_desc);
        dtt_entry->dtt = *newtype;
        dtt_entry->extent = *extent;
        parsec_hash_table_nolock_insert(dplasma_datatypes, &dtt_entry->ht_item);
    }
    parsec_hash_table_unlock_bucket(dplasma_datatypes, k);
#else
    /* New dtt */
    rc = parsec_matrix_define_datatype( newtype, oldtype,
                                        uplo, diag, m, n, ld,
                                        resized, extent);
#endif
    return rc;

}

int dplasma_cleanup_datatype( parsec_datatype_t *dtt )
{
#ifdef REUSE_ARENA_DATATYPE
    /* The key of a datatype can only be constructed with the arguments used during creation.
     * Thus, datatypes will be cleaned during callback during parsec context at fini.
     */
    (void)dtt;
    return 0;
#else
    return parsec_type_free( dtt );
#endif
}
