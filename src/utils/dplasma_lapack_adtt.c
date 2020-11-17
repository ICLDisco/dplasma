#include "dplasma_lapack_adtt.h"
#include <string.h>

/*************************************************************/
/* Management of LAPACK-TILED datatypes info and translation */
/*************************************************************/

static int dtt_key_equal(parsec_key_t keyA, parsec_key_t keyB, void *data)
{
    (void)data;
    return !strncmp((char*)keyA, (char*)keyB, LAPACK_ADT_KEY_STR_SZ);
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

static parsec_key_fn_t dtt_lapack_key_fns = {
    .key_equal = dtt_key_equal,
    .key_print = dtt_key_print,
    .key_hash = dtt_key_hash
};

parsec_hash_table_t *dplasma_datatypes_lapack_helper = NULL;

static void dplasma_datatypes_lapack_helper_release(void *item, void*cb_data)
{
    (void)cb_data;
    dplasma_datatype_lapack_helper_t * dtt_entry = (dplasma_datatype_lapack_helper_t *)item;
    assert(dplasma_datatypes_lapack_helper != NULL);
    parsec_hash_table_nolock_remove(dplasma_datatypes_lapack_helper, dtt_entry->ht_item.key);
    free(dtt_entry);
}

static void dplasma_datatypes_info_fini()
{
    if(dplasma_datatypes_lapack_helper != NULL){
        parsec_hash_table_for_all(dplasma_datatypes_lapack_helper, dplasma_datatypes_lapack_helper_release, NULL);
        parsec_hash_table_fini(dplasma_datatypes_lapack_helper);
        PARSEC_OBJ_RELEASE(dplasma_datatypes_lapack_helper);
    }
}

static int dplasma_set_dtt_to_info(const parsec_tiled_matrix_dc_t *dc, parsec_arena_datatype_t adt, int lda, int rows, int cols, int loc, int shape, int layout)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;

    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|%p|-1|-1|-1", dc, adt.opaque_dtt);
    parsec_key_t k = (uint64_t)static_desc;
    if( NULL == dplasma_datatypes_lapack_helper ) {
        dplasma_datatypes_lapack_helper = PARSEC_OBJ_NEW(parsec_hash_table_t);
        parsec_hash_table_init(dplasma_datatypes_lapack_helper,
                               offsetof(dplasma_datatype_lapack_helper_t, ht_item),
                               16, dtt_lapack_key_fns, NULL);
        parsec_context_at_fini(dplasma_datatypes_info_fini, NULL);
    } else {
        dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
        if( NULL != dtt_entry ) {
            assert(dtt_entry->lda == lda);
            return 0;
        }
    }

    dtt_entry = (dplasma_datatype_lapack_helper_t *)malloc(sizeof(dplasma_datatype_lapack_helper_t));
    dtt_entry->ht_item.key = (uint64_t)strdup(static_desc);
    dtt_entry->lda = lda;
    dtt_entry->rows = rows;
    dtt_entry->cols = cols;
    dtt_entry->loc = loc;
    dtt_entry->shape = shape;
    dtt_entry->layout = layout;
    dtt_entry->adt = adt;
    parsec_hash_table_nolock_insert(dplasma_datatypes_lapack_helper, &dtt_entry->ht_item);
    PARSEC_DEBUG_VERBOSE(27, parsec_debug_output, " insert dtt -> info lda %d rows %d cols %d loc %d shape %d dtt %p (%30s)", lda, rows, cols, loc, shape, adt.opaque_dtt, static_desc);
    return 0;
}

static int dplasma_set_info_to_dtt(const parsec_tiled_matrix_dc_t *dc, parsec_arena_datatype_t adt, int lda, int rows, int cols, int loc, int shape, int layout)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;

    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|-1|%d|%d|%d", dc, loc, shape, layout);
    parsec_key_t k = (uint64_t)static_desc;
    if( NULL == dplasma_datatypes_lapack_helper ) {
        dplasma_datatypes_lapack_helper = PARSEC_OBJ_NEW(parsec_hash_table_t);
        parsec_hash_table_init(dplasma_datatypes_lapack_helper,
                               offsetof(dplasma_datatype_lapack_helper_t, ht_item),
                               16, dtt_lapack_key_fns, NULL);
        parsec_context_at_fini(dplasma_datatypes_info_fini, NULL);
    } else {
        dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
        if( NULL != dtt_entry ) {
            assert(dtt_entry->adt.opaque_dtt == adt.opaque_dtt);
            return 0;
        }
    }

    dtt_entry = (dplasma_datatype_lapack_helper_t *)malloc(sizeof(dplasma_datatype_lapack_helper_t));
    dtt_entry->ht_item.key = (uint64_t)strdup(static_desc);
    dtt_entry->lda = lda;
    dtt_entry->rows = rows;
    dtt_entry->cols = cols;
    dtt_entry->loc = loc;
    dtt_entry->shape = shape;
    dtt_entry->layout = layout;
    dtt_entry->adt = adt;
    parsec_hash_table_nolock_insert(dplasma_datatypes_lapack_helper, &dtt_entry->ht_item);
    PARSEC_DEBUG_VERBOSE(27, parsec_debug_output, " insert dtt <- info lda %d rows %d cols %d loc %d shape %d dtt %p (%30s)", lda, rows, cols, loc, shape, adt.opaque_dtt, static_desc);
    return 0;
}

int dplasma_set_datatype_info(const parsec_tiled_matrix_dc_t *dc, parsec_arena_datatype_t adt, int lda, int rows, int cols, int loc, int shape, int layout)
{
    dplasma_set_dtt_to_info(dc, adt, lda, rows, cols, loc, shape, layout);
    dplasma_set_info_to_dtt(dc, adt, lda, rows, cols, loc, shape, layout);
    return 0;
}

int dplasma_get_datatype_info(const parsec_tiled_matrix_dc_t *dc, parsec_datatype_t dtt, int *lda, int *rows, int *cols, int *loc, int *shape, int *layout, parsec_arena_datatype_t **adt)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;
    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|%p|-1|-1|-1", dc, dtt);
    parsec_key_t k = (uint64_t)static_desc;
    assert(dplasma_datatypes_lapack_helper!= NULL);
    dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
    assert(dtt_entry != NULL);
    *lda = dtt_entry->lda;
    *rows = dtt_entry->rows;
    *cols = dtt_entry->cols;
    *loc = dtt_entry->loc;
    *shape = dtt_entry->shape;
    *layout = dtt_entry->layout;
    *adt  = &dtt_entry->adt;
    return 0;
}

int dplasma_get_datatype(const parsec_tiled_matrix_dc_t *dc, int loc, int shape, int layout, int *lda, parsec_arena_datatype_t **adt)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;
    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|-1|%d|%d|%d", dc, loc, shape, layout);
    parsec_key_t k = (uint64_t)static_desc;
    assert(dplasma_datatypes_lapack_helper!= NULL);
    dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
    assert(dtt_entry != NULL);
    *lda = dtt_entry->lda;
    *adt  = &dtt_entry->adt;
    return 0;
}

int dplasma_cleanup_datatype_info(const parsec_tiled_matrix_dc_t *dc, int loc, int shape, int layout, parsec_arena_datatype_t *adt)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;
    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    parsec_key_t k;
    assert(dplasma_datatypes_lapack_helper!= NULL);
    /* Remove the two entries from the hash table */
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|-1|%d|%d|%d", dc, loc, shape, layout);
    k = (uint64_t)static_desc;
    dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
    if(dtt_entry != NULL) {//assert(dtt_entry != NULL);
        *adt = dtt_entry->adt;
        dplasma_datatypes_lapack_helper_release((void*)dtt_entry, NULL);

        snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|%p|-1|-1|-1", dc, adt->opaque_dtt);
        k = (uint64_t)static_desc;
        dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
        if(dtt_entry != NULL){ /* same desc + dtt can be reuse for several loc & shape */
            dplasma_datatypes_lapack_helper_release((void*)dtt_entry, NULL);
        }
        return 0;
    }
    /* No found */
    return -1;
}


/***************************************************************************/
/* Management of dplasma_data_collection_t to wrap parsec_datacollecctions */
/***************************************************************************/

static parsec_data_t* data_of(parsec_data_collection_t *desc, ...)
{
    int m, n;
    va_list ap;
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;

    /* Get coordinates */
    va_start(ap, desc);
    m = (int)va_arg(ap, unsigned int);
    n = (int)va_arg(ap, unsigned int);
    va_end(ap);
    parsec_data_t* dt = ((parsec_data_collection_t*)ddc->dc_original)->data_of(((parsec_data_collection_t *)ddc->dc_original), m, n);
    parsec_data_copy_t *cp = dt->device_copies[0];

    /* Correct datatype if not default because of location */
    int loc = LOC_TYPE( ddc->dc_original, m, n);
    if(loc != DFULL_T){
        int lda, rows, cols, shape, old_loc, layout, rc;
        parsec_arena_datatype_t *adt;
        /* obtain the shape */
        rc = dplasma_get_datatype_info(ddc->dc_original, cp->dtt, &lda, &rows, &cols, &old_loc, &shape, &layout, &adt);
        assert(rc == 0);
        /* obtain appropriate dtt for that location, shape and layout */
        rc = dplasma_get_datatype(ddc->dc_original, loc, shape, layout, &lda, &adt);
        assert(rc == 0);
        PARSEC_DEBUG_VERBOSE(27, parsec_debug_output, "data_of CP %p [old type %p] loc %d -> dtt %p target_shape %d layout %d", cp, cp->dtt, loc, adt->opaque_dtt, shape, layout);
        parsec_set_data_copy_datatype(cp, adt->opaque_dtt);
    }
    return dt;
}

static parsec_data_t* data_of_key(parsec_data_collection_t *desc, parsec_data_key_t key)
{
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    return ((parsec_data_collection_t*)ddc->dc_original)->data_of_key(((parsec_data_collection_t *)ddc->dc_original), key);
}

static uint32_t rank_of(parsec_data_collection_t *desc, ...)
{
    uint32_t ret;
    va_list ap;
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    va_start(ap, desc);
    ret = ((parsec_data_collection_t*)ddc->dc_original)->rank_of(((parsec_data_collection_t *)ddc->dc_original), ap);
    va_end(ap);
    return ret;
}

static uint32_t rank_of_key(parsec_data_collection_t *desc, parsec_data_key_t key)
{
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    return ((parsec_data_collection_t*)ddc->dc_original)->rank_of_key(((parsec_data_collection_t *)ddc->dc_original), key);
}


static int32_t vpid_of(parsec_data_collection_t *desc, ...)
{
    int32_t ret;
    va_list ap;
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    va_start(ap, desc);
    ret = ((parsec_data_collection_t*)ddc->dc_original)->vpid_of(((parsec_data_collection_t *)ddc->dc_original), ap);
    va_end(ap);
    return ret;
}

static int32_t vpid_of_key(parsec_data_collection_t *desc, parsec_data_key_t key)
{
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    return ((parsec_data_collection_t*)ddc->dc_original)->vpid_of_key(((parsec_data_collection_t *)ddc->dc_original), key);
}

static int register_memory(parsec_data_collection_t *desc, parsec_device_module_t *device)
{
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    return ((parsec_data_collection_t*)ddc->dc_original)->register_memory(((parsec_data_collection_t *)ddc->dc_original), device);
}

static int unregister_memory(parsec_data_collection_t *desc, parsec_device_module_t *device)
{
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    return ((parsec_data_collection_t*)ddc->dc_original)->unregister_memory(((parsec_data_collection_t *)ddc->dc_original), device);
}

static int key_to_string(parsec_data_collection_t * desc, parsec_data_key_t datakey, char * buffer, uint32_t buffer_size)
{
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    return ((parsec_data_collection_t*)ddc->dc_original)->key_to_string(((parsec_data_collection_t *)ddc->dc_original), datakey, buffer, buffer_size);
}

dplasma_data_collection_t *dplasma_wrap_data_collection(parsec_tiled_matrix_dc_t *A)
{
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)malloc(sizeof(dplasma_data_collection_t));
    parsec_data_collection_init((parsec_data_collection_t *)ddc,
                                 A->super.nodes,
                                 A->super.myrank );

#if defined(PARSEC_PROF_TRACE)
    asprintf(&((parsec_data_collection_t *)ddc)->key_dim, "%s_dplasmawrapped", ((parsec_data_collection_t *)A)->key_dim);
#endif
    parsec_data_collection_set_key((parsec_data_collection_t*)ddc, ((parsec_data_collection_t *)A)->key_base);

    ddc->dc_original = A;
    parsec_data_collection_t * dc   = (parsec_data_collection_t *)ddc;
    dc->rank_of      = rank_of;
    dc->vpid_of      = vpid_of;
    dc->data_of      = data_of;
    dc->rank_of_key  = rank_of_key;
    dc->vpid_of_key  = vpid_of_key;
    dc->data_of_key  = data_of_key;
    dc->register_memory   = register_memory;
    dc->unregister_memory = unregister_memory;
    dc->key_to_string     = key_to_string;

    ((parsec_data_collection_t *)ddc)->default_dtt = ((parsec_data_collection_t *)A)->default_dtt;

    return ddc;
}

void dplasma_unwrap_data_collection(dplasma_data_collection_t *ddc)
{
    free(ddc);
}

