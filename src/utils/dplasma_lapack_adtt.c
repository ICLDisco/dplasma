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

/*************************************************************/
/* Static helper functions that are not protected.           */
/*************************************************************/

static void dplasma_datatypes_lapack_helper_release(void *item, void*cb_data)
{
    (void)cb_data;
    dplasma_datatype_lapack_helper_t * dtt_entry = (dplasma_datatype_lapack_helper_t *)item;
    assert(dplasma_datatypes_lapack_helper != NULL);
    parsec_hash_table_nolock_remove(dplasma_datatypes_lapack_helper, dtt_entry->ht_item.key);
    free(((char*)dtt_entry->ht_item.key));
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

/*************************************************************/
/* Static helper functions that are lock protected.          */
/*************************************************************/

static inline int
dplasma_set_helper_hash_table_atomic(void)
{
    if( NULL == dplasma_datatypes_lapack_helper ) {
        parsec_hash_table_t* new_ht = PARSEC_OBJ_NEW(parsec_hash_table_t);
        parsec_hash_table_init(new_ht,
                               offsetof(dplasma_datatype_lapack_helper_t, ht_item),
                               16, dtt_lapack_key_fns, NULL);
        if(parsec_atomic_cas_ptr(&dplasma_datatypes_lapack_helper, NULL, new_ht)) {
            parsec_context_at_fini(dplasma_datatypes_info_fini, NULL);
            return 1;
        }
        parsec_hash_table_fini(new_ht);
        PARSEC_OBJ_RELEASE(new_ht);
    }
    return 0;
}

static int
dplasma_set_dtt_to_info( const dplasma_data_collection_t *dc,
                         parsec_arena_datatype_t adt,
                         const lapack_info_t* info)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;

    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|%"PRIxPTR"|-1|-1|-1", dc, (intptr_t)adt.opaque_dtt);
    parsec_key_t k = (uint64_t)static_desc;

    (void)dplasma_set_helper_hash_table_atomic();

    parsec_hash_table_lock_bucket(dplasma_datatypes_lapack_helper, k);  /* protect access to the hash table */
    /* Find the entry */
    dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
    if( NULL != dtt_entry ) {
        assert(dtt_entry->info.lda == info->lda);
        parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);
        return 0;
    }

    dtt_entry = (dplasma_datatype_lapack_helper_t *)malloc(sizeof(dplasma_datatype_lapack_helper_t));
    dtt_entry->ht_item.key = (uint64_t)strdup(static_desc);
    dtt_entry->info = *info;
    dtt_entry->adt = adt;
    parsec_hash_table_nolock_insert(dplasma_datatypes_lapack_helper, &dtt_entry->ht_item);
    parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);
    PARSEC_DEBUG_VERBOSE(27, parsec_debug_output,
        " insert dtt -> info lda %d rows %d cols %d loc %d shape %d layout %d dtt %p arena %p (%30s)",
        info->lda, info->rows, info->cols, info->loc, info->shape, info->layout, adt.opaque_dtt, adt.arena, static_desc);
    return 0;
}

static int
dplasma_set_info_to_dtt( const dplasma_data_collection_t *dc,
                         parsec_arena_datatype_t adt,
                         const lapack_info_t* info)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;

    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|-1|%d|%d|%d", dc, info->loc, info->shape, info->layout);
    parsec_key_t k = (uint64_t)static_desc;

    (void)dplasma_set_helper_hash_table_atomic();

    parsec_hash_table_lock_bucket(dplasma_datatypes_lapack_helper, k);  /* protect access to the hash table */
    /* Find the entry */
    dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
    if( NULL != dtt_entry ) {
        assert(dtt_entry->adt.opaque_dtt == adt.opaque_dtt);
        parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);
        return 0;
    }

    dtt_entry = (dplasma_datatype_lapack_helper_t *)malloc(sizeof(dplasma_datatype_lapack_helper_t));
    dtt_entry->ht_item.key = (uint64_t)strdup(static_desc);
    dtt_entry->info = *info;
    dtt_entry->adt = adt;
    parsec_hash_table_nolock_insert(dplasma_datatypes_lapack_helper, &dtt_entry->ht_item);
    parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);
    PARSEC_DEBUG_VERBOSE(27, parsec_debug_output,
        " insert dtt <- info lda %d rows %d cols %d loc %d shape %d layout %d dtt %p arena %p (%30s)",
        info->lda, info->rows, info->cols, info->loc, info->shape, info->layout, adt.opaque_dtt, adt.arena, static_desc);
    return 0;
}

int dplasma_set_datatype_info( const dplasma_data_collection_t *dc,
                               parsec_arena_datatype_t adt,
                               const lapack_info_t* info)
{
    dplasma_set_dtt_to_info(dc, adt, info);
    dplasma_set_info_to_dtt(dc, adt, info);
    return 0;
}

int dplasma_get_info_from_datatype( const dplasma_data_collection_t *dc, parsec_datatype_t dtt,
                                    const lapack_info_t **info,
                                    const parsec_arena_datatype_t **adt)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;
    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|%"PRIxPTR"|-1|-1|-1", dc, (intptr_t)dtt);
    parsec_key_t k = (uint64_t)static_desc;
    assert(dplasma_datatypes_lapack_helper!= NULL);
    parsec_hash_table_lock_bucket(dplasma_datatypes_lapack_helper, k);  /* protect access to the hash table */
    dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
    assert(dtt_entry != NULL);
    *info = &dtt_entry->info;
    *adt  = &dtt_entry->adt;
    parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);
    return 0;
}

int dplasma_get_datatype_from_info( const dplasma_data_collection_t *dc,
                                    lapack_info_t *info,
                                    const parsec_arena_datatype_t **adt)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;
    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|-1|%d|%d|%d", dc, info->loc, info->shape, info->layout);
    parsec_key_t k = (uint64_t)static_desc;
    assert(dplasma_datatypes_lapack_helper!= NULL);
    parsec_hash_table_lock_bucket(dplasma_datatypes_lapack_helper, k);  /* protect access to the hash table */
    dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
    assert(dtt_entry != NULL);
    *info = dtt_entry->info;
    *adt = &dtt_entry->adt;
    parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);
    return 0;
}

int dplasma_cleanup_datatype_info( const dplasma_data_collection_t *dc,
                                   const lapack_info_t* info,
                                   parsec_arena_datatype_t *adt)
{
    dplasma_datatype_lapack_helper_t * dtt_entry = NULL;
    char static_desc[LAPACK_ADT_KEY_STR_SZ];
    parsec_key_t k;

    assert(dplasma_datatypes_lapack_helper!= NULL);

    /* Remove the two entries from the hash table */
    snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|-1|%d|%d|%d", dc, info->loc, info->shape, info->layout);
    k = (uint64_t)static_desc;
    parsec_hash_table_lock_bucket(dplasma_datatypes_lapack_helper, k);  /* protect access to the info bucket */
    dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
    if(dtt_entry != NULL) {
        *adt = dtt_entry->adt;
        PARSEC_DEBUG_VERBOSE(27, parsec_debug_output,
            " remove dtt <- info lda %d rows %d cols %d loc %d shape %d layout %d dtt %p arena %p (%30s)",
            info->lda, info->rows, info->cols, info->loc, info->shape, info->layout, adt->opaque_dtt, adt->arena, static_desc);
        dplasma_datatypes_lapack_helper_release((void*)dtt_entry, NULL);
        parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);  /* release the loc on the info bucket */

        snprintf( static_desc, LAPACK_ADT_KEY_STR_SZ, "%p|%"PRIxPTR"|-1|-1|-1", dc, (intptr_t)adt->opaque_dtt);
        k = (uint64_t)static_desc;
        parsec_hash_table_lock_bucket(dplasma_datatypes_lapack_helper, k);  /* protect access to the datatype bucket */
        dtt_entry = (dplasma_datatype_lapack_helper_t *)parsec_hash_table_nolock_find(dplasma_datatypes_lapack_helper, k);
        if(dtt_entry != NULL){ /* same desc + dtt can be reuse for several loc & shape */
            PARSEC_DEBUG_VERBOSE(27, parsec_debug_output,
                " remove dtt -> info lda %d rows %d cols %d loc %d shape %d layout %d dtt %"PRIxPTR" arena %p (%30s)",
                dtt_entry->info.lda, dtt_entry->info.rows, dtt_entry->info.cols, dtt_entry->info.loc, dtt_entry->info.shape, dtt_entry->info.layout,
                (intptr_t)adt->opaque_dtt, adt->arena, static_desc);
            dplasma_datatypes_lapack_helper_release((void*)dtt_entry, NULL);
        }
        parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);  /* release the lock on the datatype bucket */
        return 0;
    }
    parsec_hash_table_unlock_bucket(dplasma_datatypes_lapack_helper, k);  /* release the lock on the info bucket */
    /* No found */
    return -1;
}


/***************************************************************************/
/* Management of dplasma_data_collection_t to wrap parsec_data_collection  */
/***************************************************************************/

static parsec_data_t* data_of(parsec_data_collection_t *desc, ...)
{
    int m, n;
    va_list ap;
    int rc;
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)desc;
    (void)rc;

    /* Get coordinates */
    va_start(ap, desc);
    m = (int)va_arg(ap, unsigned int);
    n = (int)va_arg(ap, unsigned int);
    va_end(ap);
    parsec_data_t* dt = ((parsec_data_collection_t*)ddc->dc_original)->data_of(((parsec_data_collection_t *)ddc->dc_original), m, n);
    parsec_data_copy_t *cp = dt->device_copies[0];

    /* Correct datatype if not default because of location */
    int loc = dplasma_tile_location( ddc->dc_original, m, n);

    const lapack_info_t *cp_info;
    const parsec_arena_datatype_t *adt;

    /* obtain the shape of the original copy (aka device[0]) */
    rc = dplasma_get_info_from_datatype(ddc, cp->dtt, &cp_info, &adt);
    assert(rc == 0);
    if( loc != cp_info->loc ) {  /* if we are looking for a different shape/location */
        /* obtain appropriate dtt for that location, shape and layout */
        /* make a copy of the info, don't alter the one currently in the hash table */
        lapack_info_t info = *cp_info;
        info.loc = loc;
        rc = dplasma_get_datatype_from_info(ddc, &info, &adt);
        assert(rc == 0);
        if(cp->dtt != adt->opaque_dtt){
            PARSEC_DEBUG_VERBOSE(27, parsec_debug_output,
                "data_of CP %p [old type %p] loc %d -> dtt %p target_shape %d layout %d",
                cp, cp->dtt, loc, adt->opaque_dtt, info.shape, info.layout);
            dt = parsec_data_create_with_type( dt->dc,
                                               dt->key, cp->device_private, dt->nb_elts,
                                               adt->opaque_dtt);
        }
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

dplasma_data_collection_t *dplasma_wrap_data_collection(parsec_tiled_matrix_t *A)
{
    dplasma_data_collection_t * ddc = (dplasma_data_collection_t *)malloc(sizeof(dplasma_data_collection_t));
    parsec_data_collection_init((parsec_data_collection_t *)ddc,
                                 A->super.nodes,
                                 A->super.myrank );

#if defined(PARSEC_PROF_TRACE)
    if(NULL != A->super.key_dim) {
        asprintf(&ddc->super.key_dim, "%s_dplasmawrapped", A->super.key_dim);
    }
#endif
    if(NULL != A->super.key_base) {
        parsec_data_collection_set_key((parsec_data_collection_t*)ddc, A->super.key_base);
    }

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

