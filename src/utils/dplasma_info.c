/*
 * Copyright (c) 2021      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */

#include <string.h>
#include <stdio.h>
#include "utils/dplasma_info.h"

typedef struct {
    parsec_hash_table_item_t ht_item;
    char *name;
    char *value;
} dplasma_info_item_t;

static int dplasma_info_key_equal(parsec_key_t a, parsec_key_t b, void *data)
{
    char *as = (char *)a;
    char *bs = (char *)b;
    (void)data;
    return 0 == strcmp(as, bs);
}

static uint64_t dplasma_info_key_hash(parsec_key_t k, void *data)
{
    uint64_t key = 5381, c;
    char *ks = (char *)k;
    (void)data;

    while( (c = *ks++) )
        key = (key * 33) ^ c;
    return key;
}

static char *dplasma_info_key_print(char *buffer, size_t buffer_len, parsec_key_t k, void *data)
{
    char *ks = (char *)k;
    (void)data;
    strncpy(buffer, ks, buffer_len);
    return buffer;
}

static parsec_key_fn_t dplasma_info_hash_functions =
        { .key_hash  = dplasma_info_key_hash,
          .key_equal = dplasma_info_key_equal,
          .key_print = dplasma_info_key_print };

void dplasma_info_create(dplasma_info_t *info)
{
    if(NULL == info)
        return;
    *info = malloc(sizeof(struct dplasma_info_s));
    PARSEC_OBJ_CONSTRUCT( &((*info)->info_table), parsec_hash_table_t );
    parsec_hash_table_init(&((*info)->info_table), offsetof(dplasma_info_item_t, ht_item), 6,
                           dplasma_info_hash_functions, NULL);
    (*info)->nb_elts = 0;
}

void dplasma_info_set(dplasma_info_t info, const char *name, const char *value)
{
    dplasma_info_item_t *info_item;
    assert(NULL != name);
    parsec_hash_table_lock_bucket(&info->info_table, (parsec_key_t)name);
    info_item = (dplasma_info_item_t*)parsec_hash_table_nolock_find(&info->info_table, (parsec_key_t)name);
    if(NULL == info_item) {
        dplasma_info_item_t *new_item = malloc(sizeof(dplasma_info_item_t));
        new_item->name = strdup(name);
        new_item->value = strdup(value);
        new_item->ht_item.key = (parsec_key_t)new_item->name;
        parsec_hash_table_nolock_insert(&info->info_table, &new_item->ht_item);
        info->nb_elts++;
    } else {
        free(info_item->value);
        info_item->value = strdup(value);
    }
    parsec_hash_table_unlock_bucket(&info->info_table, (parsec_key_t)name);
}

void dplasma_info_get_nkeys(dplasma_info_t info, int *nkeys)
{
    *nkeys = info->nb_elts;
}

typedef struct {
    int key_index;
    char *key;
} dplasma_info_key_finder_arg_t;

static void dplasma_info_key_finder(void *item, void *cb_data)
{
    dplasma_info_key_finder_arg_t *arg = (dplasma_info_key_finder_arg_t*)cb_data;
    dplasma_info_item_t *ii = (dplasma_info_item_t*)item;
    if(arg->key_index == 0)
        strncpy(arg->key, ii->name, DPLASMA_MAX_INFO_KEY);
    arg->key_index--;
}

void dplasma_info_get_nthkey(dplasma_info_t info, int key_index, char *key)
{
    dplasma_info_key_finder_arg_t arg;
    assert(key_index >= 0);
    assert(key_index < info->nb_elts);
    arg.key_index = key_index;
    arg.key = key;
    parsec_hash_table_for_all(&info->info_table, dplasma_info_key_finder, &arg);
}

void dplasma_info_get(dplasma_info_t info, const char *key, int valuelen, char *value, int *flag)
{
    parsec_key_t k = (parsec_key_t)key;
    dplasma_info_item_t *item = (dplasma_info_item_t *)parsec_hash_table_find(&info->info_table, k);
    if(NULL == item) {
        if (NULL != flag)
            *flag = 0;
        return;
    }
    if(NULL != flag)
        *flag = 1;
    strncpy(value, item->value, valuelen);
}

void dplasma_info_delete(dplasma_info_t info, const char *key)
{
    parsec_key_t k = (parsec_key_t)key;
    dplasma_info_item_t *item;
    parsec_hash_table_lock_bucket(&info->info_table, k);
    item = (dplasma_info_item_t *)parsec_hash_table_nolock_remove(&info->info_table, k);
    if( NULL != item) {
        info->nb_elts--;
        parsec_hash_table_unlock_bucket(&info->info_table, k);
        free(item->name);
        free(item->value);
        free(item);
    } else {
        parsec_hash_table_unlock_bucket(&info->info_table, k);
    }
}

static void dplasma_info_key_remover(void *item, void *cb_data)
{
    dplasma_info_t info = (dplasma_info_t)cb_data;
    dplasma_info_item_t *ii = (dplasma_info_item_t *)item;
    parsec_hash_table_nolock_remove(&info->info_table, (parsec_key_t)ii->name);
    info->nb_elts--;
    free(ii->name);
    free(ii->value);
    free(ii);
}

void dplasma_info_free(dplasma_info_t *nfo)
{
    dplasma_info_t info = *nfo;
    parsec_hash_table_for_all(&info->info_table, dplasma_info_key_remover, info);
    PARSEC_OBJ_DESTRUCT(&info->info_table);
    free(info);
    *nfo = NULL;
}
