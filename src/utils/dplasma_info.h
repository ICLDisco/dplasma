#ifndef _DPLASMA_INFO_H
#define _DPLASMA_INFO_H

#include "parsec/class/parsec_hash_table.h"

typedef struct dplasma_info_s {
    parsec_hash_table_t info_table;
    int                 nb_elts;
} *dplasma_info_t;

#define DPLASMA_MAX_INFO_KEY  64
#define DPLASMA_MAX_INFO_VAL 128

void dplasma_info_create(dplasma_info_t *info);
void dplasma_info_set(dplasma_info_t info, const char *name, const char *value);
void dplasma_info_get_nkeys(dplasma_info_t info, int *nkeys);
void dplasma_info_get_nthkey(dplasma_info_t info, int key_index, char *key);
void dplasma_info_get(dplasma_info_t info, const char *key, int valuelen, char *value, int *flag);
void dplasma_info_delete(dplasma_info_t info, const char *key);
void dplasma_info_free(dplasma_info_t *info);

#endif