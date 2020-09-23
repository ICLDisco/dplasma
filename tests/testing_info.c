#include <stdio.h>
#include <string.h>
#include "utils/dplasma_info.h"

int main(int argc, char *argv[])
{
    dplasma_info_t info;
    int nkeys, i, present, errors = 0, foundb=0, foundc=0;
    char key[DPLASMA_MAX_INFO_KEY];
    char value[DPLASMA_MAX_INFO_VAL];

    (void)argc;
    (void)argv;

    dplasma_info_create(&info);
    dplasma_info_set(info, "KEYA", "VALUEA");
    dplasma_info_set(info, "KEYB", "VALUEB1");
    dplasma_info_set(info, "KEYC", "VALUEC");
    dplasma_info_set(info, "KEYB", "VALUEB2");

    dplasma_info_delete(info, "KEYA");
    dplasma_info_delete(info, "NONKEY");

    dplasma_info_get_nkeys(info, &nkeys);
    fprintf(stderr, "There are %d keys (should be 2: KEYB and KEYC)\n", nkeys);
    if(nkeys != 2) errors++;
    for(i = 0; i < nkeys; i++) {
        dplasma_info_get_nthkey(info, i, key);
        fprintf(stderr, "Key of index %d is '%s'\n", i, key);
        if(strcmp(key, "KEYB")) {
            if (strcmp(key, "KEYC")) {
                fprintf(stderr, "Unexpected key found (this is an error...)\n");
                errors++;
            } else {
                if (foundc) {
                    fprintf(stderr, "KEYC is found twice or more (this is an error...)\n");
                    errors++;
                } else {
                    fprintf(stderr, "'KEYC' found (this is correct)\n");
                }
                foundc = 1;
            }
        } else {
            if(foundb) {
                fprintf(stderr, "KEYB is found twice or more (this is an error...)\n");
                errors++;
            } else {
                fprintf(stderr, "'KEYB' found (this is correct)\n");
            }
            foundb=1;
        }
    }

    dplasma_info_get(info, "KEYA", DPLASMA_MAX_INFO_VAL, value, &present);
    if(present) {
        fprintf(stderr, "KEYA is present (this is an error)... With value '%s'\n", value);
        errors++;
    } else
        fprintf(stderr, "KEYA is not present (this is correct)\n");

    dplasma_info_get(info, "KEYB", DPLASMA_MAX_INFO_VAL, value, &present);
    if(present) {
        fprintf(stderr, "KEYB is present (this is correct) ");
        if(strcmp(value, "VALUEB2")) {
            fprintf(stderr, "with value '%s' (this is an error...)\n", value);
            errors++;
        } else {
            fprintf(stderr, "with value '%s' (this is correct)\n", value);
        }
    } else {
        fprintf(stderr, "KEYB is not present (this is an error...)\n");
        errors++;
    }

     dplasma_info_get(info, "KEYC", DPLASMA_MAX_INFO_VAL, value, &present);
    if(present) {
        fprintf(stderr, "KEYC is present (this is correct) ");
        if(strcmp(value, "VALUEC")) {
            fprintf(stderr, "with value '%s' (this is an error...)\n", value);
            errors++;
        } else {
            fprintf(stderr, "with value '%s' (this is correct)\n", value);
        }
    } else {
        fprintf(stderr, "KEYC is not present (this is an error...)\n");
        errors++;
    }

    dplasma_info_free(&info);

    return errors;
}
