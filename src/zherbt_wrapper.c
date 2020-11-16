/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */


#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasmaaux.h"
#include "cores/core_blas.h"

#include "parsec/private_mempool.h"

#include "zherbt_L.h"

parsec_taskpool_t *
dplasma_zherbt_New( dplasma_enum_t uplo, int IB,
                    parsec_tiled_matrix_dc_t *A,
                    parsec_tiled_matrix_dc_t *T)
{
    parsec_zherbt_L_taskpool_t *parsec_zherbt = NULL;
    parsec_memory_pool_t *pool[4];

    if( dplasmaLower != uplo ) {
        dplasma_error("dplasma_zherbt_New", "illegal value of uplo");
        return NULL;
    }

    pool[0] = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));  /* tau */
    parsec_private_memory_init( pool[0], (sizeof(dplasma_complex64_t)*T->nb) );

    pool[1] = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));  /* work */
    parsec_private_memory_init( pool[1], (sizeof(dplasma_complex64_t)*T->nb*IB) );

    pool[2] = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));  /* work for HERFB1 */
    parsec_private_memory_init( pool[2], (sizeof(dplasma_complex64_t)*T->nb*2 *T->nb) );

    pool[3] = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));  /* work for the TSMQRLR */
    parsec_private_memory_init( pool[3], (sizeof(dplasma_complex64_t)*T->nb*4 *T->nb) );

    if( dplasmaLower == uplo ) {
        parsec_zherbt = parsec_zherbt_L_new(uplo, IB,
                                          A,
                                          T,
                                          pool[3], pool[2], pool[1], pool[0]);
        dplasma_add2arena_rectangle( &parsec_zherbt->arenas_datatypes[PARSEC_zherbt_L_DEFAULT_ADT_IDX],
                                     A->mb*A->nb*sizeof(dplasma_complex64_t),
                                     PARSEC_ARENA_ALIGNMENT_SSE,
                                     parsec_datatype_double_complex_t, A->mb, A->nb, -1);
        dplasma_add2arena_rectangle( &parsec_zherbt->arenas_datatypes[PARSEC_zherbt_L_LITTLE_T_ADT_IDX],
                                     T->mb*T->nb*sizeof(dplasma_complex64_t),
                                     PARSEC_ARENA_ALIGNMENT_SSE,
                                     parsec_datatype_double_complex_t, T->mb, T->nb, -1);

    }

    return (parsec_taskpool_t*)parsec_zherbt;
}

void dplasma_zherbt_Destruct( parsec_taskpool_t *tp )
{
    parsec_zherbt_L_taskpool_t *parsec_zherbt = (parsec_zherbt_L_taskpool_t *)tp;

    if( dplasmaLower == parsec_zherbt->_g_uplo ) {

        dplasma_matrix_del2arena( &parsec_zherbt->arenas_datatypes[PARSEC_zherbt_L_DEFAULT_ADT_IDX   ] );
        dplasma_matrix_del2arena( &parsec_zherbt->arenas_datatypes[PARSEC_zherbt_L_LITTLE_T_ADT_IDX  ] );

        parsec_private_memory_fini( parsec_zherbt->_g_pool_0 );
        free( parsec_zherbt->_g_pool_0 );
        parsec_private_memory_fini( parsec_zherbt->_g_pool_1 );
        free( parsec_zherbt->_g_pool_1 );
        parsec_private_memory_fini( parsec_zherbt->_g_pool_2 );
        free( parsec_zherbt->_g_pool_2 );
        parsec_private_memory_fini( parsec_zherbt->_g_pool_3 );
        free( parsec_zherbt->_g_pool_3 );

        parsec_taskpool_free(tp);
    }
}
