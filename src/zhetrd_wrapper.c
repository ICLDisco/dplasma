/*
 * Copyright (c) 2014-2020 The University of Tennessee and The University
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

#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/private_mempool.h"
#include "parsec/data_dist/matrix/diag_band_to_rect.h"
#include "zhetrd_h2b_L.h"
#include "zhetrd_b2s.h"


/* SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
 */
int
dplasma_zhetrd( parsec_context_t* parsec,
                dplasma_enum_t uplo,
                int ib,
                parsec_tiled_matrix_dc_t* A,
                parsec_tiled_matrix_dc_t* DE,
                parsec_tiled_matrix_dc_t* T,
                int* info )
{
    parsec_zhetrd_h2b_L_taskpool_t * h2b = NULL;
    parsec_diag_band_to_rect_taskpool_t* band2rect = NULL;
    parsec_zhetrd_b2s_taskpool_t * b2s = NULL;
    parsec_memory_pool_t pool[4];

    if( uplo != dplasmaLower && uplo != dplasmaUpper ) {
        dplasma_error("DPLASMA_zhetrd", "illegal value of uplo");
        *info = -1;
        return *info;
    }

    parsec_private_memory_init( &pool[0], (sizeof(dplasma_complex64_t)*T->nb) ); /* tau */
    parsec_private_memory_init( &pool[1], (sizeof(dplasma_complex64_t)*T->nb*ib) ); /* work */
    parsec_private_memory_init( &pool[2], (sizeof(dplasma_complex64_t)*T->nb*2 *T->nb) ); /* work for HERFB1 */
    parsec_private_memory_init( &pool[3], (sizeof(dplasma_complex64_t)*T->nb*4 *T->nb) ); /* work for the TSMQRLR */

    if( dplasmaLower == uplo ) {
        h2b = parsec_zhetrd_h2b_L_new( ib, A, T, &pool[3], &pool[2], &pool[1], &pool[0] );
        dplasma_add2arena_rectangle( &h2b->arenas_datatypes[PARSEC_zhetrd_h2b_L_DEFAULT_ARENA],
                                     A->mb*A->nb*sizeof(dplasma_complex64_t),
                                     PARSEC_ARENA_ALIGNMENT_SSE,
                                     parsec_datatype_double_complex_t, A->mb, A->nb, -1 );
        dplasma_add2arena_rectangle( &h2b->arenas_datatypes[PARSEC_zhetrd_h2b_L_LITTLE_T_ARENA],
                                     T->mb*T->nb*sizeof(dplasma_complex64_t),
                                     PARSEC_ARENA_ALIGNMENT_SSE,
                                     parsec_datatype_double_complex_t, T->mb, T->nb, -1 );
#if 0
    } else {
        h2b = parsec_zhetrd_h2b_U_new( ib, A, *A, T, *T, pool[3], pool[2], pool[1], pool[0] );
        dplasma_add2arena_rectangle( &h2b->arenas_datatypes[PARSEC_zhetrd_h2b_U_DEFAULT_ARENA],
                                     A->mb*A->nb*sizeof(dplasma_complex64_t),
                                     PARSEC_ARENA_ALIGNMENT_SSE,
                                     parsec_datatype_double_complex_t, A->mb, A->nb, -1 );
        dplasma_add2arena_rectangle( &h2b->arenas_datatypes[PARSEC_zhetrd_h2b_U_LITTLE_T_ARENA],
                                     T->mb*T->nb*sizeof(dplasma_complex64_t),
                                     PARSEC_ARENA_ALIGNMENT_SSE,
                                     parsec_datatype_double_complex_t, T->mb, T->nb, -1 );
#endif
    }
    if( NULL == h2b ) {
        *info=-101; goto cleanup;
    }

    band2rect = parsec_diag_band_to_rect_new((sym_two_dim_block_cyclic_t*)A, (two_dim_block_cyclic_t*)DE,
                                             A->mt, A->nt, A->mb, A->nb, sizeof(dplasma_complex64_t));
    if( NULL == band2rect ) goto cleanup;
    dplasma_add2arena_tile( &band2rect->arenas_datatypes[PARSEC_diag_band_to_rect_DEFAULT_ARENA],
                            A->mb*A->nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, A->mb );

    b2s = parsec_zhetrd_b2s_new( DE, DE->mb-1 );
    if( NULL == b2s ) goto cleanup;
    dplasma_add2arena_rectangle( &b2s->arenas_datatypes[PARSEC_zhetrd_b2s_DEFAULT_ARENA],
                                 DE->mb*DE->nb*sizeof(dplasma_complex64_t),
                                 PARSEC_ARENA_ALIGNMENT_SSE,
                                 parsec_datatype_double_complex_t, DE->mb, DE->nb, -1 );

    parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)h2b );
    parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)band2rect );
    parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)b2s );
    dplasma_wait_until_completion(parsec);

cleanup:
    if( h2b ) parsec_taskpool_free( &h2b->super );
    if( band2rect ) parsec_taskpool_free( &band2rect->super );
    if( b2s ) parsec_taskpool_free( &b2s->super );
    parsec_private_memory_fini( &pool[0] );
    parsec_private_memory_fini( &pool[1] );
    parsec_private_memory_fini( &pool[2] );
    parsec_private_memory_fini( &pool[3] );
    return *info;
}


