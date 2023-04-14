/*
 * Copyright (c) 2011-2023 The University of Tennessee and The University
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

#include "parsec/data_dist/matrix/matrix.h"
#include "parsec/private_mempool.h"

#include "zhbrdt.h"

parsec_taskpool_t* dplasma_zhbrdt_New(parsec_tiled_matrix_t* A /* data A */)
{
    parsec_zhbrdt_taskpool_t *parsec_zhbrdt = NULL;

    parsec_zhbrdt = parsec_zhbrdt_new(A, A->mb-1);

    dplasma_add2arena_rectangle( &parsec_zhbrdt->arenas_datatypes[PARSEC_zhbrdt_DEFAULT_ADT_IDX],
                                 (A->nb)*(A->mb)*sizeof(dplasma_complex64_t), 16,
                                 parsec_datatype_double_complex_t,
                                 A->mb, A->nb, -1 );

    return (parsec_taskpool_t*)parsec_zhbrdt;
}

void dplasma_zhbrdt_Destruct( parsec_taskpool_t *tp )
{
    parsec_zhbrdt_taskpool_t *zhbrdt_tp = (parsec_zhbrdt_taskpool_t*)tp;
    dplasma_matrix_del2arena( &zhbrdt_tp->arenas_datatypes[PARSEC_zhbrdt_DEFAULT_ADT_IDX] );
    parsec_taskpool_free(tp);
}

