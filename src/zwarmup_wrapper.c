/*
 * Copyright (c) 2022      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * $COPYRIGHT
 *
 * @precisions normal z -> s d c
 *
 */

#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasma/types_lapack.h"
#include "dplasmaaux.h"
#include "parsec/utils/zone_malloc.h"

#include "cores/dplasma_plasmatypes.h"
#include <parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h>
#include "zwarmup.h"

#define MAX_SHAPES 1

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zwarmup_New - Generates the taskpool that loads all data to GPU once
 *
 * WARNING: The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A.
 *          On exit A is unmodified.
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_zwarmup_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_zwarmup
 * @sa dplasma_zwarmup_Destruct
 *r @sa dplasma_cwarmup_New
 * @sa dplasma_dwarmup_New
 * @sa dplasma_swarmup_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_zwarmup_New(parsec_tiled_matrix_t *A)
{
    parsec_taskpool_t *tp = NULL;
    dplasma_data_collection_t * ddc_A = dplasma_wrap_data_collection(A);

    int shape = 0;
    dplasma_setup_adtt_all_loc( ddc_A,
                                parsec_datatype_double_complex_t,
                                PARSEC_MATRIX_FULL/*uplo*/, 1/*diag:for PARSEC_MATRIX_UPPER or PARSEC_MATRIX_LOWER types*/,
                                &shape);

    tp = (parsec_taskpool_t*)parsec_zwarmup_new( ddc_A );
    assert(shape == MAX_SHAPES);
    return tp;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zwarmup_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_zwarmup_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_zwarmup_New
 * @sa dplasma_zwarmup
 *
 ******************************************************************************/
void
dplasma_zwarmup_Destruct( parsec_taskpool_t *tp )
{
    parsec_zwarmup_taskpool_t *parsec_zwarmup = (parsec_zwarmup_taskpool_t *)tp;
    dplasma_clean_adtt_all_loc(parsec_zwarmup->_g_ddescA, MAX_SHAPES);
    dplasma_data_collection_t * ddc_A = parsec_zwarmup->_g_ddescA;

    parsec_taskpool_free(tp);
    /* free the dplasma_data_collection_t */
    dplasma_unwrap_data_collection(ddc_A);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa dplasma_zwarmup_New
 * @sa dplasma_zwarmup_Destruct
 * @sa dplasma_cwarmup
 * @sa dplasma_dwarmup
 * @sa dplasma_swarmup
 *
 ******************************************************************************/
int
dplasma_zwarmup( parsec_context_t *parsec,
                parsec_tiled_matrix_t *A )
{
    parsec_taskpool_t *parsec_zwarmup = NULL;

    parsec_zwarmup = dplasma_zwarmup_New( A );

    if ( parsec_zwarmup != NULL ) {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zwarmup);
        dplasma_wait_until_completion(parsec);
        dplasma_zwarmup_Destruct( parsec_zwarmup );
    }

    return DPLASMA_SUCCESS;
}

