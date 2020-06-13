/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 *
 * @precisions normal z -> c d s
 *
 */

#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasmaaux.h"

#include "zprint.h"

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zprint - Print a matrix tile by tile.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] uplo
 *          Specifies which part of the matrix is printed
 *          = dplasmaUpper: Upper part of A;
 *          = dplasmaLower: Lower part of A;
 *          = dplasmaUpperLower: ALL elements of A.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A to generate. Any tiled matrix
 *          descriptor can be used.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa dplasma_cprint
 * @sa dplasma_dprint
 * @sa dplasma_sprint
 *
 ******************************************************************************/
int dplasma_zprint( parsec_context_t *parsec,
                    dplasma_enum_t uplo,
                    const parsec_tiled_matrix_dc_t *A)
{
    parsec_zprint_taskpool_t* tp;

    /* Check input arguments */
    if ((uplo != dplasmaLower) &&
        (uplo != dplasmaUpper) &&
        (uplo != dplasmaUpperLower))
    {
        dplasma_error("dplasma_zplghe", "illegal value of type");
        return -3;
    }

    tp = parsec_zprint_new( uplo, A);

    if (tp != NULL) {
        /* Default type */
        dplasma_add2arena_tile( &tp->arenas_datatypes[PARSEC_zprint_DEFAULT_ARENA],
                                A->mb*A->nb*sizeof(dplasma_complex64_t),
                                PARSEC_ARENA_ALIGNMENT_SSE,
                                parsec_datatype_double_complex_t, A->mb );

        parsec_context_add_taskpool(parsec, (parsec_taskpool_t*)tp);
        dplasma_wait_until_completion(parsec);

        dplasma_matrix_del2arena( &tp->arenas_datatypes[PARSEC_zprint_DEFAULT_ARENA] );
        parsec_taskpool_free( &tp->super );
        return 0;
    }
    return -101;
}
