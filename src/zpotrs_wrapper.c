/*
 * Copyright (c) 2010-2018 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
 *
 * @precisions normal z -> s d c
 *
 */

#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasmaaux.h"

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zpotrs - Solves a system of linear equations A * X = B with a
 * symmetric positive definite (or Hermitian positive definite in the complex
 * case) matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
 * computed by dplasma_zpotrf().
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] uplo
 *          = dplasmaUpper: Upper triangle of A is referenced;
 *          = dplasmaLower: Lower triangle of A is referenced.
 *
 * @param[in] A
 *          Descriptor of the distributed factorized matrix A.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, B is overwritten by the solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa dplasma_zpotrs_New
 * @sa dplasma_zpotrs_Destruct
 * @sa dplasma_cpotrs
 * @sa dplasma_dpotrs
 * @sa dplasma_spotrs
 *
 ******************************************************************************/
int
dplasma_zpotrs( parsec_context_t *parsec,
                dplasma_enum_t uplo,
                const parsec_tiled_matrix_dc_t* A,
                parsec_tiled_matrix_dc_t* B )
{
    /* Check input arguments */
    if (uplo != dplasmaUpper && uplo != dplasmaLower) {
        dplasma_error("dplasma_zpotrs", "illegal value of uplo");
        return -1;
    }

#ifdef PARSEC_COMPOSITION
    parsec_taskpool_t *parsec_ztrsm1 = NULL;
    parsec_taskpool_t *parsec_ztrsm2 = NULL;

    if ( uplo == dplasmaUpper ) {
      parsec_ztrsm1 = dplasma_ztrsm_New(dplasmaLeft, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0, A, B);
      parsec_ztrsm2 = dplasma_ztrsm_New(dplasmaLeft, uplo, dplasmaNoTrans,   dplasmaNonUnit, 1.0, A, B);
    } else {
      parsec_ztrsm1 = dplasma_ztrsm_New(dplasmaLeft, uplo, dplasmaNoTrans,   dplasmaNonUnit, 1.0, A, B);
      parsec_ztrsm2 = dplasma_ztrsm_New(dplasmaLeft, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0, A, B);
    }

    parsec_context_add_taskpool( parsec, parsec_ztrsm1 );
    parsec_context_add_taskpool( parsec, parsec_ztrsm2 );

    dplasma_wait_until_completion( parsec );

    dplasma_ztrsm_Destruct( parsec_ztrsm1 );
    dplasma_ztrsm_Destruct( parsec_ztrsm2 );
#else
    if ( uplo == dplasmaUpper ) {
      dplasma_ztrsm( parsec, dplasmaLeft, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0, A, B );
      dplasma_ztrsm( parsec, dplasmaLeft, uplo, dplasmaNoTrans,   dplasmaNonUnit, 1.0, A, B );
    } else {
      dplasma_ztrsm( parsec, dplasmaLeft, uplo, dplasmaNoTrans,   dplasmaNonUnit, 1.0, A, B );
      dplasma_ztrsm( parsec, dplasmaLeft, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0, A, B );
    }
#endif
    return 0;
}

