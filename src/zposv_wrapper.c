/*
 * Copyright (c) 2010-2022 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
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
 * dplasma_zposv - Computes the solution to a system of linear equations A * X =
 * B, where A is an N-by-N symmetric positive definite (or Hermitian positive
 * definite in the complex case) matrix and X and B are N-by-NRHS matrices.  The
 * Cholesky decomposition is used to factor A as
 *
 *    \f[ A = \{_{L\times L^H, if uplo = dplasmaLower}^{U^H\times U, if uplo = dplasmaUpper} \f]
 *
 * where U is an upper triangular matrix and  L is a lower triangular matrix.
 * The factored form of A is then used to solve the system of equations A * X = B.
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
 * @param[in,out] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten with the factorized
 *          matrix.
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
 * @sa dplasma_zposv_New
 * @sa dplasma_zposv_Destruct
 * @sa dplasma_cposv
 * @sa dplasma_dposv
 * @sa dplasma_sposv
 *
 ******************************************************************************/
int
dplasma_zposv( parsec_context_t *parsec,
               dplasma_enum_t uplo,
               parsec_tiled_matrix_t *A,
               parsec_tiled_matrix_t *B )
{
    int info;

    /* Check input arguments */
    if (uplo != dplasmaUpper && uplo != dplasmaLower) {
        dplasma_error("dplasma_zposv", "illegal value of uplo");
        return -1;
    }

#ifdef PARSEC_COMPOSITION
    parsec_taskpool_t *parsec_ztrsm1 = NULL;
    parsec_taskpool_t *parsec_ztrsm2 = NULL;
    parsec_taskpool_t *parsec_zpotrf = NULL;

    parsec_zpotrf = dplasma_zpotrf_New(uplo, A, &info);
    if ( uplo == dplasmaUpper ) {
      parsec_ztrsm1 = dplasma_ztrsm_New(dplasmaLeft, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0, A, B);
      parsec_ztrsm2 = dplasma_ztrsm_New(dplasmaLeft, uplo, dplasmaNoTrans,   dplasmaNonUnit, 1.0, A, B);
    } else {
      parsec_ztrsm1 = dplasma_ztrsm_New(dplasmaLeft, uplo, dplasmaNoTrans,   dplasmaNonUnit, 1.0, A, B);
      parsec_ztrsm2 = dplasma_ztrsm_New(dplasmaLeft, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0, A, B);
    }

    parsec_context_add_taskpool( parsec, parsec_zpotrf );
    parsec_context_add_taskpool( parsec, parsec_ztrsm1 );
    parsec_context_add_taskpool( parsec, parsec_ztrsm2 );

    dplasma_wait_until_completion( parsec );

    dplasma_zpotrf_Destruct( parsec_zpotrf );
    dplasma_ztrsm_Destruct( parsec_ztrsm1 );
    dplasma_ztrsm_Destruct( parsec_ztrsm2 );
#else
    info = dplasma_zpotrf( parsec, uplo, A);
    if ( info == 0 ) {
      if ( uplo == dplasmaUpper ) {
        dplasma_ztrsm( parsec, dplasmaLeft, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0, A, B );
        dplasma_ztrsm( parsec, dplasmaLeft, uplo, dplasmaNoTrans,   dplasmaNonUnit, 1.0, A, B );
      } else {
        dplasma_ztrsm( parsec, dplasmaLeft, uplo, dplasmaNoTrans,   dplasmaNonUnit, 1.0, A, B );
        dplasma_ztrsm( parsec, dplasmaLeft, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0, A, B );
      }
    }
#endif
    return info;
}
