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
 * dplasma_zpotri - Computes the inverse of a complex Hermitian positive definite
 * matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
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
 * @param[in,out] A
 *          Descriptor of the distributed factorized matrix A.
 *
 *******************************************************************************
 *
 * @return
 *          \retval 0 on success.
 *          \retval -i if the ith parameters is incorrect.
 *          \retval >0 if i, the leading minor of order i of A is not
 *               positive definite, so the factorization could not be
 *               completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa dplasma_zpotri_New
 * @sa dplasma_zpotri_Destruct
 * @sa dplasma_cpotri
 * @sa dplasma_dpotri
 * @sa dplasma_spotri
 *
 ******************************************************************************/
int
dplasma_zpotri( parsec_context_t *parsec,
                dplasma_enum_t uplo,
                parsec_tiled_matrix_t* A )
{
    int info = 0;
    /* Check input arguments */
    if (uplo != dplasmaUpper && uplo != dplasmaLower) {
        dplasma_error("dplasma_zpotri", "illegal value of uplo");
        return -1;
    }

#ifdef PARSEC_COMPOSITION
    parsec_taskpool_t *parsec_ztrtri = NULL;
    parsec_taskpool_t *parsec_zlauum = NULL;

    parsec_ztrtri = dplasma_ztrtri_New(uplo, dplasmaNonUnit, A, &info );
    parsec_zlauum = dplasma_zlauum_New(uplo, A );

    parsec_context_add_taskpool( parsec, parsec_ztrtri );
    parsec_context_add_taskpool( parsec, parsec_zlauum );

    dplasma_wait_until_completion( parsec );

    dplasma_ztrtri_Destruct( parsec_ztrtri );
    dplasma_zlauum_Destruct( parsec_zlauum );
#else
    info = dplasma_ztrtri( parsec, uplo, dplasmaNonUnit, A );
    dplasma_zlauum( parsec, uplo, A );
#endif
    return info;
}

