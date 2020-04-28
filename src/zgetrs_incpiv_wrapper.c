/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "dplasma.h"
#include "dplasmaaux.h"
#include "cores/core_blas.h"

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zgetrs - Solves a system of linear equations A * X = B with a general
 * square matrix A using the LU factorization with incremental pivoting strategy
 * computed by dplasma_zgetrf_incpiv().
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          conjugate transposed:
 *          = dplasmaNoTrans:   A is transposed;
 *          = dplasmaTrans:     A is not transposed;
 *          = dplasmaConjTrans: A is conjugate transposed.
 *          Currently only dplasmaNoTrans is supported.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A to be factorized.
 *          On entry, The factorized matrix through dplasma_zgetrf_incpiv_New()
 *          routine.  Elements on and above the diagonal are the elements of
 *          U. Elements belowe the diagonal are NOT the classic L, but the L
 *          factors obtaines by succesive pivoting.
 *
 * @param[in] L
 *          Descriptor of the matrix L distributed exactly as the A matrix.
 *           - If IPIV != NULL, L.mb defines the IB parameter of the tile LU
 *          algorithm. This matrix must be of size A.mt * L.mb - by - A.nt *
 *          L.nb, with L.nb == A.nb.
 *          On entry, contains auxiliary information required to solve the
 *          system and generated by dplasma_zgetrf_inciv_New().
 *           - If IPIV == NULL, pivoting information are stored within
 *          L. (L.mb-1) defines the IB parameter of the tile LU algorithm. This
 *          matrix must be of size A.mt * L.mb - by - A.nt * L.nb, with L.nb =
 *          A.nb, and L.mb = ib+1.
 *          The first A.mb elements contains the IPIV information, the leftover
 *          contains auxiliary information required to solve the system.
 *
 * @param[in] IPIV
 *          Descriptor of the IPIV matrix. Should be distributed exactly as the
 *          A matrix. This matrix must be of size A.m - by - A.nt with IPIV.mb =
 *          A.mb and IPIV.nb = 1.
 *          On entry, contains the pivot indices of the successive row
 *          interchanged performed during the factorization.
 *          If IPIV == NULL, rows interchange information is stored within L.
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
 * @sa dplasma_zgetrs_New
 * @sa dplasma_zgetrs_Destruct
 * @sa dplasma_cgetrs
 * @sa dplasma_dgetrs
 * @sa dplasma_sgetrs
 *
 ******************************************************************************/
int
dplasma_zgetrs_incpiv(parsec_context_t *parsec,
                      dplasma_enum_t trans,
                      parsec_tiled_matrix_dc_t *A,
                      parsec_tiled_matrix_dc_t *L,
                      parsec_tiled_matrix_dc_t *IPIV,
                      parsec_tiled_matrix_dc_t *B)
{
    /* Check input arguments */
    if (trans != dplasmaNoTrans) {
        dplasma_error("dplasma_zgetrs", "only dplasmaNoTrans supported");
        return -1;
    }

#ifdef PARSEC_COMPOSITION
    parsec_taskpool_t *parsec_ztrsmpl = NULL;
    parsec_taskpool_t *parsec_ztrsm   = NULL;

    parsec_ztrsmpl = dplasma_ztrsmpl_New(A, L, IPIV, B);
    parsec_ztrsm   = dplasma_ztrsm_New(dplasmaLeft, dplasmaUpper, dplasmaNoTrans, dplasmaNonUnit, 1.0, A, B);

    parsec_context_add_taskpool( parsec, parsec_ztrsmpl );
    parsec_context_add_taskpool( parsec, parsec_ztrsm   );

    dplasma_wait_until_completion( parsec );

    dplasma_ztrsm_Destruct( parsec_ztrsmpl );
    dplasma_ztrsm_Destruct( parsec_ztrsm   );
#else
    dplasma_ztrsmpl(parsec, A, L, IPIV, B );
    dplasma_ztrsm( parsec, dplasmaLeft, dplasmaUpper, dplasmaNoTrans, dplasmaNonUnit, 1.0, A, B );
#endif
    return 0;
}

