/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "dplasma.h"
#include "cores/core_blas.h"

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_zgesv - Solves a system of linear equations A * X = B with a general
 * square matrix A using the LU factorization with incremental pivoting strategy
 * computed by dplasma_zgetrf_incpiv().
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A to be factorized.
 *          On entry, describes the M-by-N matrix A.
 *          On exit, elements on and above the diagonal are the elements of
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
 * @sa dplasma_zgesv_New
 * @sa dplasma_zgesv_Destruct
 * @sa dplasma_cgesv
 * @sa dplasma_dgesv
 * @sa dplasma_sgesv
 *
 ******************************************************************************/
int
dplasma_zgesv_incpiv( parsec_context_t *parsec,
                      parsec_tiled_matrix_dc_t *A,
                      parsec_tiled_matrix_dc_t *L,
                      parsec_tiled_matrix_dc_t *IPIV,
                      parsec_tiled_matrix_dc_t *B )
{
    int info;

#ifdef PARSEC_COMPOSITION
    parsec_taskpool_t *parsec_zgetrf  = dplasma_zgetrf_incpiv_New(A, L, IPIV, &info);
    parsec_taskpool_t *parsec_ztrsmpl = dplasma_ztrsmpl_New(A, L, IPIV, B);
    parsec_taskpool_t *parsec_ztrsm   = dplasma_ztrsm_New(dplasmaLeft, dplasmaUpper, dplasmaNoTrans, dplasmaNonUnit, 1.0, A, B);

    parsec_context_add_taskpool( parsec, parsec_zgetrf  );
    parsec_context_add_taskpool( parsec, parsec_ztrsmpl );
    parsec_context_add_taskpool( parsec, parsec_ztrsm   );

    dplasma_wait_until_completion( parsec );

    dplasma_zgetrf_incpiv_Destruct( parsec_zgetrf  );
    dplasma_ztrsmpl_Destruct( parsec_ztrsmpl );
    dplasma_ztrsm_Destruct( parsec_ztrsm   );
#else
    info = dplasma_zgetrf_incpiv(parsec, A, L, IPIV );
    if( info == 0 ) {
        dplasma_zgetrs_incpiv(parsec, dplasmaNoTrans, A, L, IPIV, B );
    }
#endif

    return info;
}
