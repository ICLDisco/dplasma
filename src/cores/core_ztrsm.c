/**
 *
 * @file core_ztrsm.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "core_blas.h"

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_ztrsm - Computes triangular solve A*X = B or X*A = B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  A*X = B
 *          = PlasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] M
 *          The order of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of right hand sides, i.e., the number of columns of the matrix B. N >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading M-by-M upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading M-by-M
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in,out] B
 *          On entry, the M-by-N right hand side matrix B.
 *          On exit, if return value = 0, the M-by-N solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 ******************************************************************************/
void CORE_ztrsm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t *B, int LDB)
{
    cblas_ztrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB);
}
