/**
 *
 * @file core_zlaset.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "core_blas.h"

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zlaset - Sets the elements of the matrix A on the diagonal
 *  to beta and on the off-diagonals to alpha
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = PlasmaUpper: Upper part of A is set;
 *          = PlasmaLower: Lower part of A is set;
 *          = PlasmaUpperLower: ALL elements of A are set.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the off-diagonal elements are to be set.
 *
 * @param[in] beta
 *         The constant to which the diagonal elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaset = PCORE_zlaset
#define CORE_zlaset PCORE_zlaset
#endif
void CORE_zlaset(PLASMA_enum uplo, int M, int N,
                 PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta,
                 PLASMA_Complex64_t *A, int LDA)
{
    LAPACKE_zlaset_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, alpha, beta, A, LDA);
}

