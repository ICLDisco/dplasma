/**
 *
 * @file core_zgemv.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "parsec/parsec_config.h"
#include "dplasma.h"
#include "cores/dplasma_cores.h"
#include "dplasma_zcores.h"

/***************************************************************************//**
 *
 * @ingroup CORE_parsec_complex64_t
 *
 * CORE_zgemv performs one of the matrix-vector operations
 *
 *    y := alpha*A*x    + beta*y,   or
 *    y := alpha*A**T*x + beta*y,   or
 *    y := alpha*A**H*x + beta*y,
 *
 * where alpha and beta are scalars, x and y are vectors,
 * and A is an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg PlasmaNoTrans:   y := alpha*A*x    + beta*y
 *         @arg PlasmaTrans:     y := alpha*A**T*x + beta*y
 *         @arg PlasmaConjTrans: y := alpha*A**H*x + beta*y
 *
 * @param[in] m
 *         Number of rows of matrix A.
 *
 * @param[in] n
 *         Number of columns of matrix A.
 *
 * @param[in] alpha
 *         Scalar alpha.
 *
 * @param[in] A
 *         On entry, m by n matrix A. Dimension (lda,n).
 *
 * @param[in] lda
 *         Leading dimension of array A. lda >= max(1,m).
 *
 * @param[in] x
 *         On entry, vector x.
 *         If trans == PlasmaNoTrans, the n vector x has dimension 1 + (n-1)*abs(incx).
 *         Else, the m vector x has dimension 1 + (m-1)*abs(incx).
 *
 * @param[in] incx
 *         Increment between elements of x. incx must not be zero.
 *
 * @param[in] beta
 *         Scalar beta.
 *
 * @param[in,out] y
 *         On entry, vector y. On exit, y  is overwritten by updated vector y.
 *         If trans == PlasmaNoTrans, the m vector y has dimension 1 + (m-1)*abs(incy).
 *         Else, the n vector y has dimension 1 + (n-1)*abs(incy).
 *
 * @param[in] incy
 *         Increment between elements of y. incy must not be zero.
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemv = PCORE_zgemv
#define CORE_zgemv PCORE_zgemv
#endif
void CORE_zgemv(PLASMA_enum trans,
                int m, int n,
                parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                                          const parsec_complex64_t *x, int incx,
                parsec_complex64_t beta,        parsec_complex64_t *y, int incy)
{
    cblas_zgemv(
        CblasColMajor,
        (CBLAS_TRANSPOSE)trans,
        m, n,
        CBLAS_SADDR(alpha), A, lda,
        x, incx,
        CBLAS_SADDR(beta), y, incy );
}
