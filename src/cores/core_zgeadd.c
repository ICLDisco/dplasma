/**
 *
 * @file core_zgeadd.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/**
 ******************************************************************************
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zgeadd adds two matrices together as in PBLAS pzgeadd.
 *
 *       B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          = PlasmaNoTrans:   op(A) = A
 *          = PlasmaTrans:     op(A) = A'
 *          = PlasmaConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices op(A) and B.
 *
 * @param[in] N
 *          Number of columns of the matrices op(A) and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans = PlasmaNoTrans, LDA-by-M
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,k), with k=M, if
 *          trans = PlasmaNoTrans, and k=N otherwise.
 *
 * @param[in] beta
 *          Scalar factor of B.
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *          On exit, B = alpha * op(A) + beta * B
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeadd = PCORE_zgeadd
#define CORE_zgeadd PCORE_zgeadd
#endif
int CORE_zgeadd(PLASMA_enum trans, int M, int N,
                      PLASMA_Complex64_t  alpha,
                const PLASMA_Complex64_t *A, int LDA,
                      PLASMA_Complex64_t  beta,
                      PLASMA_Complex64_t *B, int LDB)
{
    int i, j;

    if ((trans != PlasmaNoTrans) &&
        (trans != PlasmaTrans)   &&
        (trans != PlasmaConjTrans))
    {
        coreblas_error(1, "illegal value of trans");
        return -1;
    }

    if (M < 0) {
        coreblas_error(2, "Illegal value of M");
        return -2;
    }
    if (N < 0) {
        coreblas_error(3, "Illegal value of N");
        return -3;
    }
    if ( ((trans == PlasmaNoTrans) && (LDA < max(1,M)) && (M > 0)) ||
         ((trans != PlasmaNoTrans) && (LDA < max(1,N)) && (N > 0)) )
    {
        coreblas_error(6, "Illegal value of LDA");
        return -6;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        coreblas_error(8, "Illegal value of LDB");
        return -8;
    }

    switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PlasmaConjTrans:
        for (j=0; j<N; j++, A++) {
            for(i=0; i<M; i++, B++) {
                *B = beta * (*B) + alpha * conj(A[LDA*i]);
            }
            B += LDB-M;
        }
        break;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    case PlasmaTrans:
        for (j=0; j<N; j++, A++) {
            for(i=0; i<M; i++, B++) {
                *B = beta * (*B) + alpha * A[LDA*i];
            }
            B += LDB-M;
        }
        break;

    case PlasmaNoTrans:
    default:
        for (j=0; j<N; j++) {
            for(i=0; i<M; i++, B++, A++) {
                *B = beta * (*B) + alpha * (*A);
            }
            A += LDA-M;
            B += LDB-M;
        }
    }
    return PLASMA_SUCCESS;
}
