/**
 *
 * @file core_zgeqrt.c
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
#include <lapacke.h>
#include "parsec/parsec_config.h"
#include "dplasma.h"
#include "cores/dplasma_cores.h"
#include "dplasma_zcores.h"

/***************************************************************************//**
 *
 * @ingroup CORE_parsec_complex64_t
 *
 *  CORE_zgeqrt computes a QR factorization of a complex M-by-N tile A:
 *  A = Q * R.
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k), where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
 *  and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the min(M,N)-by-N upper trapezoidal tile R (R is
 *         upper triangular if M >= N); the elements below the diagonal,
 *         with the array TAU, represent the unitary tile Q as a
 *         product of elementary reflectors (see Further Details).
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int CORE_zgeqrt(int M, int N, int IB,
                parsec_complex64_t *A, int LDA,
                parsec_complex64_t *T, int LDT,
                parsec_complex64_t *TAU,
                parsec_complex64_t *WORK)
{
    int i, k, sb;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ((IB < 0) || ( (IB == 0) && ((M > 0) && (N > 0)) )) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDA < coreblas_imax(1,M)) && (M > 0)) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ((LDT < coreblas_imax(1,IB)) && (IB > 0)) {
        coreblas_error(7, "Illegal value of LDT");
        return -7;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    k = coreblas_imin(M, N);

    for(i = 0; i < k; i += IB) {
        sb = coreblas_imin(IB, k-i);

        LAPACKE_zgeqr2_work(LAPACK_COL_MAJOR, M-i, sb,
                            &A[LDA*i+i], LDA, &TAU[i], WORK);

        LAPACKE_zlarft_work(LAPACK_COL_MAJOR,
            lapack_const(PlasmaForward),
            lapack_const(PlasmaColumnwise),
            M-i, sb,
            &A[LDA*i+i], LDA, &TAU[i],
            &T[LDT*i], LDT);

        if (N > i+sb) {
            LAPACKE_zlarfb_work(
                LAPACK_COL_MAJOR,
                lapack_const(PlasmaLeft),
                lapack_const(PlasmaConjTrans),
                lapack_const(PlasmaForward),
                lapack_const(PlasmaColumnwise),
                M-i, N-i-sb, sb,
                &A[LDA*i+i],      LDA,
                &T[LDT*i],        LDT,
                &A[LDA*(i+sb)+i], LDA,
                WORK, N-i-sb);
        }
    }
    return PLASMA_SUCCESS;
}
