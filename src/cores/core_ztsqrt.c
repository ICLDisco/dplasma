/**
 *
 * @file core_ztsqrt.c
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
#include "parsec/parsec_config.h"
#include "dplasma.h"
#include "dplasma_zcores.h"
#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 * CORE_ztsqrt computes a QR factorization of a rectangular matrix
 * formed by coupling a complex N-by-N upper triangular tile A1
 * on top of a complex M-by-N tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of columns of the tile A2. M >= 0.
 *
 * @param[in] N
 *         The number of rows of the tile A1.
 *         The number of columns of the tiles A1 and A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the N-by-N tile A1.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the N-by-N upper trapezoidal tile R;
 *         the elements below the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N tile A2.
 *         On exit, all the elements with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M).
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

int CORE_ztsqrt(int M, int N, int IB,
                PLASMA_Complex64_t *A1, int LDA1,
                PLASMA_Complex64_t *A2, int LDA2,
                PLASMA_Complex64_t *T, int LDT,
                PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK)
{
    static PLASMA_Complex64_t zone  = 1.0;
    static PLASMA_Complex64_t zzero = 0.0;

    PLASMA_Complex64_t alpha;
    int i, ii, sb;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (IB < 0) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDA2 < coreblas_imax(1,M)) && (M > 0)) {
        coreblas_error(8, "Illegal value of LDA2");
        return -8;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    for(ii = 0; ii < N; ii += IB) {
        sb = coreblas_imin(N-ii, IB);
        for(i = 0; i < sb; i++) {
            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate
             * A( II*IB+I:M, II*IB+I )
             */
            LAPACKE_zlarfg_work(M+1, &A1[LDA1*(ii+i)+ii+i], &A2[LDA2*(ii+i)], 1, &TAU[ii+i]);

            if (ii+i+1 < N) {
                /*
                 * Apply H( II*IB+I ) to A( II*IB+I:M, II*IB+I+1:II*IB+IB ) from the left
                 */
                alpha = -conj(TAU[ii+i]);
                cblas_zcopy(
                    sb-i-1,
                    &A1[LDA1*(ii+i+1)+(ii+i)], LDA1,
                    WORK, 1);
#ifdef COMPLEX
                LAPACKE_zlacgv_work(sb-i-1, WORK, 1);
#endif
                cblas_zgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaConjTrans,
                    M, sb-i-1,
                    CBLAS_SADDR(zone), &A2[LDA2*(ii+i+1)], LDA2,
                    &A2[LDA2*(ii+i)], 1,
                    CBLAS_SADDR(zone), WORK, 1);
#ifdef COMPLEX
                LAPACKE_zlacgv_work(sb-i-1, WORK, 1 );
#endif
                cblas_zaxpy(
                    sb-i-1, CBLAS_SADDR(alpha),
                    WORK, 1,
                    &A1[LDA1*(ii+i+1)+ii+i], LDA1);
#ifdef COMPLEX
                LAPACKE_zlacgv_work(sb-i-1, WORK, 1 );
#endif
                cblas_zgerc(
                    CblasColMajor, M, sb-i-1, CBLAS_SADDR(alpha),
                    &A2[LDA2*(ii+i)], 1,
                    WORK, 1,
                    &A2[LDA2*(ii+i+1)], LDA2);
            }
            /*
             * Calculate T
             */
            alpha = -TAU[ii+i];
            cblas_zgemv(
                CblasColMajor, (CBLAS_TRANSPOSE)PlasmaConjTrans, M, i,
                CBLAS_SADDR(alpha), &A2[LDA2*ii], LDA2,
                &A2[LDA2*(ii+i)], 1,
                CBLAS_SADDR(zzero), &T[LDT*(ii+i)], 1);

            cblas_ztrmv(
                CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                (CBLAS_TRANSPOSE)PlasmaNoTrans, (CBLAS_DIAG)PlasmaNonUnit, i,
                &T[LDT*ii], LDT,
                &T[LDT*(ii+i)], 1);

            T[LDT*(ii+i)+i] = TAU[ii+i];
        }
        if (N > ii+sb) {
            CORE_ztsmqr(
                PlasmaLeft, PlasmaConjTrans,
                sb, N-(ii+sb), M, N-(ii+sb), IB, IB,
                &A1[LDA1*(ii+sb)+ii], LDA1,
                &A2[LDA2*(ii+sb)], LDA2,
                &A2[LDA2*ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, sb);
        }
    }
    return PLASMA_SUCCESS;
}
