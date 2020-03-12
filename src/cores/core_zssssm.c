/**
 *
 * @file core_zssssm.c
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
 *  CORE_zssssm applies the LU factorization update from a complex
 *  matrix formed by a lower triangular IB-by-K tile L1 on top of a
 *  M2-by-K tile L2 to a second complex matrix formed by a M1-by-N1
 *  tile A1 on top of a M2-by-N2 tile A2 (N1 == N2).
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M1
 *         The number of rows of the tile A1.  M1 >= 0.
 *
 * @param[in] N1
 *         The number of columns of the tile A1.  N1 >= 0.
 *
 * @param[in] M2
 *         The number of rows of the tile A2 and of the tile L2.
 *         M2 >= 0.
 *
 * @param[in] N2
 *         The number of columns of the tile A2.  N2 >= 0.
 *
 * @param[in] K
 *         The number of columns of the tiles L1 and L2.  K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is updated by the application of L (L1 L2).
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is updated by the application of L (L1 L2).
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M2).
 *
 * @param[in] L1
 *         The IB-by-K lower triangular tile as returned by
 *         CORE_ztstrf.
 *
 * @param[in] LDL1
 *         The leading dimension of the array L1.  LDL1 >= max(1,IB).
 *
 * @param[in] L2
 *         The M2-by-K tile as returned by CORE_ztstrf.
 *
 * @param[in] LDL2
 *         The leading dimension of the array L2.  LDL2 >= max(1,M2).
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_ztstrf.
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/

int CORE_zssssm(int M1, int N1, int M2, int N2, int K, int IB,
                PLASMA_Complex64_t *A1, int LDA1,
                PLASMA_Complex64_t *A2, int LDA2,
                const PLASMA_Complex64_t *L1, int LDL1,
                const PLASMA_Complex64_t *L2, int LDL2,
                const int *IPIV)
{
    static PLASMA_Complex64_t zone  = 1.0;
    static PLASMA_Complex64_t mzone =-1.0;

    int i, ii, sb;
    int im, ip;

    /* Check input arguments */
    if (M1 < 0) {
        coreblas_error(1, "Illegal value of M1");
        return -1;
    }
    if (N1 < 0) {
        coreblas_error(2, "Illegal value of N1");
        return -2;
    }
    if (M2 < 0) {
        coreblas_error(3, "Illegal value of M2");
        return -3;
    }
    if (N2 < 0) {
        coreblas_error(4, "Illegal value of N2");
        return -4;
    }
    if (K < 0) {
        coreblas_error(5, "Illegal value of K");
        return -5;
    }
    if (IB < 0) {
        coreblas_error(6, "Illegal value of IB");
        return -6;
    }
    if (LDA1 < coreblas_imax(1,M1)) {
        coreblas_error(8, "Illegal value of LDA1");
        return -8;
    }
    if (LDA2 < coreblas_imax(1,M2)) {
        coreblas_error(10, "Illegal value of LDA2");
        return -10;
    }
    if (LDL1 < coreblas_imax(1,IB)) {
        coreblas_error(12, "Illegal value of LDL1");
        return -12;
    }
    if (LDL2 < coreblas_imax(1,M2)) {
        coreblas_error(14, "Illegal value of LDL2");
        return -14;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    ip = 0;

    for(ii = 0; ii < K; ii += IB) {
        sb = coreblas_imin(K-ii, IB);

        for(i = 0; i < sb; i++) {
            im = IPIV[ip]-1;

            if (im != (ii+i)) {
                im = im - M1;
                cblas_zswap(N1, &A1[ii+i], LDA1, &A2[im], LDA2);
            }
            ip = ip + 1;
        }

        cblas_ztrsm(
            CblasColMajor, CblasLeft, CblasLower,
            CblasNoTrans, CblasUnit,
            sb, N1, CBLAS_SADDR(zone),
            &L1[LDL1*ii], LDL1,
            &A1[ii], LDA1);

        cblas_zgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans,
            M2, N2, sb,
            CBLAS_SADDR(mzone), &L2[LDL2*ii], LDL2,
            &A1[ii], LDA1,
            CBLAS_SADDR(zone), A2, LDA2);
    }
    return PLASMA_SUCCESS;
}
