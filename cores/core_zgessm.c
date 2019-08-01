/**
 *
 * @file core_zgessm.c
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
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zgessm applies the factors L computed by CORE_zgetrf_incpiv to
 *  a complex M-by-N tile A.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] K
 *         The number of columns of the tile L. K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_zgetrf_incpiv.
 *
 * @param[in] L
 *         The M-by-K lower triangular tile.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,M).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, updated by the application of L.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
int CORE_zgessm(int M, int N, int K, int IB,
                const int *IPIV,
                const parsec_complex64_t *L, int LDL,
                parsec_complex64_t *A, int LDA)
{
    static parsec_complex64_t zone  =  1.0;
    static parsec_complex64_t mzone = -1.0;
    static int                ione  =  1;

    int i, sb;
    int tmp, tmp2;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (K < 0) {
        coreblas_error(3, "Illegal value of K");
        return -3;
    }
    if (IB < 0) {
        coreblas_error(4, "Illegal value of IB");
        return -4;
    }
    if ((LDL < coreblas_imax(1,M)) && (M > 0)) {
        coreblas_error(7, "Illegal value of LDL");
        return -7;
    }
    if ((LDA < coreblas_imax(1,M)) && (M > 0)) {
        coreblas_error(9, "Illegal value of LDA");
        return -9;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (K == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    for(i = 0; i < K; i += IB) {
        sb = coreblas_imin(IB, K-i);
        /*
         * Apply interchanges to columns I*IB+1:IB*( I+1 )+1.
         */
        tmp  = i+1;
        tmp2 = i+sb;
        LAPACKE_zlaswp_work(LAPACK_COL_MAJOR, N, A, LDA, tmp, tmp2, IPIV, ione);
        /*
         * Compute block row of U.
         */
        cblas_ztrsm(
            CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
            sb, N, CBLAS_SADDR(zone),
            &L[LDL*i+i], LDL,
            &A[i], LDA );

        if (i+sb < M) {
        /*
        * Update trailing submatrix.
        */
        cblas_zgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans,
            M-(i+sb), N, sb,
            CBLAS_SADDR(mzone), &L[LDL*i+(i+sb)], LDL,
            &A[i], LDA,
            CBLAS_SADDR(zone), &A[i+sb], LDA );
        }
    }
    return PLASMA_SUCCESS;
}
