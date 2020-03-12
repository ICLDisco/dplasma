/**
 *
 * @file core_ztrasm.c
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
#include <math.h>
#include "core_blas.h"

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_ztrasm - Computes the sums of the absolute values of elements in a same
 *  row or column in a triangular matrix.
 *  This function is an auxiliary function to triangular matrix norm computations.
 *
 *******************************************************************************
 *
 * @param[in] storev
 *          Specifies whether the sums are made per column or row.
 *          = PlasmaColumnwise: Computes the sum on each column
 *          = PlasmaRowwise:    Computes the sum on each row
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular
 *          = PlasmaUpper: Upper triangle of A is referenced;
 *          = PlasmaLower: Lower triangle of A is referenced.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          A is a M-by-N matrix.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,M).
 *
 * @param[in,out] work
 *          Array of dimension M if storev = PlasmaRowwise; N otherwise.
 *          On exit, contains the sums of the absolute values per column or row
 *          added to the input values.
 *
 ******************************************************************************/
void CORE_ztrasm(PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N,
                 const PLASMA_Complex64_t *A, int lda, double *work)
{
    const PLASMA_Complex64_t *tmpA;
    int i, j, imax;
    int idiag = (diag == PlasmaUnit) ? 1 : 0;

    /*
     * PlasmaUpper / PlasmaColumnwise
     */
    if  (uplo == PlasmaUpper ) {
        M = coreblas_imin(M, N);

        if (storev == PlasmaColumnwise) {
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                imax = coreblas_imin(j+1-idiag, M);

                if ( j < M )
                    work[j] += idiag;

                for (i = 0; i < imax; i++) {
                    work[j] += cabs(*tmpA);
                    tmpA++;
                }
            }
        }
        /*
         * PlasmaUpper / PlasmaRowwise
         */
        else {
            if (diag == PlasmaUnit) {
                for (i = 0; i < M; i++) {
                    work[i] += 1.;
                }
            }
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                imax = coreblas_imin(j+1-idiag, M);

                for (i = 0; i < imax; i++) {
                    work[i] += cabs(*tmpA);
                    tmpA++;
                }
            }
        }
    } else {
        N = coreblas_imin(M, N);

        /*
         * PlasmaLower / PlasmaColumnwise
         */
        if (storev == PlasmaColumnwise) {
            for (j = 0; j < N; j++) {
                tmpA = A + j * (lda+1) + idiag;

                work[j] += idiag;
                for (i = j+idiag; i < M; i++) {
                    work[j] += cabs(*tmpA);
                    tmpA++;
                }
            }
        }
        /*
         * PlasmaLower / PlasmaRowwise
         */
        else {
            if (diag == PlasmaUnit) {
                for (i = 0; i < N; i++) {
                    work[i] += 1.;
                }
            }
            for (j = 0; j < N; j++) {
                tmpA = A + j * (lda+1) + idiag;

                for (i = j+idiag; i < M; i++) {
                    work[i] += cabs(*tmpA);
                    tmpA++;
                }
            }
        }
    }
}
