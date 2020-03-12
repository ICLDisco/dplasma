/**
 *
 * @file core_dzasum.c
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
#include "dplasma_dcores.h"

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_dzasum - Computes the sums of the absolute values of elements in a same
 *  row or column.
 *  This function is an auxiliary function to norm computations.
 *
 *******************************************************************************
 *
 * @param[in] storev
 *          Specifies whether the sums are made per column or row.
 *          = PlasmaColumnwise: Computes the sum on each column
 *          = PlasmaRowwise:    Computes the sum on each row
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular or general
 *          = PlasmaUpperLower: All matrix A is referenced;
 *          = PlasmaUpper: Upper triangle of A is referenced;
 *          = PlasmaLower: Lower triangle of A is referenced.
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
 * @param[out] work
 *          Array of dimension M if storev = PlasmaRowwise; N otherwise.
 *          On exit, contains the sums of the absolute values per column or row.
 *
 ******************************************************************************/
void CORE_dzasum(PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                 const PLASMA_Complex64_t *A, int lda, double *work)
{
    const PLASMA_Complex64_t *tmpA;
    double *tmpW, sum, abs;
    int i,j;

    switch (uplo) {
    case PlasmaUpper:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda);
            sum = 0.0;
            for (i = 0; i < j; i++) {
                abs      = cabs(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum + cabs(*tmpA);
        }
        break;
    case PlasmaLower:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda)+j;

            sum = 0.0;
            work[j] += cabs(*tmpA);

            tmpA++;
            for (i = j+1; i < M; i++) {
                abs      = cabs(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum;
        }
        break;
    case PlasmaUpperLower:
    default:
        if (storev == PlasmaColumnwise) {
            for (j = 0; j < N; j++) {
                /* work[j] += cblas_dzasum(M, &(A[j*lda]), 1); */
                tmpA = A+(j*lda);
                for (i = 0; i < M; i++) {
                    work[j] +=  cabs(*tmpA);
                    tmpA++;
                }
            }
        }
        else {
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                tmpW = work;
                for (i = 0; i < M; i++) {
                    /* work[i] += cabs( A[j*lda+i] );*/
                    *tmpW += cabs( *tmpA );
                    tmpA++; tmpW++;
                }
            }
        }
    }
}
