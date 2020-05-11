/**
 *
 * @file core_zpltmg_hankel.c
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

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zpltmg_hankel is a kernel used in Hankel matrix generation
 *
 *  See http://en.wikipedia.org/wiki/Hankel_matrix
 *
 *  Hankel matrix
 *
 *  In linear algebra, a Hankel matrix (or catalecticant matrix), named after
 *  Hermann Hankel, is a square matrix with constant skew-diagonals (positive
 *  sloping diagonals), e.g.:
 *
 *  \f[ \begin{bmatrix}
 *  a & b & c & d & e \\
 *  b & c & d & e & f \\
 *  c & d & e & f & g \\
 *  d & e & f & g & h \\
 *  e & f & g & h & i \\
 *  \end{bmatrix}
 *  \f].
 *
 *  A(i,j) = A(i-1,j+1)
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be initialized.
 *            = PlasmaUpperLower: All the matrix A
 *            = PlasmaUpper: Upper triangular part
 *            = PlasmaLower: Lower triangular part
 *
 * @param[in] M
 *         The number of rows of the tile A to initialize. M >= 2.
 *
 * @param[in] N
 *         The number of columns of the tile A to initialize. N >= 0.
 *
 * @param[out] A
 *         On entry, the M-by-N tile to be initialized.
 *
 * @param[in] LDA
 *         The leading dimension of the tile A. LDA >= max(1,M).
 *
 * @param[in] m0
 *         The index of the first row of tile A in the full matrix. m0 >= 0.
 *
 * @param[in] n0
 *         The index of the first column of tile A in the full matrix. n0 >= 0.
 *
 * @param[in] nb
 *         The size of the V1 and V2 vectors
 *
 * @param[in] V1
 *          Workspace of size nb, that contains the first column of the tile
 *
 * @param[in] V2
 *          Workspace of size nb), that contains the last column of the tile
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpltmg_hankel = PCORE_zpltmg_hankel
#define CORE_zpltmg_hankel PCORE_zpltmg_hankel
#endif
int CORE_zpltmg_hankel( PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t *A, int LDA,
                        int m0, int n0, int nb,
                        const PLASMA_Complex64_t *V1,
                        const PLASMA_Complex64_t *V2 )
{
    int p, i, j, ii, jj;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(2, "Illegal value of M");
        return -2;
    }
    if (N < 0) {
        coreblas_error(3, "Illegal value of N");
        return -3;
    }
    if ((LDA < max(1,M)) && (M > 0)) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if (m0 < 0) {
        coreblas_error(6, "Illegal value of m0");
        return -6;
    }
    if (n0 < 0) {
        coreblas_error(7, "Illegal value of n0");
        return -7;
    }
    if (nb < 0) {
        coreblas_error(8, "Illegal value of nb");
        return -8;
    }

    switch ( uplo ) {
    case PlasmaUpper:
        for (j=0, jj=n0; j<N; j++, jj++) {
            for (i=0, ii=m0; i<M; i++, ii++) {
                if (ii > jj)
                    continue;

                p = i + j;

                if ( p < nb ){
                    A[LDA*j + i] = V1[ p ];
                } else {
                    A[LDA*j + i] = V2[ p%nb ];
                }
            }
        }
        break;

    case PlasmaLower:
        for (j=0, jj=n0; j<N; j++, jj++) {
            for (i=0, ii=m0; i<M; i++, ii++) {
                if (ii < jj)
                    continue;

                p = i + j;

                if ( p < nb ){
                    A[LDA*j + i] = V1[ p ];
                } else {
                    A[LDA*j + i] = V2[ p%nb ];
                }
            }
        }
        break;

    default:
        for (j=0, jj=n0; j<N; j++, jj++) {
            for (i=0, ii=m0; i<M; i++, ii++) {
                p = i + j;

                if ( p < nb ){
                    A[LDA*j + i] = V1[ p ];
                } else {
                    A[LDA*j + i] = V2[ p%nb ];
                }
            }
        }
    }

    return PLASMA_SUCCESS;
}
