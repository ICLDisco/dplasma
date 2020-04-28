/**
 *
 * @file core_zpltmg_chebvand.c
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
#include <lapacke.h>
#include <math.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zpltmg_chebvand is a kernel used in Vandermonde-like matrix generation
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999859
 *
 *  Vandermonde-like matrix for the Chebyshev polynomials
 *
 *  Produces the (primal) Chebyshev Vandermonde matrix based on the vector of
 *  points p, which define where the Chebyshev polynomial is calculated.
 *
 *  If seed != 0, C(i,j) = Ti – 1(p(j)) where Ti – 1 is the Chebyshev
 *  polynomial of degree i – 1, and p is a vector of N equally spaced points on
 *  the interval [0,1].
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A to initialize. M >= 2.
 *
 * @param[in] N
 *         The number of columns of the tile A to initialize. N >= 0.
 *
 * @param[out] A
 *         On entry, the M-by-N tile to be initialized.
 *         On exit, each element of A is defined by:
 *         A(i,j) = Ti – 1(p(j)) where Ti – 1 is the Chebyshev polynomial of
 *         degree i – 1
 *
 * @param[in] LDA
 *         The leading dimension of the tile A. LDA >= max(1,M).
 *
 * @param[in] gN
 *         The global number of columns of the full matrix, A is belonging to. gN >= (n0+gN).
 *
 * @param[in] m0
 *         The index of the first row of tile A in the full matrix. m0 >= 0.
 *
 * @param[in] n0
 *         The index of the first column of tile A in the full matrix. n0 >= 0.
 *
 * @param[in] gN
 *         The global number of columns of the full matrix, A is belonging to. gN >= (n0+gN).
 *
 * @param[in,out] W
 *          Workspace of size 2-by-N, that contains the N triplets:
 *               ( A( m0-2, j), A(m0-1, j) )
 *          On entry, if m == 0, W is uinitialized, otherwise contains the data
 *          described above.
 *          On exit, contains the triplets ( A(m0+M-2, j), A(m0+M-1, j) )
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpltmg_chebvand = PCORE_zpltmg_chebvand
#define CORE_zpltmg_chebvand PCORE_zpltmg_chebvand
#endif
int CORE_zpltmg_chebvand( int M, int N, PLASMA_Complex64_t *A, int LDA,
                          int gN, int m0, int n0,
                          PLASMA_Complex64_t *W )
{
    PLASMA_Complex64_t step;
    int i, j, jj;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ((LDA < max(1,M)) && (M > 0)) {
        coreblas_error(4, "Illegal value of LDA");
        return -4;
    }
    if (m0 < 0) {
        coreblas_error(6, "Illegal value of m0");
        return -6;
    }
    if (n0 < 0) {
        coreblas_error(7, "Illegal value of n0");
        return -7;
    }
    if (gN < n0+N) {
        coreblas_error(5, "Illegal value of gN");
        return -5;
    }

    step = (PLASMA_Complex64_t)1. / (gN - 1.);

    /* Initialize W if required */
    if (m0 == 0) {
        for (j=0, jj=n0; j<N; j++, jj++) {
            W[2*j  ] = 1.;
            W[2*j+1] = jj * step;
        }

        if ( M == 1 ) {
            LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A',
                                 1, N, W, 2, A, LDA );
            return PLASMA_SUCCESS;
        }

        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A',
                             2, N, W, 2, A, LDA );

        M -= 2;
        A += 2;
    }

    /* Case NB=1, W contains row 0 and 1 and M should be 1 */
    if (m0 == 1) {
        if (M != 1) {
            coreblas_error(1, "Illegal value of M for m0 = 1");
            return -1;
        }

        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A',
                             1, N, W+1, 2, A, LDA );
        return PLASMA_SUCCESS;
    }

    for (j=0, jj=n0; j<N; j++, jj++) {
        /* First line */
        if (M > 0) {
            A[LDA*j] = 2. * jj * step * W[j*2 + 1] - W[j*2  ];
        }

        /* Second line */
        if (M > 1) {
            A[LDA*j+1] = 2. * jj * step * A[LDA*j  ] - W[j*2+1];
        }

        for (i=2; i<M; i++) {
            A[LDA*j+i] = 2. * jj * step * A[LDA*j+i-1] - A[LDA*j+i-2];
        }
    }

    if ( M == 1 ) {
        cblas_zcopy( N, W+1, 2, W, 2 );
        cblas_zcopy( N, A+M-1, LDA, W+1, 2 );
    }
    else {
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A',
                             2, N, A + M - 2, LDA, W, 2 );
    }

    return PLASMA_SUCCESS;
}
