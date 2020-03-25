/**
 *
 * @file core_zpltmg_toeppd.c
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
#include "parsec/parsec_config.h"
#include "dplasma.h"
#include "dplasma_cores.h"
#include "dplasma_zcores.h"
#include "core_zblas.h"

#undef REAL
#define COMPLEX

#define pi (3.1415926535897932384626433832795028841971693992)

/***************************************************************************/
/**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zpltmg_toeppd1 is the first kernel used in toeppd matrix generation.
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000272
 *
 *  A toeppd matrix is an n-by-n symmetric, positive semi-definite (SPD)
 *  Toeplitz matrix composed of the sum of m rank 2 (or, for certain theta,
 *  rank 1) SPD Toeplitz matrices. Specifically,
 *
 *  T = w(1)*T(theta(1)) + ... + w(m)*T(theta(m))
 *
 *  where T(theta(k)) has (i,j) element cos(2*pi*theta(k)*(i-j)).
 *
 *  In this matrix generation: w = rand(m,1), and theta = rand(m,1).
 *
 *  This kernel generates a portion of size 2-by-m of the full W and theta
 *  vectors.
 *
 *******************************************************************************
 *
 * @param[in] gM
 *         The size of the full vectors W and theta. gM >= M+m0.
 *
 * @param[in] m0
 *         Index of the first element of W, in the full vector. m0 >= 0
 *
 * @param[in] M
 *         The number of elements to generate for w and theta vector. M >= 0.
 *
 * @param[out] W
 *          An 2-by-M matrix.
 *          On exit, the first row contains the walue of w[m0;m0+M]
 *          The second row contains the vector 2*pi*theta[m0;m0+M]
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all call to this routines generating the w and theta vectors.
 *
 ******************************************************************************/
void CORE_zpltmg_toeppd1( int gM, int m0, int M, parsec_complex64_t *W,
                          unsigned long long int seed )
{
    int i;

    /* Store the transpose to have the couple of values used together */
    CORE_zplrnt( 2, M, W, 2, gM, 0, m0, seed );

    for(i=0; i<M; i++, W+=2) {
        /* Translate random value to interval [0,1] */
        W[0] += 0.5;

        /* Initialize second value of each couple to 2 * Pi * (rand() + 0.5) */
#ifdef COMPLEX
        W[1] = 2. * pi * ( W[1] + 0.5 + I * 0.5 );
#else
        W[1] = 2. * pi * ( W[1] + 0.5 );
#endif
    }
}



/***************************************************************************/
/**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zpltmg_toeppd2 is the first kernel used in toeppd matrix generation.
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000272
 *
 *  A toeppd matrix is an n-by-n symmetric, positive semi-definite (SPD)
 *  Toeplitz matrix composed of the sum of m rank 2 (or, for certain theta,
 *  rank 1) SPD Toeplitz matrices. Specifically,
 *
 *  T = w(1)*T(theta(1)) + ... + w(m)*T(theta(m))
 *
 *  where T(theta(k)) has (i,j) element cos(2*pi*theta(k)*(i-j)).
 *
 *  In this matrix generation: w = rand(m,1), and theta = rand(m,1).
 *
 *  This kernel adds to the tile A the local sum of:
 *        w(1)*T(theta(1)) + ... + w(K) * T(theta(K))
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A. M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A. N >= 0.
 *
 * @param[in] K
 *         The number of matrices W() * T(theta()) to apply.
 *
 * @param[in] m0
 *         The index of the first row of tile A in the full matrix. m0 >= 0.
 *
 * @param[in] n0
 *         The index of the first column of tile A in the full matrix. n0 >= 0.
 *
 * @param[in] W
 *          The 2-by-K array that stores the values of W and 2*pi*Theta.
 *          W being stored on the first row, 2*pi*theta on the second.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be initialized with a partial sum of the
 *         Toeppliz matrices.
 *         On exit, the M-by-N tile update with the sum of K extra Toepplitz matrices
 *
 * @param[in] LDA
 *         The leading dimension of the tile A. LDA >= max(1,M).
 *
 ******************************************************************************/
void CORE_zpltmg_toeppd2( int M, int N, int K, int m0, int n0,
                          const parsec_complex64_t *W,
                                parsec_complex64_t *A, int LDA )
{
    const parsec_complex64_t *tmpW;
    int i, j, k, ii, jj;

    for (j=0, jj=n0; j<N; j++, jj++) {
        for (i=0, ii=m0; i<M; i++, ii++, A++) {
            tmpW = W;
            for(k=0; k<K; k++, tmpW+=2) {
                *A += creal(tmpW[0]) * ccos( tmpW[1] * (double)(ii-jj) );
            }
        }
        A += (LDA - M);
    }
}
