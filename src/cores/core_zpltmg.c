/**
 *
 * @file core_zpltmg.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Ichitaro Yamazaki
 * @author Julien Herrmann
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/

#include <math.h>
#include <lapacke.h>
#include "common.h"

#define pi (3.1415926535897932384626433832795028841971693992)

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zpltmg initialize a tile of a random matrix from the MatLab
 *  gallery configured with the default parameters, and a few other
 *  specific matrices.
 *
 *******************************************************************************
 *
 * @param[in] mtxtype
 *          Possible types are: PlasmaMatrixRandom, PlasmaMatrixHadamard,
 *          PlasmaMatrixParter, PlasmaMatrixRis, PlasmaMatrixKms,
 *          PlasmaMatrixMoler, PlasmaMatrixCompan, PlasmaMatrixRiemann,
 *          PlasmaMatrixLehmer, PlasmaMatrixMinij, PlasmaMatrixDorr,
 *          PlasmaMatrixDemmel, PlasmaMatrixInvhess, PlasmaMatrixCauchy,
 *          PlasmaMatrixHilb, PlasmaMatrixLotkin, PlasmaMatrixOrthog,
 *          PlasmaMatrixWilkinson, PlasmaMatrixFoster, PlasmaMatrixWright,
 *          PlasmaMatrixLangou
 *          (See further in the code for more details)
 *
 * @param[in] M
 *         The number of rows of the tile A. M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A. N >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be initialized.
 *         On exit, the tile initialized in the mtxtype format.
 *
 * @param[in] LDA
 *         The leading dimension of the tile A. LDA >= max(1,M).
 *
 * @param[in] gM
 *         The global number of rows of the full matrix, A is belonging to. gM >= (m0+M).
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
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpltmg = PCORE_zpltmg
#define CORE_zpltmg PCORE_zpltmg
#define CORE_zplrnt  PCORE_zplrnt
void
CORE_zplrnt( int M, int N, PLASMA_Complex64_t *A, int LDA,
             int gM, int m0, int n0,
             unsigned long long int seed );
#endif
int CORE_zpltmg( PLASMA_enum mtxtype,
                  int M, int N, PLASMA_Complex64_t *A, int LDA,
                  int gM, int gN, int m0, int n0,
                  unsigned long long int seed )
{
    int i, j;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(2, "Illegal value of M");
        return -3;
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
        coreblas_error(8, "Illegal value of m0");
        return -8;
    }
    if (n0 < 0) {
        coreblas_error(9, "Illegal value of n0");
        return -9;
    }
    if (gM < m0+M) {
        coreblas_error(6, "Illegal value of gM");
        return -6;
    }
    if (gN < n0+N) {
        coreblas_error(7, "Illegal value of gN");
        return -7;
    }

    /* Quick return */
    if ((M == 0) || (N == 0))
        return PLASMA_SUCCESS;

    switch( mtxtype ) {
    case PlasmaMatrixRandom:
    {
        CORE_zplrnt( M, N, A, LDA, gM, m0, n0, seed );
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/hadamard.html
     *
     * Initialize the tile A to create the Hadamard matrix of order gN.
     *
     * Hadamard matrices are matrices of 1's and -1's whose columns are orthogonal,
     *
     *   H'*H = gN*I
     *
     *   where [gN gN]=size(H) and I = eye(gN,gN) ,.
     *
     *  They have applications in several different areas, including
     *  combinatorics, signal processing, and numerical analysis.
     *
     *  An n-by-n Hadamard matrix with n > 2 exists only if rem(n,4) =
     *  0. This function handles only the cases where n is a power of
     *  2.
     */
    case PlasmaMatrixHadamard:
    {
        int tmp, nbone;

        /* Extra parameters check */
        if (gM != gN) {
            coreblas_error(6, "Illegal value of gM (Matrix must be square)");
            return -6;
        }

        tmp = gM;
        while ( tmp > 1 ) {
            if( tmp % 2 != 0 ) {
                coreblas_error(6, "Illegal value of gM (Matrix dimension must be a power of 2)");
                return -6;
            }
            tmp /= 2;
        }

        for (j=0; j<N; j++) {
            for (i=0; i<M; i++) {
                tmp = ((m0 + i) & (n0 + j));
                nbone = 0;
                while ( tmp != 0 )
                {
                    nbone += ( tmp & 1 );
                    tmp >>= 1;
                }
                A[j*LDA+i] = (PLASMA_Complex64_t)(1. - 2. * ( nbone % 2 ));
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000116
     *
     * Toeplitz matrix with singular values near pi.
     * Returns the tile A, such that the elment of the matrix are 1/(i-j+0.5).
     *
     * C is a Cauchy matrix and a Toeplitz matrix. Most of the
     * singular values of C are very close to pi.
     *
     */
    case PlasmaMatrixParter:
    {
        PLASMA_Complex64_t tmp;

        if (gM != gN) {
            coreblas_error(6, "Illegal value of gM (Matrix must be square for Parter)");
            return -6;
        }

        tmp = (PLASMA_Complex64_t)( .5 + m0 - n0 );
        for (j=0; j<N; j++) {
            for (i=0; i<M; i++) {
                A[j*LDA+i] = (PLASMA_Complex64_t)1. / (PLASMA_Complex64_t)( tmp + i - j );
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000243
     *
     * Symmetric Hankel matrix
     * Returns a symmetric gN-by-gN Hankel matrix with elements.
     *
     *     A(i,j) = 0.5/(n-i-j+1.5)
     *
     * The eigenvalues of A cluster around π/2 and –π/2. This matrix
     * was invented by F.N. Ris.
     *
     */
     case PlasmaMatrixRis:
    {
        PLASMA_Complex64_t tmp;

        if (gM != gN) {
            coreblas_error(6, "Illegal value of gM (Matrix must be square for RIS)");
            return -6;
        }

        tmp = (PLASMA_Complex64_t)( gM - m0 - n0 - 0.5 );
        for (j=0; j<N; j++) {
            for (i=0; i<M; i++) {
                A[j*LDA+i] = (PLASMA_Complex64_t).5 / (PLASMA_Complex64_t)( tmp - i - j );
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000026
     *
     * Kac-Murdock-Szego Toeplitz matrix
     *
     * Returns the n-by-n Kac-Murdock-Szego Toeplitz matrix such that
     * A(i,j) = rho^(abs(i-j)), for real rho.
     *
     * For complex rho, the same formula holds except that elements
     * below the diagonal are conjugated. rho defaults to 0.5.
     *
     * The KMS matrix A has these properties:
     *     - An LDL' factorization with L inv(gallery('triw',n,-rho,1))',
     *       and D(i,i) (1-abs(rho)^2)*eye(n), except D(1,1) = 1.
     *     - Positive definite if and only if 0 < abs(rho) < 1.
     *     - The inverse inv(A) is tridiagonal.Symmetric Hankel matrix
     *
     * In this function, rho is set to 0.5 and cannot be changed.
     */
    case PlasmaMatrixKms:
    {
        PLASMA_Complex64_t rho;

        if (gM != gN) {
            coreblas_error(6, "Illegal value of gM (Matrix must be square for KMS)");
            return -6;
        }

        rho = .5;
        for (j=0; j<N; j++) {
            for (i=0; i<M; i++) {
                A[j*LDA+i] = (PLASMA_Complex64_t)( cpow( rho, fabs( (double)( m0 + i - n0 - j ) ) ) );
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000074
     *
     * Symmetric positive definite matrix
     *
     * Returns the symmetric positive definite n-by-n matrix U'*U,
     * where U = gallery('triw',n,alpha).
     *
     * For the default alpha = -1, A(i,j) = min(i,j)-2, and A(i,i) =
     * i. One of the eigenvalues of A is small.
     *
     */
    case PlasmaMatrixMoler:
    {
        int ii, jj;
        for (j=0,jj=n0; j<N; j++,jj++) {
            for (i=0,ii=m0; i<M; i++,ii++) {
                if ( ii == jj ) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)( ii + 1. );
                } else {
                    A[j*LDA+i] = (PLASMA_Complex64_t)( min( ii, jj ) - 1. );
                }
            }
        }
    }
    break;

    /*
     *  See http://www.mathworks.fr/fr/help/matlab/ref/compan.html
     *
     *  Companion matrix
     *
     *  A = compan(u) returns the corresponding companion matrix whose first row is
     *  -u(2:n)/u(1), where u is a vector of polynomial coefficients. The
     *  eigenvalues of compan(u) are the roots of the polynomial.
     *
     */
    case PlasmaMatrixCompan:
    {
        int jj;

        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', M, N, 0., 0., A, LDA);

        /* First row */
        if ( m0 == 0 ) {
            PLASMA_Complex64_t v0;
            /* Get V0 */
            CORE_zplrnt( 1, 1, &v0, 1, 1, 1, 0, seed );
            v0 = 1. / v0 ;

            /* Initialize random vector */
            CORE_zplrnt( 1, N, A, LDA, 1, 1, n0, seed );
            cblas_zscal( N, CBLAS_SADDR(v0), A, LDA);

            /* Restore A(0,0) */
            if ( n0 == 0 )
                A[0] = 0.;
        }

        /* Sub diagonal*/
        for (j=0,jj=n0; j<N; j++,jj++)
        {
            i = jj + 1 - m0;
            if ( ( i > 0 ) && (i < M) )
            {
                A[j*LDA+i] = 1.;
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000232
     *
     * Matrix associated with the Riemann hypothesis
     *
     * Returns an n-by-n matrix for which the Riemann hypothesis is
     * true if and only if for every eps > 0.
     *
     * The Riemann matrix is defined by:
     *
     *    A = B(2:n+1,2:n+1)
     *
     *    where B(i,j) = i-1 if i divides j, and B(i,j) = -1 otherwise.
     *
     * The Riemann matrix has these properties:
     *     - Each eigenvalue e(i) satisfies abs(e(i)) <= m-1/m, where m = n+1.
     *     - i <= e(i) <= i+1 with at most m-sqrt(m) exceptions.
     *     - All integers in the interval (m/3, m/2] are eigenvalues.
     *
     */
    case PlasmaMatrixRiemann:
    {
        int ii, jj;
        for (j=0,jj=n0+2; j<N; j++,jj++) {
            for (i=0,ii=m0+2; i<M; i++,ii++) {
                if ( jj%ii == 0 ) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)( ii - 1. );
                } else {
                    A[j*LDA+i] = (PLASMA_Complex64_t)( -1. );
                }
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000049
     *
     * Symmetric positive definite matrix
     *
     * Returns the symmetric positive definite n-by-n matrix such that
     * A(i,j) = i/j for j >= i.
     *
     * The Lehmer matrix A has these properties:
     *     - A is totally nonnegative.
     *     - The inverse inv(A) is tridiagonal and explicitly known.
     *     - The order n <= cond(A) <= 4*n*n.Matrix associated with the
     *       Riemann hypothesis
     *
     */
     case PlasmaMatrixLehmer:
    {
        int ii, jj;
        for (j=0,jj=n0+1; j<N; j++,jj++) {
            for (i=0,ii=m0+1; i<M; i++,ii++) {
                if ( jj >= ii ) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)( ii ) / (PLASMA_Complex64_t)( jj );
                } else {
                    A[j*LDA+i] = (PLASMA_Complex64_t)( jj ) / (PLASMA_Complex64_t)( ii );
                }
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000066
     *
     * Symmetric positive definite matrix
     *
     * Returns the n-by-n symmetric positive definite matrix with
     * A(i,j) = min(i,j).
     *
     * The minij matrix has these properties:
     *     - The inverse inv(A) is tridiagonal and equal to -1 times the
     *       second difference matrix, except its (n,n) element is 1.
     *     - Givens' matrix, 2*A-ones(size(A)), has tridiagonal inverse
     *       and eigenvalues 0.5*sec((2*r-1)*pi/(4*n))^2, where r=1:n.
     *     - (n+1)*ones(size(A))-A has elements that are max(i,j) and a
     *       tridiagonal inverse.
     *
     */
    case PlasmaMatrixMinij:
    {
        int ii, jj;
        for (j=0,jj=n0+1; j<N; j++,jj++) {
            for (i=0,ii=m0+1; i<M; i++,ii++) {
                A[j*LDA+i] = (PLASMA_Complex64_t) min( ii, jj );
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999936
     *
     * Diagonally dominant, ill-conditioned, tridiagonal matrix
     *
     * Returns the n-by-n matrix, row diagonally dominant, tridiagonal
     * matrix that is ill-conditioned for small nonnegative values of
     * theta. The default value of theta is 0.01. The Dorr matrix
     * itself is the same as gallery('tridiag',c,d,e).
     *
     */
    case PlasmaMatrixDorr:
    {
        PLASMA_Complex64_t theta = 0.01;
        PLASMA_Complex64_t h     = 1. / ( gN + 1. );
        PLASMA_Complex64_t term  = theta / ( h * h );
        int jj;
        int half = (gN+1) / 2;

        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', M, N, 0., 0., A, LDA);

        /* First half */
        for (j=0, jj=n0; (j<N) && (jj<half); j++, jj++) {
            i = jj - m0;
            if ( ( i < -1 ) || (i > M) )
                continue;

            /* Over the diagonal */
            if (i > 0)
                A[j*LDA + i-1] = - term - (0.5 - jj*h)/h;

            /* Diagonal */
            if ( i >= M )
                return PLASMA_SUCCESS;

            if ( i >=0 )
                A[j*LDA + i] = 2. * term + (0.5 - (jj+1) * h) / h;

            /* Below the diagonal */
            if (i+1 < M) {
                if (jj+1 == half)
                    A[j*LDA + i+1] = - term + (0.5 - (jj+2)*h)/h;
                else
                    A[j*LDA + i+1] = - term;
            }
        }

        /* Second half */
        for (; j<N; j++,jj++) {
            i = jj - m0;
            if ( ( i < -1 ) || (i > M))
                continue;

            if (i > 0) {
                if (jj == half)
                    A[j*LDA + i-1] = - term - (0.5 - jj*h)/h;
                else
                    A[j*LDA + i-1] = - term;
            }

            if ((i >=0) && (i < M))
                A[j*LDA + i] = 2. * term - (0.5 - (jj+1)*h)/h;

            if (i+1 < M)
                A[j*LDA + i+1] = - term + (0.5 - (jj+2)*h)/h;
        }
    }
    break;

    /*
     * See [1] J. Demmel, Applied Numerical Linear Algebra, SIAM,
     *         Philadelphia, 1997
     *
     * Returns a matrix defined by:
     *    A = D * ( I + 1e-7* rand(n)), where D = diag(10^{14*(0:n-1)/n})
     *
     */
    case PlasmaMatrixDemmel:
    {
        PLASMA_Complex64_t dii;
        int ii, jj;

        /* Randomize the matrix */
        CORE_zplrnt( M, N, A, LDA, gM, m0, n0, seed );


        for (j=0,jj=n0; j<N; j++,jj++) {
            for (i=0,ii=m0; i<M; i++,ii++) {
                dii = cpow( 10. , 14. * ii / gM );
                A[j*LDA+i] *= dii * ( (jj == ii) ? 1. : 1.e-7 );
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000000
     *
     * Inverse of an upper Hessenberg matrix
     *
     * A = gallery('invhess',x,y), where x is a length n vector and y
     * is a length n-1 vector, returns the matrix whose lower triangle
     * agrees with that of ones(n,1)*x' and whose strict upper
     * triangle agrees with that of [1 y]*ones(1,n).
     *
     * The matrix is nonsingular if x(1) ~= 0 and x(i+1) ~= y(i) for
     * all i, and its inverse is an upper Hessenberg matrix. Argument
     * y defaults to -x(1:n-1).
     *
     * If x is a scalar, invhess(x) is the same as invhess(1:x).
     *
     * Here: gallery('invhess', gM)
     */
    case PlasmaMatrixInvhess:
    {
        int ii, jj;
        for (j=0,jj=n0+1; j<N; j++,jj++) {
            for (i=0,ii=m0+1; i<M; i++,ii++) {
                if ( jj <= ii ) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)( jj );
                } else {
                    A[j*LDA+i] = (PLASMA_Complex64_t)( -ii );
                }
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1019317
     *
     * Cauchy matrix
     *
     * Returns an n-by-n matrix C such that, C(i,j) = 1/(i + j).
     *
     */
    case PlasmaMatrixCauchy:
    {
        int ii, jj;
        for (j=0,jj=n0+1; j<N; j++,jj++) {
            for (i=0,ii=m0+1; i<M; i++,ii++) {
                A[j*LDA+i] = (PLASMA_Complex64_t)( 1. ) / (PLASMA_Complex64_t)( ii+jj );
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/hilb.html
     *
     * Hilbert Matrix
     *
     * The Hilbert matrix is a notable example of a poorly conditioned
     * matrix. The elements of the Hilbert matrices are:
     *   H(i,j) = 1/(i * + j – 1).
     *
     */
    case PlasmaMatrixHilb:
    {
        PLASMA_Complex64_t tmp = (PLASMA_Complex64_t)( m0 + n0 + 1. );
        for (j=0; j<N; j++) {
            for (i=0; i<M; i++) {
                A[j*LDA+i] = (PLASMA_Complex64_t)1. / (PLASMA_Complex64_t)( tmp + i + j );
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000062
     *
     * Lotkin matrix
     *
     * Returns the Hilbert matrix with its first row altered to all
     * ones. The Lotkin matrix A is nonsymmetric, ill-conditioned, and
     * has many negative eigenvalues of small magnitude. Its inverse
     * has integer entries and is known explicitly.
     *
     */
    case PlasmaMatrixLotkin:
    {
        PLASMA_Complex64_t tmp = (PLASMA_Complex64_t)( m0 + n0 + 1. );
        if (m0 == 0) {
            for (j=0; j<N; j++) {
                A[j*LDA] = (PLASMA_Complex64_t)1.;
                for (i=1; i<M; i++) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)1. / (PLASMA_Complex64_t)( tmp + i + j );
                }
            }
        } else {
            for (j=0; j<N; j++) {
                for (i=0; i<M; i++) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)1. / (PLASMA_Complex64_t)( tmp + i + j );
                }
            }
        }
    }
    break;

    /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000083
     *
     * Orthogonal and nearly orthogonal matrices
     *
     * Returns the matrix Q of order n, such that:
     *    Q(i,j) = sqrt(2/(n+1)) * sin(i*j*pi/(n+1))
     *
     * Symmetric eigenvector matrix for second difference matrix.
     *
     */
    case PlasmaMatrixOrthog: /* Default: k=1 */
    {
        PLASMA_Complex64_t sqrtn = (PLASMA_Complex64_t) sqrt( 2. / (gN+1.) );
        double scale = pi / (double)(gN+1.);

        int ii, jj;
        for (j=0,jj=n0+1; j<N; j++,jj++) {
            for (i=0,ii=m0+1; i<M; i++,ii++) {
                A[j*LDA+i] = sqrtn * (PLASMA_Complex64_t) sin( (double)ii * (double)jj * scale );
            }
        }
    }
    break;

     /*
     * See http://www.mathworks.fr/fr/help/matlab/ref/wilkinson.html
     *
     * Wilkinson's eigenvalue test matrix
     *
     * Returns one of J. H. Wilkinson's eigenvalue test matrices. It
     * is a symmetric, tridiagonal matrix with pairs of nearly, but
     * not exactly, equal eigenvalues.
     *
     */
   case PlasmaMatrixWilkinson:
    {
        if (gM != gN) {
            coreblas_error(6, "Illegal value of gM (Matrix must be square for Wilkinson)");
            return -6;
        }

        int ii, jj;
        for (j=0,jj=n0; j<N; j++,jj++) {
            for (i=0,ii=m0; i<M; i++,ii++) {
                if (ii == jj) {
                    PLASMA_Complex64_t tmp = (PLASMA_Complex64_t)(( (gN - 1 - ii) < ii ) ? gN - 1 - ii : ii );
                     A[j*LDA+i] = (PLASMA_Complex64_t)(gN - 2. * tmp - 1.) / 2.;
                }
                else if ( (ii == jj+1) || (ii == jj-1) ) {
                     A[j*LDA+i] = (PLASMA_Complex64_t)1.;
                }
                else {
                     A[j*LDA+i] = (PLASMA_Complex64_t)0.;
                }
           }
       }
    }
    break;

    /*
     * See [1] L. V. Foster, Gaussian Elimination with Partial
     *         Pivoting Can Fail in Practice, SIAM J. Matrix
     *         Anal. Appl., 15 (1994), pp. 1354-1362.
     *
     *     [2] L. V. Foster, The growth factor and efficiency of
     *         Gaussian elimination with rook pivoting,
     *         J. Comput. Appl. Math., 86 (1997), pp. 177-194
     *
     * A pathological case for LU with gaussian elimination.
     *
     */
    case PlasmaMatrixFoster: /* Default: k=h=c=1 */
    {
        double k=1., h=1., c=1.;

        int ii, jj;
        for (j=0,jj=n0; j<N; j++,jj++) {
            for (i=0,ii=m0; i<M; i++,ii++) {

                if (ii == jj) {
                    if (jj == 0)
                        A[j*LDA+i] = (PLASMA_Complex64_t)1.;
                    else if (jj == gN-1)
                        A[j*LDA+i] = (PLASMA_Complex64_t)(1. - (1. / c) - (k*h)/2. );
                    else
                        A[j*LDA+i] = (PLASMA_Complex64_t)(1. - (k*h)/2. );
                }
                else if (jj == 0) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)(-k*h/2.);
                }
                else if (jj == gN-1) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)(-1./c);
                }
                else if (ii > jj) {
                    A[j*LDA+i] = (PLASMA_Complex64_t)(-k*h);
                }
                else {
                    A[j*LDA+i] = (PLASMA_Complex64_t)0.;
                }
            }
        }
    }
    break;

    /*
     * See [3] S. J. Wright, A collection of problems for which
     *         Gaussian elimination with partial pivoting is unstable,
     *         SIAM J. SCI. STATIST. COMPUT., 14 (1993), pp. 231-238.
     *
     * A pathological case for LU with gaussian elimination.
     *
     */
    /*
     * Default: h=0.01, M=[-10 -19, 19 30]. Then,
     *   exp(h*M)=[0.9048 0.8270, 1.2092 1.3499]
     */
    case PlasmaMatrixWright:
    {
        int ii, jj;
        for (j=0,jj=n0; j<N; j++,jj++) {
            for (i=0,ii=m0; i<M; i++,ii++) {

                if (ii == jj)
                    A[j*LDA+i] = (PLASMA_Complex64_t)1.;
                else if ((ii == jj + 2) && (jj % 2 == 0))
                    A[j*LDA+i] = (PLASMA_Complex64_t)(-0.9048);
                else if ((ii == jj + 3) && (jj % 2 == 0))
                    A[j*LDA+i] = (PLASMA_Complex64_t)(-1.2092);
                else if ((ii == jj + 2) && (jj % 2 == 1))
                    A[j*LDA+i] = (PLASMA_Complex64_t)(-0.8270);
                else if ((ii == jj + 3) && (jj % 2 == 1))
                    A[j*LDA+i] = (PLASMA_Complex64_t)(-1.3499);
                else if ((jj == gM-2) && (ii == 0))
                    A[j*LDA+i] = (PLASMA_Complex64_t)1.;
                else if ((jj == gM-1) && (ii == 1))
                    A[j*LDA+i] = (PLASMA_Complex64_t)1.;
                else
                    A[j*LDA+i] = (PLASMA_Complex64_t)0.;
            }
        }
    }
    break;

    /*
     * Generates a pathological case for LU with gaussian elimination.
     *
     * Returns a random matrix on which, the columns from N/4 to N/2
     * are scaled down by eps.
     * These matrices fails on LU with partial pivoting, but Hybrid
     * LU-QR algorithms manage to recover the scaled down columns.
     *
     */
    case PlasmaMatrixLangou:
    {
        PLASMA_Complex64_t eps = (PLASMA_Complex64_t) LAPACKE_dlamch_work( 'e' );
        int ii, jj, mm, minMN;
        int firstcol, lastcol;

         /* Create random tile */
         CORE_zplrnt( M, N, A, LDA, gM, m0, n0, seed );

         /* Scale down the columns gN/4 to gN/2 below the diagonal */
         minMN = min( gM, gN );
         firstcol = minMN / 4;
         lastcol  = minMN / 2;

         if ( (m0 >= n0) && ((n0+N) >= firstcol ) && (n0 < lastcol) ) {

             jj = max( n0, firstcol );
             j = jj - n0;

             for (; j<N && jj < lastcol; j++,jj++) {

                 ii = max( m0, jj );
                 i  = ii - m0;
                 mm = M - i;
                 cblas_zscal( mm, CBLAS_SADDR(eps), A + j*LDA + i, 1);
             }
         }
    }
    break;

    default:
        coreblas_error(1, "Illegal value of mtxtype");
        return -1;
    }

    return PLASMA_SUCCESS;
}
