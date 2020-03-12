/**
 *
 * @file core_zlantr.c
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
#include "core_blas.h"

#define LAPACKE_CORRECT_DLANTR

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zlantr returns the value
 *
 *     zlantr = ( max(abs(A(i,j))), NORM = PlasmaMaxNorm
 *              (
 *              ( norm1(A),         NORM = PlasmaOneNorm
 *              (
 *              ( normI(A),         NORM = PlasmaInfNorm
 *              (
 *              ( normF(A),         NORM = PlasmaFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(A(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = PlasmaMaxNorm: Max norm
 *          = PlasmaOneNorm: One norm
 *          = PlasmaInfNorm: Infinity norm
 *          = PlasmaFrobeniusNorm: Frobenius norm
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *          If uplo == PlasmaUpper, M <= N. When M = 0, CORE_zlantr returns 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *          If uplo == PlasmaLower, N <= M. When N = 0, CORE_zlantr returns 0.
 *
 * @param[in] A
 *          The LDA-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in,out] work
 *          Array of dimension (MAX(1,LWORK)), where LWORK >= M when norm =
 *          PlasmaInfNorm, or LWORK >= N when norm = PlasmaOneNorm; otherwise,
 *          work is not referenced.
 *
 * @param[out] normA
 *          On exit, normA is the norm of matrix A.
 *
 ******************************************************************************/
void CORE_zlantr(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *work, double *normA)
{
#if defined(LAPACKE_CORRECT_DLANTR)
    *normA = LAPACKE_zlantr_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        lapack_const(uplo),
        lapack_const(diag),
        M, N, A, LDA, work);
#else
    const PLASMA_Complex64_t *tmpA;
    double value;
    int i, j, imax;
    int idiag = (diag == PlasmaUnit) ? 1 : 0;

    if ( coreblas_imin(M, N) == 0 ) {
        *normA = 0;
        return;
    }

    switch ( norm ) {
    case PlasmaMaxNorm:
        if ( diag == PlasmaUnit ) {
            *normA = 1.;
        } else {
            *normA = 0.;
        }

        if ( uplo == PlasmaUpper ) {
            M = coreblas_imin(M, N);
            for (j = 0; j < N; j++) {
                tmpA = A+(j*LDA);
                imax = coreblas_imin(j+1-idiag, M);

                for (i = 0; i < imax; i++) {
                    value = cabs( *tmpA );
                    *normA = ( value > *normA ) ? value : *normA;
                    tmpA++;
                }
            }
        } else {
            N = coreblas_imin(M, N);
            for (j = 0; j < N; j++) {
                tmpA = A + j * (LDA+1) + idiag;

                for (i = j+idiag; i < M; i++) {
                    value = cabs( *tmpA );
                    *normA = ( value > *normA ) ? value : *normA;
                    tmpA++;
                }
            }
        }
        break;

    case PlasmaOneNorm:
        CORE_ztrasm( PlasmaColumnwise, uplo, diag, M, N,
                     A, LDA, work );
        if ( uplo == PlasmaLower )
            N = coreblas_imin(M,N);

        *normA = 0;
        for (i = 0; i < N; i++) {
            *normA = ( work[i] > *normA ) ? work[i] : *normA;
        }
        break;

    case PlasmaInfNorm:
        CORE_ztrasm( PlasmaRowwise, uplo, diag, M, N,
                     A, LDA, work );
        if ( uplo == PlasmaUpper )
            M = coreblas_imin(M,N);

        *normA = 0;
        for (i = 0; i < M; i++) {
            *normA = ( work[i] > *normA ) ? work[i] : *normA;
        }
        break;

    case PlasmaFrobeniusNorm:
    {
        double scale = 0.;
        double sumsq = 1.;
        CORE_ztrssq( uplo, diag, M, N,
                     A, LDA, &scale, &sumsq );

        *normA = scale * sqrt( sumsq );
    }
    break;
    default:
        coreblas_error(1, "Illegal value of norm");
        return;
    }
#endif
}
