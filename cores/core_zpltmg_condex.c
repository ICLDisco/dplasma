/**
 *
 * @file core_zpltmg_condex.c
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
#include <string.h>
#include <lapacke.h>
#include "parsec/parsec_config.h"
#include "dplasma.h"
#include "cores/dplasma_cores.h"
#include "dplasma_zcores.h"

/***************************************************************************//**
 *
 * @ingroup CORE_parsec_complex64_t
 *
 *  CORE_zpltmg_condexq generates the Q used in condex matrix generation
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999898
 *  gallery('condex',n,4,100)
 *
 *  Returns a "counter-example" matrix to a condition estimator. It has order n
 *  and scalar parameter theta (default 100).
 *
 *  LAPACK (RCOND): It is the inverse of this matrix that is a counter-example.
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the matrix Q used in condex generation. M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A to be generated. N >= 0.
 *
 * @param[out] Q
 *         On entry, the M-by-3 matrix to be initialized.
 *         On exit, the housholder reflectors required for condex generation.
 *
 * @param[in] LDQ
 *         The leading dimension of the matrix Q. LDQ >= max(1,M).
 *
 ******************************************************************************/
void CORE_zpltmg_condexq( int M, int N, parsec_complex64_t *Q, int LDQ )
{
    parsec_complex64_t tau[3];
    parsec_complex64_t *tQ = Q;
    int i;

    /* First column is [ 1 ... 1 ] */
    for( i=0; i < M; i++, tQ++ )
        *tQ = (parsec_complex64_t)1.0;

    /* Second column is [1 0 0 ... 0] */
    tQ = Q + LDQ;
    *tQ = (parsec_complex64_t)1.;
    tQ++;
    memset( tQ, 0, (M-1) * sizeof(parsec_complex64_t) );

    /* Third column is  (-1)^i * (1. + i / (N-1)) */
    tQ = Q + 2 * LDQ;
    for( i=0; i<M; i++, tQ++ )
        *tQ = (parsec_complex64_t)( cpow( -1.0, (double)i ) * (1.0 + (double)i/(N-1) ) );

    /* Generate orthogonal projector */
    LAPACKE_zgeqrf( LAPACK_COL_MAJOR, M, 3,    Q, LDQ, tau );
    LAPACKE_zungqr( LAPACK_COL_MAJOR, M, 3, 3, Q, LDQ, tau );

    return;
}
