/**
 *
 * @file core_zgessq.c
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
#include "cores/dplasma_cores.h"
#include "dplasma_zcores.h"

#define COMPLEX

#define UPDATE( __nb, __value )                                         \
    if (__value != 0. ){                                                \
        if ( *scale < __value ) {                                       \
            *sumsq = __nb + (*sumsq) * ( *scale / __value ) * ( *scale / __value ); \
            *scale = __value;                                           \
        } else {                                                        \
            *sumsq = *sumsq + __nb * ( __value / *scale ) *  ( __value / *scale ); \
        }                                                               \
    }

/*****************************************************************************
 *
 * @ingroup CORE_parsec_complex64_t
 *
 *  CORE_zgessq returns the values scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *                     i,j
 *
 * with i from 0 to M-1 and j form 0 to N-1. The value of sumsq is
 * assumed to be at least unity and the value of ssq will then satisfy
 *
 *    1.0 .le. ssq .le. ( sumsq + 2*m*n ).
 *
 * scale is assumed to be non-negative and scl returns the value
 *
 *    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
 *           i
 *
 * scale and sumsq must be supplied in SCALE and SUMSQ respectively.
 * SCALE and SUMSQ are overwritten by scl and ssq respectively.
 *
 * The routine makes only one pass through the tile A.
 * See also LAPACK zlassq.f
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of rows in the tile A.
 *
 *  @param[in] N
 *          The number of columns in the tile A.
 *
 *  @param[in] A
 *          The M-by-N matrix on which to compute the norm.
 *
 *  @param[in] LDA
 *          The leading dimension of the tile A. LDA >= max(1,M).
 *
 *  @param[in,out] scale
 *          On entry, the value  scale  in the equation above.
 *          On exit, scale is overwritten with the value scl.
 *
 *  @param[in,out] sumsq
 *          On entry, the value  sumsq  in the equation above.
 *          On exit, SUMSQ is overwritten with the value ssq.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval -k, the k-th argument had an illegal value
 *
 */
int CORE_zgessq(int M, int N,
                const parsec_complex64_t *A, int LDA,
                double *scale, double *sumsq)
{
    int i, j;
    double tmp;
    double *ptr;

    for(j=0; j<N; j++) {
        ptr = (double*) ( A + j * LDA );
        for(i=0; i<M; i++, ptr++) {
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );

#ifdef COMPLEX
            ptr++;
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );
#endif
        }
    }
    return PLASMA_SUCCESS;
}
