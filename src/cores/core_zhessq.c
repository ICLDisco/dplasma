/**
 *
 * @file core_zhessq.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 **/
#include <math.h>
#include <lapacke.h>
#include "common.h"
#include "sumsq_update.h"

/*****************************************************************************
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_zhessq returns the values scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *                     i,j
 *
 * with i from 0 to N-1 and j form 0 to N-1. The value of sumsq is
 * assumed to be at least unity and the value of ssq will then satisfy
 *
 *    1.0 .le. ssq .le. ( sumsq + 2*n*n ).
 *
 * scale is assumed to be non-negative and scl returns the value
 *
 *    scl = max( scale, abs( real( A( i, j ) ) ), abs( aimag( A( i, j ) ) ) ),
 *          i,j
 *
 * scale and sumsq must be supplied in SCALE and SUMSQ respectively.
 * SCALE and SUMSQ are overwritten by scl and ssq respectively.
 *
 * The routine makes only one pass through the tile triangular part of the
 * hermitian tile A defined by uplo.
 * See also LAPACK zlassq.f
 *
 *******************************************************************************
 *
 *  @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the hermitian matrix A is to be referenced as follows:
 *          = PlasmaLower:     Only the lower triangular part of the
 *                             hermitian matrix A is to be referenced.
 *          = PlasmaUpper:     Only the upper triangular part of the
 *                             hermitian matrix A is to be referenced.
 *
 *  @param[in] N
 *          The number of columns and rows in the tile A.
 *
 *  @param[in] A
 *          The N-by-N matrix on which to compute the norm.
 *
 *  @param[in] LDA
 *          The leading dimension of the tile A. LDA >= max(1,N).
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
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zhessq = PCORE_zhessq
#define CORE_zhessq PCORE_zhessq
#endif
int CORE_zhessq(PLASMA_enum uplo, int N,
                const PLASMA_Complex64_t *A, int LDA,
                double *scale, double *sumsq)
{
    int i, j;
    double *ptr;

    if ( uplo == PlasmaUpper ) {
        for(j=0; j<N; j++) {
            ptr = (double*) ( A + j * LDA );

            for(i=0; i<j; i++, ptr++) {
		sumsq_update( 2, scale, sumsq, ptr );
		
#if defined(PRECISION_z) || defined(PRECISION_c)
		ptr++;
		sumsq_update( 2, scale, sumsq, ptr );
#endif
            }

            /* Diagonal */
	    sumsq_update( 1, scale, sumsq, ptr );
		
#if defined(PRECISION_z) || defined(PRECISION_c)
	    ptr++;
#endif
        }
    } else {

        for(j=0; j<N; j++) {
            ptr = (double*) ( A + j * LDA + j);

            /* Diagonal */
	    sumsq_update( 1, scale, sumsq, ptr );
	    ptr++;
		
#if defined(PRECISION_z) || defined(PRECISION_c)
	    ptr++;
#endif

            for(i=j+1; i<N; i++, ptr++) {
		sumsq_update( 2, scale, sumsq, ptr );
		
#if defined(PRECISION_z) || defined(PRECISION_c)
		ptr++;
		sumsq_update( 2, scale, sumsq, ptr );
#endif
            }
        }
    }
    return PLASMA_SUCCESS;
}
