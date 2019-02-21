/**
 *
 * @file core_dsyrfb.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:30 2016
 *
 **/
#include <lapacke.h>
#include "parsec/parsec_config.h"
#include "dplasma.h"
#include "cores/dplasma_cores.h"
#include "dplasma_dcores.h"

#undef COMPLEX
#define REAL

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dsyrfb overwrites the symmetric complex N-by-N tile C with
 *
 *    Q**T*C*Q
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_dgeqrt. Only PlasmaLower supported!
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *         @arg PlasmaLower : the upper part of the symmetric matrix C
 *                            is not referenced.
 *         @arg PlasmaUpper : the lower part of the symmetric matrix C
 *                            is not referenced (not supported).
 *
 * @param[in] n
 *         The number of rows/columns of the tile C.  N >= 0.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q. K >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] nb
 *         The blocking size.  NB >= 0.
 *
 * @param[in] A
 *         The i-th column must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_dgeqrt in the first k columns of its array argument A.
 *
 * @param[in] lda
 *         The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the symmetric N-by-N tile C.
 *         On exit, C is overwritten by Q**T*C*Q.
 *
 * @param[in] ldc
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] ldwork
 *         The dimension of the array WORK. LDWORK >= max(1,N);
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsyrfb = PCORE_dsyrfb
#define CORE_dsyrfb PCORE_dsyrfb
#define CORE_dormlq PCORE_dormlq
#define CORE_dormqr PCORE_dormqr
int  CORE_dormlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const double *V, int LDV,
                 const double *T, int LDT,
                 double *C, int LDC,
                 double *WORK, int LDWORK);
int  CORE_dormqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const double *V, int LDV,
                 const double *T, int LDT,
                 double *C, int LDC,
                 double *WORK, int LDWORK);
#endif
int CORE_dsyrfb( PLASMA_enum uplo, int n,
                 int k, int ib, int nb,
                 const double *A, int lda,
                 const double *T, int ldt,
                 double *C, int ldc,
                 double *WORK, int ldwork )
{
    double tmp;
    int i, j;

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        coreblas_error(1, "Illegal value of uplo");
        return -1;
    }
    if (n < 0) {
        coreblas_error(2, "Illegal value of n");
        return -2;
    }
    if (k < 0) {
        coreblas_error(3, "Illegal value of k");
        return -3;
    }
    if (ib < 0) {
        coreblas_error(4, "Illegal value of ib");
        return -4;
    }
    if (nb < 0) {
        coreblas_error(5, "Illegal value of nb");
        return -5;
    }
    if ( (lda < max(1,n)) && (n > 0) ) {
        coreblas_error(7, "Illegal value of lda");
        return -7;
    }
    if ( (ldt < max(1,ib)) && (ib > 0) ) {
        coreblas_error(9, "Illegal value of ldt");
        return -9;
    }
    if ( (ldc < max(1,n)) && (n > 0) ) {
        coreblas_error(11, "Illegal value of ldc");
        return -11;
    }

    if (uplo == PlasmaLower) {
        /* Rebuild the symmetric block: WORK <- C */
        for (j = 0; j < n; j++) {
            *(WORK + j + j * ldwork) =  *(C + ldc*j + j);
            for (i = j+1; i < n; i++){
                tmp = *(C + i + j*ldc);
                *(WORK + i + j * ldwork) = tmp;
                *(WORK + j + i * ldwork) = ( tmp );
            }
        }

        /* Left */
        CORE_dormqr(PlasmaLeft, PlasmaTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Right */
        CORE_dormqr(PlasmaRight, PlasmaNoTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);

        /*
         * Copy back the final result to the lower part of C
         */
        LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, lapack_const(PlasmaLower), n, n, WORK, ldwork, C, ldc );
    }
    else {
        /* Rebuild the symmetric block: WORK <- C */
        for (j = 0; j < n; j++) {
            for (i = 0; i < j; i++){
                tmp = *(C + i + j*ldc);
                *(WORK + i + j * ldwork) = tmp;
                *(WORK + j + i * ldwork) = ( tmp );
            }
            *(WORK + j + j * ldwork) =  *(C + ldc*j + j);
        }

        /* Right */
        CORE_dormlq(PlasmaRight, PlasmaTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Left */
        CORE_dormlq(PlasmaLeft, PlasmaNoTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);

        /*
         * Copy back the final result to the upper part of C
         */
        LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, lapack_const(PlasmaUpper), n, n, WORK, ldwork, C, ldc );
    }
    return 0;
}
