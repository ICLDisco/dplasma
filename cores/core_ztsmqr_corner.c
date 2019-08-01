/**
 *
 * @file core_ztsmqr_corner.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "parsec/parsec_config.h"
#include "dplasma.h"
#include "cores/dplasma_cores.h"
#include "dplasma_zcores.h"

#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup dplasma_cores_complex64
 *
 *  CORE_ztsmqr_corner: see CORE_ztsmqr
 *
 * This kernel applies left and right transformations as depicted below:
 * |I -VT'V'| * | A1 A2'| * |I - VTV'|
 *              | A2 A3 |
 * where A1 and A3 are symmetric matrices.
 * Only the lower part is referenced.
 * This is an adhoc implementation, can be further optimized...
 *
 *******************************************************************************
 *
 * @param[in] m1
 *         The number of rows of the tile A1. m1 >= 0.
 *
 * @param[in] n1
 *         The number of columns of the tile A1. n1 >= 0.
 *
 * @param[in] m2
 *         The number of rows of the tile A2. m2 >= 0.
 *
 * @param[in] n2
 *         The number of columns of the tile A2. n2 >= 0.
 *
 * @param[in] m3
 *         The number of rows of the tile A3. m3 >= 0.
 *
 * @param[in] n3
 *         The number of columns of the tile A3. n3 >= 0.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] ib
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] nb
 *         The blocking size.  NB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. lda1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. lda2 >= max(1,M2).
 *
 * @param[in,out] A3
 *         On entry, the m3-by-n3 tile A3.
 *
 * @param[in] lda3
 *         The leading dimension of the tile A3. lda3 >= max(1,m3).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_ZTSQRT in the first k columns of its array argument V.
 *
 * @param[in] ldv
 *         The leading dimension of the array V. ldv >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. ldt >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size
 *             LDWORK-by-N1 if side == PlasmaLeft
 *             LDWORK-by-IB if side == PlasmaRight
 *
 * @param[in] ldwork
 *         The leading dimension of the array WORK.
 *             LDWORK >= max(1,IB) if side == PlasmaLeft
 *             LDWORK >= max(1,M1) if side == PlasmaRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int CORE_ztsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        parsec_complex64_t *A1, int lda1,
                        parsec_complex64_t *A2, int lda2,
                        parsec_complex64_t *A3, int lda3,
                        const parsec_complex64_t *V, int ldv,
                        const parsec_complex64_t *T, int ldt,
                        parsec_complex64_t *WORK, int ldwork)
{
    int i, j;
    PLASMA_enum side, trans;

    if ( m1 != n1 ) {
        coreblas_error(1, "Illegal value of M1, N1");
        return -1;
    }

    /*  Rebuild the symmetric block: WORK <- A1 */
    for (j = 0; j < n1; j++)
        for (i = j; i < m1; i++){
            *(WORK + i + j*ldwork) = *(A1 + i + j*lda1);
            if (i > j){
                *(WORK + j + i*ldwork) =  conj( *(WORK + i + j*ldwork) );
            }
        }

    /*  Copy the transpose of A2: WORK+nb*ldwork <- A2' */
    for (j = 0; j < n2; j++)
        for (i = 0; i < m2; i++){
            *(WORK + j + (i + nb) * ldwork) = conj( *(A2 + i + j*lda2) );
        }

    side  = PlasmaLeft;
    trans = PlasmaConjTrans;

    /*  Left application on |A1| */
    /*                      |A2| */
    CORE_ztsmqr(side, trans, m1, n1, m2, n2, k, ib,
                WORK, ldwork, A2, lda2,
                V, ldv, T, ldt,
                WORK + 3*nb*ldwork, ldwork);

    /*  Rebuild the symmetric block: WORK+2*nb*ldwork <- A3 */
    for (j = 0; j < n3; j++)
        for (i = j; i < m3; i++){
            *(WORK + i + (j + 2*nb) * ldwork) = *(A3 + i + j*lda3);
            if (i != j){
                *(WORK + j + (i + 2*nb) * ldwork) =  conj( *(WORK + i + (j + 2*nb) * ldwork) );
            }
        }
    /*  Left application on | A2'| */
    /*                      | A3 | */
    CORE_ztsmqr(side, trans, n2, m2, m3, n3, k, ib,
                WORK+nb*ldwork, ldwork, WORK+2*nb*ldwork, ldwork,
                V, ldv, T, ldt,
                WORK + 3*nb*ldwork, ldwork);

    side  = PlasmaRight;
    trans = PlasmaNoTrans;

    /*  Right application on | A1 A2' | */
    CORE_ztsmqr(side, trans, m1, n1, n2, m2, k, ib,
                WORK, ldwork, WORK+nb*ldwork, ldwork,
                V, ldv, T, ldt,
                WORK + 3*nb*ldwork, ldwork);

    /*  Copy back the final result to the lower part of A1 */
    /*  A1 = WORK */
    for (j = 0; j < n1; j++)
        for (i = j; i < m1; i++)
            *(A1 + i + j*lda1) = *(WORK + i + j*ldwork);

    /*  Right application on | A2 A3 | */
    CORE_ztsmqr(side, trans, m2, n2, m3, n3, k, ib,
                A2, lda2, WORK+2*nb*ldwork, ldwork,
                V,  ldv,  T, ldt,
                WORK + 3*nb*ldwork, ldwork);

    /*  Copy back the final result to the lower part of A3 */
    /*  A3 = WORK+2*nb*ldwork */
    for (j = 0; j < n3; j++)
        for (i = j; i < m3; i++)
            *(A3 + i + j*lda3) = *(WORK + i + (j+ 2*nb) * ldwork);

    return PLASMA_SUCCESS;
}
