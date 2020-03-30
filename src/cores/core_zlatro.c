/**
 *
 * @file core_zlatro.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup dplasma_core_complex64
 *
 *  CORE_zlatro transposes a m-by-n matrix out of place.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          = PlasmaUpper: the upper triangle of A and the lower triangle of B
 *          are referenced.
 *          = PlasmaLower: the lower triangle of A and the upper triangle of B
 *          are referenced.
 *          = PlasmaUpperLower: All A and B are referenced.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          conjugate transposed:
 *          = PlasmaNoTrans:   B is a copy of A (equivalent to zlacpy);
 *          = PlasmaTrans:     B is the transpose of A;
 *          = PlasmaConjTrans: B is the conjugate transpose of A.
 *
 * @param[in] M
 *         Number of rows of the matrix A and number of columns of the matrix B, if trans == Pasma[Conj]Trans.
 *         Number of rows of the matrix A and the matrix B, if trans == PasmaNoTrans.
 *
 * @param[in] N
 *         Number of columns of the matrix A and number of rows of the matrix B, if trans == Pasma[Conj]Trans.
 *         Number of columns of the matrix A and of the matrix B, if trans == PlasmaNoTrans.
 *
 * @param[in] A
 *         Matrix of size LDA-by-N, if trans == Pasma[Conj]Trans.
 *         Matrix of size LDA-by-M, if trans == PasmaNoTrans.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.
 *         LDA >= max(1,M), if trans == Pasma[Conj]Trans.
 *         LDA >= max(1,N), if trans == PasmaNoTrans.
 *
 * @param[out] B
 *         Matrix of size LDB-by-M, if trans == Pasma[Conj]Trans.
 *         Matrix of size LDB-by-N, if trans == PasmaNoTrans.
 *
 * @param[in] LDB
 *         The leading dimension of the array B.
 *         LDB >= max(1,N), if trans == Pasma[Conj]Trans.
 *         LDB >= max(1,M), if trans == PasmaNoTrans.
 *
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlatro = PCORE_zlatro
#define CORE_zlatro PCORE_zlatro
#endif
int CORE_zlatro(PLASMA_enum uplo, PLASMA_enum trans,
                int M, int N,
                const PLASMA_Complex64_t *A, int LDA,
                      PLASMA_Complex64_t *B, int LDB)
{
    int i, j;

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower) && (uplo != PlasmaUpperLower) ) {
        coreblas_error(1, "Illegal value of uplo");
        return -1;
    }
    if ((trans != PlasmaConjTrans) && (trans != PlasmaNoTrans) && (trans != PlasmaTrans) ) {
        coreblas_error(2, "Illegal value of trans");
        return -2;
    }
    if (M < 0) {
        coreblas_error(3, "Illegal value of M");
        return -3;
    }
    if (N < 0) {
        coreblas_error(4, "Illegal value of N");
        return -4;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        coreblas_error(6, "Illegal value of LDA");
        return -6;
    }
    if ( (LDB < max(1,N)) && (N > 0) ) {
        coreblas_error(8, "Illegal value of LDB");
        return -8;
    }

    if (trans == PlasmaNoTrans) {
        CORE_zlacpy(uplo, M, N, A, LDA, B, LDB);
    }
    else {
        if (trans == PlasmaConjTrans) {
            if(uplo == PlasmaUpper) {
                for(j=0; j<N; j++)
                    for(i=0; i<min(j+1,M); i++)
                        B[j+i*LDB] = conj(A[i+j*LDA]);
            }
            else if(uplo == PlasmaLower) {
                for(j=0;j<N;j++)
                    for(i=j;i<M;i++)
                        B[j+i*LDB] = conj(A[i+j*LDA]);
            }
            else {
                for(j=0;j<N;j++)
                    for(i=0;i<M;i++)
                        B[j+i*LDB] = conj(A[i+j*LDA]);
            }
        }
        else {
            if(uplo==PlasmaUpper) {
                for(j=0;j<N;j++)
                    for(i=0;i<min(j+1,M);i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
            else if(uplo==PlasmaLower) {
                for(j=0;j<N;j++)
                    for(i=j;i<M;i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
            else {
                for(j=0;j<N;j++)
                    for(i=0;i<M;i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
        }
    }

    return PLASMA_SUCCESS;
}
