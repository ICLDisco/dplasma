/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasma/types_lapack.h"
#include "dplasmaaux.h"
#include "cores/core_blas.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

#include "ztrmm_LLN.h"
#include "ztrmm_LLT.h"
#include "ztrmm_LUN.h"
#include "ztrmm_LUT.h"
#include "ztrmm_RLN.h"
#include "ztrmm_RLT.h"
#include "ztrmm_RUN.h"
#include "ztrmm_RUT.h"

#define MAX_SHAPES 2

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_ztrmm_New - Generates parsec taskpool to compute:
 *
 *  B = alpha*op( A )*B or B = alpha*B*op( A ).
 *
 *  WARNING: The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = dplasmaLeft:  A*X = B
 *          = dplasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = dplasmaUpper: Upper triangle of A is stored;
 *          = dplasmaLower: Lower triangle of A is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = dplasmaNoTrans:   A is transposed;
 *          = dplasmaTrans:     A is not transposed;
 *          = dplasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = dplasmaNonUnit: A is non unit;
 *          = dplasmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          Descriptor of the triangular matrix A of size N-by-N.
 *          The triangular matrix A. If uplo = dplasmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = dplasmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = dplasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          Descriptor of the N-by-NRHS right hand side B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_ztrmm_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_ztrmm
 * @sa dplasma_ztrmm_Destruct
 * @sa dplasma_ctrmm_New
 * @sa dplasma_dtrmm_New
 * @sa dplasma_strmm_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_ztrmm_New( dplasma_enum_t side,  dplasma_enum_t uplo,
                   dplasma_enum_t trans, dplasma_enum_t diag,
                   dplasma_complex64_t alpha,
                   const parsec_tiled_matrix_t *A,
                   parsec_tiled_matrix_t *B )
{
    parsec_taskpool_t *parsec_trmm = NULL;
    dplasma_data_collection_t * ddc_A = dplasma_wrap_data_collection((parsec_tiled_matrix_t*)A);
    dplasma_data_collection_t * ddc_B = dplasma_wrap_data_collection((parsec_tiled_matrix_t*)B);

    /* Check input arguments */
    if (side != dplasmaLeft && side != dplasmaRight) {
        dplasma_error("dplasma_ztrmm_New", "illegal value of side");
        return NULL /*-1*/;
    }
    if (uplo != dplasmaUpper && uplo != dplasmaLower) {
        dplasma_error("dplasma_ztrmm_New", "illegal value of uplo");
        return NULL /*-2*/;
    }
    if (trans != dplasmaConjTrans && trans != dplasmaNoTrans && trans != dplasmaTrans ) {
        dplasma_error("dplasma_ztrmm_New", "illegal value of trans");
        return NULL /*-3*/;
    }
    if (diag != dplasmaUnit && diag != dplasmaNonUnit) {
        dplasma_error("dplasma_ztrmm_New", "illegal value of diag");
        return NULL /*-4*/;
    }

    if ( side == dplasmaLeft ) {
        if ( uplo == dplasmaLower ) {
            if ( trans == dplasmaNoTrans ) {
                parsec_trmm = (parsec_taskpool_t*)parsec_ztrmm_LLN_new(
                    side, uplo, trans, diag, alpha,
                    ddc_A, ddc_B);
            } else { /* trans =! dplasmaNoTrans */
                parsec_trmm = (parsec_taskpool_t*)parsec_ztrmm_LLT_new(
                    side, uplo, trans, diag, alpha,
                    ddc_A, ddc_B);
            }
        } else { /* uplo = dplasmaUpper */
            if ( trans == dplasmaNoTrans ) {
                parsec_trmm = (parsec_taskpool_t*)parsec_ztrmm_LUN_new(
                    side, uplo, trans, diag, alpha,
                    ddc_A, ddc_B);
            } else { /* trans =! dplasmaNoTrans */
                parsec_trmm = (parsec_taskpool_t*)parsec_ztrmm_LUT_new(
                    side, uplo, trans, diag, alpha,
                    ddc_A, ddc_B);
            }
        }
    } else { /* side == dplasmaRight */
        if ( uplo == dplasmaLower ) {
            if ( trans == dplasmaNoTrans ) {
                parsec_trmm = (parsec_taskpool_t*)parsec_ztrmm_RLN_new(
                    side, uplo, trans, diag, alpha,
                    ddc_A, ddc_B);
            } else { /* trans =! dplasmaNoTrans */
                parsec_trmm = (parsec_taskpool_t*)parsec_ztrmm_RLT_new(
                    side, uplo, trans, diag, alpha,
                    ddc_A, ddc_B);
            }
        } else { /* uplo = dplasmaUpper */
            if ( trans == dplasmaNoTrans ) {
                parsec_trmm = (parsec_taskpool_t*)parsec_ztrmm_RUN_new(
                    side, uplo, trans, diag, alpha,
                    ddc_A, ddc_B);
            } else { /* trans =! dplasmaNoTrans */
                parsec_trmm = (parsec_taskpool_t*)parsec_ztrmm_RUT_new(
                    side, uplo, trans, diag, alpha,
                    ddc_A, ddc_B);
            }
        }
    }

    /* When supporting LAPACK we can't assume both matrixes have the same layout, e.g. LDA.
     * Therefore, generate types for both.
     */
    int shape = 0;
    dplasma_setup_adtt_all_loc( ddc_A,
                                parsec_datatype_double_complex_t,
                                PARSEC_MATRIX_FULL/*uplo*/, 1/*diag:for PARSEC_MATRIX_UPPER or PARSEC_MATRIX_LOWER types*/,
                                &shape);

    dplasma_setup_adtt_all_loc( ddc_B,
                                parsec_datatype_double_complex_t,
                                PARSEC_MATRIX_FULL/*uplo*/, 1/*diag:for PARSEC_MATRIX_UPPER or PARSEC_MATRIX_LOWER types*/,
                                &shape);

    assert(shape == MAX_SHAPES);

    return parsec_trmm;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_ztrmm_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_ztrmm_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_ztrmm_New
 * @sa dplasma_ztrmm
 *
 ******************************************************************************/
void
dplasma_ztrmm_Destruct( parsec_taskpool_t *tp )
{
    parsec_ztrmm_LLN_taskpool_t *otrmm = (parsec_ztrmm_LLN_taskpool_t *)tp;
    dplasma_clean_adtt_all_loc(otrmm->_g_ddescA, MAX_SHAPES);
    dplasma_clean_adtt_all_loc(otrmm->_g_ddescB, MAX_SHAPES);
    dplasma_data_collection_t * ddc_A = otrmm->_g_ddescA;
    dplasma_data_collection_t * ddc_B = otrmm->_g_ddescB;

    parsec_taskpool_free(tp);

    /* free the dplasma_data_collection_t */
    dplasma_unwrap_data_collection(ddc_A);
    dplasma_unwrap_data_collection(ddc_B);

}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_ztrmm - Computes:
 *
 *  B = alpha*op( A )*B or B = alpha*B*op( A ).
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = dplasmaLeft:  A*X = B
 *          = dplasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = dplasmaUpper: Upper triangle of A is stored;
 *          = dplasmaLower: Lower triangle of A is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = dplasmaNoTrans:   A is transposed;
 *          = dplasmaTrans:     A is not transposed;
 *          = dplasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = dplasmaNonUnit: A is non unit;
 *          = dplasmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          Descriptor of the triangular matrix A of size N-by-N.
 *          The triangular matrix A. If uplo = dplasmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = dplasmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = dplasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          Descriptor of the N-by-NRHS right hand side B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa dplasma_ztrmm
 * @sa dplasma_ztrmm_Destruct
 * @sa dplasma_ctrmm_New
 * @sa dplasma_dtrmm_New
 * @sa dplasma_strmm_New
 *
 ******************************************************************************/
int
dplasma_ztrmm( parsec_context_t *parsec,
               dplasma_enum_t side,  dplasma_enum_t uplo,
               dplasma_enum_t trans, dplasma_enum_t diag,
               dplasma_complex64_t alpha,
               const parsec_tiled_matrix_t *A,
               parsec_tiled_matrix_t *B)
{
    parsec_taskpool_t *parsec_ztrmm = NULL;

    /* Check input arguments */
    if (side != dplasmaLeft && side != dplasmaRight) {
        dplasma_error("dplasma_ztrmm", "illegal value of side");
        return -1;
    }
    if (uplo != dplasmaUpper && uplo != dplasmaLower) {
        dplasma_error("dplasma_ztrmm", "illegal value of uplo");
        return -2;
    }
    if (trans != dplasmaConjTrans && trans != dplasmaNoTrans && trans != dplasmaTrans ) {
        dplasma_error("dplasma_ztrmm", "illegal value of trans");
        return -3;
    }
    if (diag != dplasmaUnit && diag != dplasmaNonUnit) {
        dplasma_error("dplasma_ztrmm", "illegal value of diag");
        return -4;
    }

    if ( (A->m != A->n) ||
         (( side == dplasmaLeft )  && (A->n != B->m)) ||
         (( side == dplasmaRight ) && (A->n != B->n)) ) {
        dplasma_error("dplasma_ztrmm_New", "illegal matrix A");
        return -6;
    }

    parsec_ztrmm = dplasma_ztrmm_New(side, uplo, trans, diag, alpha, A, B);

    if ( parsec_ztrmm != NULL )
    {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_ztrmm);
        dplasma_wait_until_completion(parsec);
        dplasma_ztrmm_Destruct( parsec_ztrmm );
        return 0;
    }
    else {
        return -101;
    }
}
