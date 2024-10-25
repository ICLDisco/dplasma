/*
 * Copyright (c) 2010-2022 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 *
 * @precisions normal z -> s d c
 *
 */

#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasmaaux.h"

#include "ztrsm_LLN.h"
#include "ztrsm_LLT.h"
#include "ztrsm_LUN.h"
#include "ztrsm_LUT.h"
#include "ztrsm_RLN.h"
#include "ztrsm_RLT.h"
#include "ztrsm_RUN.h"
#include "ztrsm_RUT.h"

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_ztrsm_New - Generates parsec taskpool to compute triangular solve
 *     op( A ) * X = B or X * op( A ) = B
 *  WARNING: The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = dplasmaLeft:  op( A ) * X = B
 *          = dplasmaRight: X * op( A ) = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          = dplasmaUpper: Upper triangle of A is stored;
 *          = dplasmaLower: Lower triangle of A is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          conjugate transposed:
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
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          Descriptor of the triangular matrix A of size N-by-N.
 *          If uplo = dplasmaUpper, the leading N-by-N upper triangular part of
 *          the array A contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced. If uplo = dplasmaLower,
 *          the leading N-by-N lower triangular part of the array A contains the
 *          lower triangular matrix, and the strictly upper triangular part of A
 *          is not referenced. If diag = dplasmaUnit, the diagonal elements of A
 *          are also not referenced and are assumed to be 1.
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
 *          destroy with dplasma_ztrsm_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_ztrsm
 * @sa dplasma_ztrsm_Destruct
 * @sa dplasma_ctrsm_New
 * @sa dplasma_dtrsm_New
 * @sa dplasma_strsm_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_ztrsm_New( dplasma_enum_t side,  dplasma_enum_t uplo,
                   dplasma_enum_t trans, dplasma_enum_t diag,
                   dplasma_complex64_t alpha,
                   const parsec_tiled_matrix_t *A,
                   parsec_tiled_matrix_t *B )
{
    parsec_taskpool_t *parsec_trsm = NULL;

    if ( side == dplasmaLeft ) {
        if ( uplo == dplasmaLower ) {
            if ( trans == dplasmaNoTrans ) {
                parsec_trsm = (parsec_taskpool_t*)parsec_ztrsm_LLN_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            } else { /* trans =! dplasmaNoTrans */
                parsec_trsm = (parsec_taskpool_t*)parsec_ztrsm_LLT_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            }
        } else { /* uplo = dplasmaUpper */
            if ( trans == dplasmaNoTrans ) {
                parsec_trsm = (parsec_taskpool_t*)parsec_ztrsm_LUN_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            } else { /* trans =! dplasmaNoTrans */
                parsec_trsm = (parsec_taskpool_t*)parsec_ztrsm_LUT_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            }
        }
    } else { /* side == dplasmaRight */
        if ( uplo == dplasmaLower ) {
            if ( trans == dplasmaNoTrans ) {
                parsec_trsm = (parsec_taskpool_t*)parsec_ztrsm_RLN_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            } else { /* trans =! dplasmaNoTrans */
                parsec_trsm = (parsec_taskpool_t*)parsec_ztrsm_RLT_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            }
        } else { /* uplo = dplasmaUpper */
            if ( trans == dplasmaNoTrans ) {
                parsec_trsm = (parsec_taskpool_t*)parsec_ztrsm_RUN_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            } else { /* trans =! dplasmaNoTrans */
                parsec_trsm = (parsec_taskpool_t*)parsec_ztrsm_RUT_new(
                    side, uplo, trans, diag, alpha,
                    A,
                    B);
            }
        }
    }

    dplasma_add2arena_tile( &((parsec_ztrsm_LLN_taskpool_t*)parsec_trsm)->arenas_datatypes[PARSEC_ztrsm_LLN_DEFAULT_ADT_IDX],
                            A->mb*A->nb*sizeof(dplasma_complex64_t),
                            PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t, A->mb );

    return parsec_trsm;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_ztrsm_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_ztrsm_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_ztrsm_New
 * @sa dplasma_ztrsm
 *
 ******************************************************************************/
void
dplasma_ztrsm_Destruct( parsec_taskpool_t *tp )
{
    parsec_ztrsm_LLN_taskpool_t *otrsm = (parsec_ztrsm_LLN_taskpool_t *)tp;

    dplasma_matrix_del2arena( &otrsm->arenas_datatypes[PARSEC_ztrsm_LLN_DEFAULT_ADT_IDX] );
    parsec_taskpool_free(tp);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_ztrsm - Computes triangular solve
 *     op( A ) * X = B or X * op( A ) = B
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = dplasmaLeft:  op( A ) * X = B
 *          = dplasmaRight: X * op( A ) = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          = dplasmaUpper: Upper triangle of A is stored;
 *          = dplasmaLower: Lower triangle of A is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          conjugate transposed:
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
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          Descriptor of the triangular matrix A of size N-by-N.
 *          If uplo = dplasmaUpper, the leading N-by-N upper triangular part of
 *          the array A contains the upper triangular matrix, and the strictly
 *          lower triangular part of A is not referenced. If uplo = dplasmaLower,
 *          the leading N-by-N lower triangular part of the array A contains the
 *          lower triangular matrix, and the strictly upper triangular part of A
 *          is not referenced. If diag = dplasmaUnit, the diagonal elements of A
 *          are also not referenced and are assumed to be 1.
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
 * @sa dplasma_ztrsm_New
 * @sa dplasma_ztrsm_Destruct
 * @sa dplasma_ctrsm
 * @sa dplasma_dtrsm
 * @sa dplasma_strsm
 *
 ******************************************************************************/
int
dplasma_ztrsm( parsec_context_t *parsec,
               dplasma_enum_t side,  dplasma_enum_t uplo,
               dplasma_enum_t trans, dplasma_enum_t diag,
               dplasma_complex64_t alpha,
               const parsec_tiled_matrix_t *A,
               parsec_tiled_matrix_t *B)
{
    parsec_taskpool_t *parsec_ztrsm = NULL;

    /* Check input arguments */
    if (side != dplasmaLeft && side != dplasmaRight) {
        dplasma_error("dplasma_ztrsm", "illegal value of side");
        return -1;
    }
    if (uplo != dplasmaUpper && uplo != dplasmaLower) {
        dplasma_error("dplasma_ztrsm", "illegal value of uplo");
        return -2;
    }
    if (trans != dplasmaConjTrans && trans != dplasmaNoTrans && trans != dplasmaTrans ) {
        dplasma_error("dplasma_ztrsm", "illegal value of trans");
        return -3;
    }
    if (diag != dplasmaUnit && diag != dplasmaNonUnit) {
        dplasma_error("dplasma_ztrsm", "illegal value of diag");
        return -4;
    }

    if ( (A->m != A->n) ||
         (( side == dplasmaLeft )  && (A->n != B->m)) ||
         (( side == dplasmaRight ) && (A->n != B->n)) ) {
        dplasma_error("dplasma_ztrsm", "illegal matrix A");
        return -6;
    }

    parsec_ztrsm = dplasma_ztrsm_New(side, uplo, trans, diag, alpha, A, B);

    if ( parsec_ztrsm != NULL )
    {
        parsec_context_add_taskpool( parsec, parsec_ztrsm );
        dplasma_wait_until_completion( parsec );
        dplasma_ztrsm_Destruct( parsec_ztrsm );
        return 0;
    }
    else {
        return -101;
    }
}
