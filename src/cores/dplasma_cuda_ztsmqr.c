/**
 *
 * @file dplasma_cuda_ztsmqr.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.7.1
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * Copyright (c) 2026      NVIDIA Corporation.  All rights reserved.
 **/
/*
 * @precisions normal z -> c d s
 */
#include <assert.h>
#include <stdio.h>

#include <cuda_runtime_api.h>

#include "common.h"
#include "dplasmaaux_cuda.h"

int
dplasma_cuda_zparfb(PLASMA_enum side, PLASMA_enum trans,
                    PLASMA_enum direct, PLASMA_enum storev,
                    int M1, int N1,
                    int M2, int N2,
                    int K, int L,
                    PLASMA_Complex64_t *A1, int LDA1,
                    PLASMA_Complex64_t *A2, int LDA2,
                    const PLASMA_Complex64_t *V, int LDV,
                    const PLASMA_Complex64_t *T, int LDT,
                    PLASMA_Complex64_t *WORK, int LDWORK,
                    PLASMA_Complex64_t *WORKC, int LDWORKC,
                    cublasHandle_t handle,
                    cudaStream_t stream)
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex zzero = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex zone  = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex mzone = make_cuDoubleComplex(-1.0, 0.0);
#else
    double zzero = 0.0;
    double zone  = 1.0;
    double mzone = -1.0;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    cublasStatus_t cublas_status;
    cudaError_t cuda_status;
    int j;
    (void)L;

    cublas_status = cublasSetStream(handle, stream);
    DPLASMA_CUBLAS_CHECK_STATUS("cublasSetStream ", cublas_status,
                                { return PLASMA_ERR_UNEXPECTED; });

    /* Check input arguments */
    if ((side != PlasmaLeft) && (side != PlasmaRight)) {
        return -1;
    }
    if ((trans != PlasmaNoTrans) && (trans != PlasmaConjTrans)) {
        return -2;
    }
    if ((direct != PlasmaForward) && (direct != PlasmaBackward)) {
        return -3;
    }
    if ((storev != PlasmaColumnwise) && (storev != PlasmaRowwise)) {
        return -4;
    }
    if (M1 < 0) {
        return -5;
    }
    if (N1 < 0) {
        return -6;
    }
    if ((M2 < 0) ||
        ( (side == PlasmaRight) && (M1 != M2) ) ) {
        return -7;
    }
    if ((N2 < 0) ||
        ( (side == PlasmaLeft) && (N1 != N2) ) ) {
        return -8;
    }
    if (K < 0) {
        return -9;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0))
        return PLASMA_SUCCESS;

    if (direct == PlasmaForward) {

        if (side == PlasmaLeft) {

            /*
             * Column or Rowwise / Forward / Left
             * ----------------------------------
             *
             * Form  H * A  or  H' * A  where  A = ( A1 )
             *                                     ( A2 )
             */

            /*
             * W = A1 + V' * A2:
             *      W = A1
             *      W = W + V' * A2
             *
             */
            cuda_status = cudaMemcpy2DAsync(WORK, LDWORK * sizeof(PLASMA_Complex64_t),
                                            A1,   LDA1   * sizeof(PLASMA_Complex64_t),
                                            K * sizeof(PLASMA_Complex64_t), N1,
                                            cudaMemcpyDeviceToDevice, stream);
            PARSEC_CUDA_CHECK_ERROR("cudaMemcpy2DAsync ", cuda_status,
                                    { return PLASMA_ERR_UNEXPECTED; });

            cublas_status = cublasZgemm_v2(handle,
                                           dplasma_cublas_op(PlasmaConjTrans), CUBLAS_OP_N,
                                           K, N1, M2,
                                           &zone,
                                           (const cuDoubleComplex*)V     /* K*M2  */, LDV,
                                           (const cuDoubleComplex*)A2    /* M2*N1 */, LDA2,
                                           &zone,
                                           (cuDoubleComplex*)WORK        /* K*N1  */, LDWORK);
            DPLASMA_CUBLAS_CHECK_STATUS("cublasZgemm_v2 ", cublas_status,
                                        { return PLASMA_ERR_UNEXPECTED; });

            if (WORKC == NULL) {
                /* W = op(T) * W */
                cublas_status = cublasZtrmm_v2(handle,
                                               CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER,
                                               dplasma_cublas_op(trans), CUBLAS_DIAG_NON_UNIT,
                                               K, N2,
                                               &zone,
                                               (const cuDoubleComplex*)T, LDT,
                                               (const cuDoubleComplex*)WORK, LDWORK,
                                               (cuDoubleComplex*)WORK, LDWORK);
                DPLASMA_CUBLAS_CHECK_STATUS("cublasZtrmm_v2 ", cublas_status,
                                            { return PLASMA_ERR_UNEXPECTED; });

                /* A1 = A1 - W = A1 - op(T) * W */
                for(j = 0; j < N1; j++) {
                    cublas_status = cublasZaxpy_v2(handle, K, &mzone,
                                                   (const cuDoubleComplex*)(WORK + LDWORK*j), 1,
                                                   (cuDoubleComplex*)(A1 + LDA1*j), 1);
                    DPLASMA_CUBLAS_CHECK_STATUS("cublasZaxpy_v2 ", cublas_status,
                                                { return PLASMA_ERR_UNEXPECTED; });
                }

                /* A2 = A2 - op(V) * W  */
                cublas_status = cublasZgemm_v2(handle,
                                               CUBLAS_OP_N, CUBLAS_OP_N,
                                               M2, N2, K,
                                               &mzone,
                                               (const cuDoubleComplex*)V    /* M2*K  */, LDV,
                                               (const cuDoubleComplex*)WORK /* K*N2  */, LDWORK,
                                               &zone,
                                               (cuDoubleComplex*)A2         /* m2*N2 */, LDA2);
                DPLASMA_CUBLAS_CHECK_STATUS("cublasZgemm_v2 ", cublas_status,
                                            { return PLASMA_ERR_UNEXPECTED; });

            } else {
                /* Wc = V * op(T) */
                cublas_status = cublasZgemm_v2(handle,
                                               CUBLAS_OP_N, dplasma_cublas_op(trans),
                                               M2, K, K,
                                               &zone,  (const cuDoubleComplex*)V,     LDV,
                                                       (const cuDoubleComplex*)T,     LDT,
                                               &zzero, (cuDoubleComplex*)WORKC, LDWORKC);
                DPLASMA_CUBLAS_CHECK_STATUS("cublasZgemm_v2 ", cublas_status,
                                            { return PLASMA_ERR_UNEXPECTED; });

                /* A1 = A1 - opt(T) * W */
                cublas_status = cublasZgemm_v2(handle,
                                               dplasma_cublas_op(trans), CUBLAS_OP_N,
                                               K, N1, K,
                                               &mzone, (const cuDoubleComplex*)T,    LDT,
                                                       (const cuDoubleComplex*)WORK, LDWORK,
                                               &zone,  (cuDoubleComplex*)A1,   LDA1);
                DPLASMA_CUBLAS_CHECK_STATUS("cublasZgemm_v2 ", cublas_status,
                                            { return PLASMA_ERR_UNEXPECTED; });

                /* A2 = A2 - Wc * W */
                cublas_status = cublasZgemm_v2(handle,
                                               CUBLAS_OP_N, CUBLAS_OP_N,
                                               M2, N2, K,
                                               &mzone, (const cuDoubleComplex*)WORKC, LDWORKC,
                                                       (const cuDoubleComplex*)WORK,  LDWORK,
                                               &zone,  (cuDoubleComplex*)A2,    LDA2);
                DPLASMA_CUBLAS_CHECK_STATUS("cublasZgemm_v2 ", cublas_status,
                                            { return PLASMA_ERR_UNEXPECTED; });
            }
        }
        else {
            /*
             * Column or Rowwise / Forward / Right
             * -----------------------------------
             *
             * Form  H * A  or  H' * A  where A  = ( A1 A2 )
             *
             */
            fprintf(stderr, "Not implemented (Column or Rowwise / Forward / Right)");
            return PLASMA_ERR_NOT_SUPPORTED;
        }
    }
    else {
        fprintf(stderr, "Not implemented (Backward / Left or Right)");
        return PLASMA_ERR_NOT_SUPPORTED;
    }

    return PLASMA_SUCCESS;
}

int
dplasma_cuda_ztsmqr( PLASMA_enum side, PLASMA_enum trans,
                     int M1, int N1,
                     int M2, int N2,
                     int K, int IB,
                     PLASMA_Complex64_t *A1, int LDA1,
                     PLASMA_Complex64_t *A2, int LDA2,
                     const PLASMA_Complex64_t *V, int LDV,
                     const PLASMA_Complex64_t *T, int LDT,
                     PLASMA_Complex64_t *WORK, int LDWORK,
                     PLASMA_Complex64_t *WORKC, int LDWORKC,
                     cublasHandle_t handle,
                     cudaStream_t stream)
{
    int i, i1, i3;
    int NQ, NW;
    int kb;
    int ic = 0;
    int jc = 0;
    int mi = M1;
    int ni = N1;
    int rc;

    /* Check input arguments */
    if ((side != PlasmaLeft) && (side != PlasmaRight)) {
        return -1;
    }

    /* NQ is the order of Q */
    if (side == PlasmaLeft) {
        NQ = M2;
        NW = IB;
    }
    else {
        NQ = N2;
        NW = M1;
    }

    if ((trans != PlasmaNoTrans) && (trans != PlasmaConjTrans)) {
        return -2;
    }
    if (M1 < 0) {
        return -3;
    }
    if (N1 < 0) {
        return -4;
    }
    if ( (M2 < 0) ||
         ( (M2 != M1) && (side == PlasmaRight) ) ){
        return -5;
    }
    if ( (N2 < 0) ||
         ( (N2 != N1) && (side == PlasmaLeft) ) ){
        return -6;
    }
    if ((K < 0) ||
        ( (side == PlasmaLeft)  && (K > M1) ) ||
        ( (side == PlasmaRight) && (K > N1) ) ) {
        return -7;
    }
    if (IB < 0) {
        return -8;
    }
    if (LDA1 < max(1,M1)){
        return -10;
    }
    if (LDA2 < max(1,M2)){
        return -12;
    }
    if (LDV < max(1,NQ)){
        return -14;
    }
    if (LDT < max(1,IB)){
        return -16;
    }
    if (LDWORK < max(1,NW)){
        return -18;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    if (((side == PlasmaLeft)  && (trans != PlasmaNoTrans))
        || ((side == PlasmaRight) && (trans == PlasmaNoTrans))) {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ((K-1) / IB)*IB;
        i3 = -IB;
    }

    for(i = i1; (i > -1) && (i < K); i += i3) {
        kb = min(IB, K-i);

        if (side == PlasmaLeft) {
            /*
             * H or H' is applied to C(i:m,1:n)
             */
            mi = M1 - i;
            ic = i;
        }
        else {
            /*
             * H or H' is applied to C(1:m,i:n)
             */
            ni = N1 - i;
            jc = i;
        }
        /*
         * Apply H or H' (NOTE: CORE_zparfb used to be CORE_ztsrfb)
         */
        rc = dplasma_cuda_zparfb( side, trans, PlasmaForward, PlasmaColumnwise,
                                  mi, ni, M2, N2, kb, 0,
                                  A1 + LDA1*jc+ic, LDA1,
                                  A2, LDA2,
                                  V + LDV*i, LDV,
                                  T + LDT*i, LDT,
                                  WORK, LDWORK, WORKC, LDWORKC,
                                  handle, stream );
        if (PLASMA_SUCCESS != rc) {
            return rc;
        }
    }
    return PLASMA_SUCCESS;
}
