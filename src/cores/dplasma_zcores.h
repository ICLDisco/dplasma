/*
 * Copyright (c) 2011-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#ifndef _DPLASMA_Z_CORES_H
#define _DPLASMA_Z_CORES_H

#include "cores/core_blas.h"

int blgchase_ztrdv2(int NT, int N, int NB,
                   PLASMA_Complex64_t *A1, PLASMA_Complex64_t *A2,
                   PLASMA_Complex64_t *V1, PLASMA_Complex64_t *TAU1,
                   PLASMA_Complex64_t *V2, PLASMA_Complex64_t *TAU2,
                   int sweep, int id, int blktile);



#if defined(DPLASMA_HAVE_CUDA)
#include <cuda.h>
#include <cuda_runtime_api.h>

int dplasma_cuda_zparfb(PLASMA_enum side, PLASMA_enum trans,
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
                        cudaStream_t stream);

int dplasma_cuda_ztsmqr( PLASMA_enum side, PLASMA_enum trans,
                         int M1, int N1,
                         int M2, int N2,
                         int K, int IB,
                         PLASMA_Complex64_t *A1, int LDA1,
                         PLASMA_Complex64_t *A2, int LDA2,
                         const PLASMA_Complex64_t *V, int LDV,
                         const PLASMA_Complex64_t *T, int LDT,
                         PLASMA_Complex64_t *WORK, int LDWORK,
                         PLASMA_Complex64_t *WORKC, int LDWORKC,
                         cudaStream_t stream);
#endif /* defined(DPLASMA_HAVE_CUDA) */

#endif /* _DPLASMA_Z_CORES_ */
