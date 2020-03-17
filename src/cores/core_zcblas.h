/**
 * Copyright (c) 2019-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Imported from:
 *
 * @file core_zcblas.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 **/
#ifndef _PLASMA_CORE_ZCBLAS_H_
#define _PLASMA_CORE_ZCBLAS_H_

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/
void CORE_clag2z(int m, int n,
                 const PLASMA_Complex32_t *A, int lda,
                 PLASMA_Complex64_t *B, int ldb);
void CORE_zlag2c(int m, int n,
                 const PLASMA_Complex64_t *A, int lda,
                 PLASMA_Complex32_t *B, int ldb, int *info);

#ifdef __cplusplus
}
#endif

#endif
