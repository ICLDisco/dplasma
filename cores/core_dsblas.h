/**
 * Copyright (c) 2019      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Imported from:
 *
 * @file core_dsblas.h
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
 * @generated ds Fri Apr  1 11:02:20 2016
 *
 **/
#ifndef _PLASMA_CORE_DSBLAS_H_
#define _PLASMA_CORE_DSBLAS_H_

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/
void CORE_slag2d(int m, int n,
                 const float *A, int lda,
                 double *B, int ldb);
void CORE_dlag2s(int m, int n,
                 const double *A, int lda,
                 float *B, int ldb, int *info);

#ifdef __cplusplus
}
#endif

#endif
