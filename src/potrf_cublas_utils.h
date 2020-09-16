/*
 * Copyright (c) 2020-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 */
#ifndef DPLASMA_POTRF_CUBLAS_UTILS_H
#define DPLASMA_POTRF_CUBLAS_UTILS_H

#if defined(DPLASMA_HAVE_CUDA)
#include <cublas.h>
#include <cusolverDn.h>

typedef struct {
  char         *tmpmem;
  void         *memory;
  int           lwork;
} dplasma_potrf_workspace_t;

typedef cusolverStatus_t (*cublas_spotrf_v2_t) (
        cusolverDnHandle_t handle, cublasFillMode_t uplo,
        int n, float *A, int lda,
        float *Workspace, int Lwork, int *devInfo );

typedef cusolverStatus_t  (*cublas_dpotrf_v2_t) (
        cusolverDnHandle_t handle, cublasFillMode_t uplo,
        int n, double *A, int lda,
        double *Workspace, int Lwork, int *devInfo );

typedef cusolverStatus_t (*cublas_cpotrf_v2_t) (
        cusolverDnHandle_t handle, cublasFillMode_t uplo,
        int n, cuComplex *A, int lda,
        cuComplex *Workspace, int Lwork, int *devInfo );

typedef cusolverStatus_t (*cublas_zpotrf_v2_t) (
        cusolverDnHandle_t handle, cublasFillMode_t uplo,
        int n, cuDoubleComplex *A, int lda,
        cuDoubleComplex *Workspace, int Lwork, int *devInfo );

#endif

#endif //DPLASMA_POTRF_CUBLAS_UTILS_H
