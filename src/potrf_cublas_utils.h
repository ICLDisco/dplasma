/*
 * Copyright (c) 2020-2023 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 */
#ifndef DPLASMA_POTRF_CUBLAS_UTILS_H
#define DPLASMA_POTRF_CUBLAS_UTILS_H

#if defined(DPLASMA_HAVE_CUDA)

typedef struct {
  char         *tmpmem;
  void         *memory;
  int           lwork;
  void*         params;
  size_t        host_size;
  void*         host_buffer;
} dplasma_potrf_workspace_t;

#endif

#endif //DPLASMA_POTRF_CUBLAS_UTILS_H
