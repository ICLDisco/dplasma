/*
 * Copyright (c) 2020-2023 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 */
#ifndef DPLASMA_POTRF_GPU_WORKSPACES_H
#define DPLASMA_POTRF_GPU_WORKSPACES_H

typedef struct {
  char         *tmpmem;
  void         *memory;
  int           lwork;
  void*         params;
  size_t        host_size;
  void*         host_buffer;
} dplasma_potrf_gpu_workspaces_t;

#endif //DPLASMA_POTRF_GPU_WORKSPACES_H
