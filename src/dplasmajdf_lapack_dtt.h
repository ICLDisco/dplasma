/*
 * Copyright (c) 2020-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation. All rights
 *                         reserved.
 *
 * $COPYRIGHT
 *
 */
#ifndef INCLUDE_DPLASMA_LAPACK_DTT_H
#define INCLUDE_DPLASMA_LAPACK_DTT_H

#include "dplasma/config.h"
#include <parsec.h>
#include <parsec/mca/device/device_gpu.h>
#include "dplasma/types.h"
#include "dplasma/types_lapack.h"

/* DON'T CHANGE SHAPE */
#define DPLASMA_SHAPE_SAME      -1

/* Obtain location on matrix.
 */
static inline
int LOC(const parsec_tiled_matrix_t *dM, int m, int n)
{
    return dplasma_tile_location(dM, m, n);
}

/* Obtain LDA of the current datacopy.
 */
static inline
int LDA_internal(const dplasma_data_collection_t *ddc, parsec_data_copy_t *cp)
{
    int rc;
    lapack_info_t info;
    parsec_arena_datatype_t *adt;
    (void)rc;
    /* obtain the lda of this dc->dtt */
    rc = dplasma_get_info_from_datatype(ddc, cp->dtt, &info, &adt);
    assert(rc == 0);
    PARSEC_DEBUG_VERBOSE(4, parsec_debug_output,
        "CP %p [%p] [type %p] lda %d", cp, cp->device_private, cp->dtt, info.lda);
    return info.lda;
}
#define LDA(ddc, FLOW_NAME)\
    LDA_internal(ddc, _f_##FLOW_NAME)

/* Obtain the appropriate parsec_arena_datatype_t for a given location and shape in the datacollection
 */
static inline
parsec_arena_datatype_t* ADTT_DC(const dplasma_data_collection_t *ddc, int loc, int target_shape, int target_layout)
{
    int rc;
    lapack_info_t info;
    parsec_arena_datatype_t *adt;
    (void)rc;

    info.loc = loc;
    info.shape = target_shape;
    info.layout = target_layout;

    /* obtain the dtt for the location & layout with the target_shape */
    rc = dplasma_get_datatype_from_info(ddc, &info, &adt);
    assert(rc == 0);

    PARSEC_DEBUG_VERBOSE(8, parsec_debug_output,
        "LOC %d target_shape %d target_layout %d-> dtt %p ", loc, target_shape, target_layout, adt->opaque_dtt);
    return adt;
}

/* Obtain the appropriate parsec_arena_datatype_t that represents the type of the datacopy with the given shape
 */
static parsec_arena_datatype_t adt_null_dc;
static inline
parsec_arena_datatype_t* ADTT_CP(parsec_data_copy_t *cp, const dplasma_data_collection_t *ddc, int target_loc, int target_shape)
{
    int rc;
    lapack_info_t info;
    parsec_arena_datatype_t *adt;
    (void)rc;

    if(cp == NULL){
      /* this flow is not actually originating an output dep */
      /* this case happens because the JDF generated code obtains the ADT before
       * evaluating the guards for the successors during iterate successors.
       */
      return &adt_null_dc;
    }

    /* obtain the location & layout of this dc->dtt */
    rc = dplasma_get_info_from_datatype(ddc, cp->dtt, &info, &adt);
    assert(rc == 0);

    if(( info.shape == target_shape )||(target_shape == DPLASMA_SHAPE_SAME)){
      PARSEC_DEBUG_VERBOSE(8, parsec_debug_output,
        "CP %p [type %p] -> target_shape %d target_loc %d dtt %p",
         cp, cp->dtt, target_shape, target_loc, adt->opaque_dtt);
      return adt;
    }

    info.loc = target_loc;
    info.shape = target_shape;
    /* obtain the equivalent dtt for the same location with the target_shape */
    rc = dplasma_get_datatype_from_info(ddc, &info, &adt);
    assert(rc == 0);
    PARSEC_DEBUG_VERBOSE(8, parsec_debug_output,
        "CP %p [type %p] loc %d layout %d -> dtt %p target_shape %d",
         cp, cp->dtt, target_loc, info.layout, adt->opaque_dtt, target_shape);
    return adt;
}

/* Obtain number of rows give datacopy and datacollection */
static inline
void ADTT_INFO_internal(parsec_data_copy_t *cp, const dplasma_data_collection_t *ddc, int *lda, int *rows, int *cols)
{
    int rc;
    lapack_info_t info;
    parsec_arena_datatype_t *adt;
    (void)rc;
    assert(cp != NULL);
    /* obtain the location & layout of this dc->dtt */
    rc = dplasma_get_info_from_datatype(ddc, cp->dtt, &info, &adt);
    assert(rc == 0);
    *lda = info.lda;
    *rows = info.rows;
    *cols = info.cols;
}

#define ADTT_INFO(ddc, FLOW_NAME, lda, rows, cols)\
    ADTT_INFO_internal(_f_##FLOW_NAME, ddc, lda, rows, cols)



/* Functions to transfer data in and out of the GPU.
 * Assuming a full tiled has been allocated on the GPU (mb*nb*size(elem))
 */
#if defined(DPLASMA_HAVE_CUDA)
int
dplasma_cuda_lapack_stage_in(parsec_gpu_task_t *gtask,
                uint32_t flow_mask,
                parsec_gpu_exec_stream_t *gpu_stream);

int
dplasma_cuda_lapack_stage_out(parsec_gpu_task_t *gtask,
                 uint32_t flow_mask,
                 parsec_gpu_exec_stream_t *gpu_stream);
#endif /* defined(DPLASMA_HAVE_CUDA) */
#if defined(DPLASMA_HAVE_HIP)
int
dplasma_hip_lapack_stage_in(parsec_gpu_task_t *gtask,
                uint32_t flow_mask,
                parsec_gpu_exec_stream_t *gpu_stream);

int
dplasma_hip_lapack_stage_out(parsec_gpu_task_t *gtask,
                 uint32_t flow_mask,
                 parsec_gpu_exec_stream_t *gpu_stream);
#endif /* defined(DPLASMA_HAVE_HIP) */

#endif /* INCLUDE_DPLASMA_LAPACK_DTT_H */

