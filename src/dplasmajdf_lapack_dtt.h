#ifndef _DPLASMAJDF_LAPACK_DTT_H_
#define _DPLASMAJDF_LAPACK_DTT_H_

#include "dplasma/types.h"
#include "dplasma/types_lapack.h"

/* DON'T CHANGE SHAPE */
#define SAME      -1

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
    const lapack_info_t *info;
    const parsec_arena_datatype_t *adt;
    (void)rc;
    /* obtain the lda of this dc->dtt */
    rc = dplasma_get_info_from_datatype(ddc, cp->dtt, &info, &adt);
    assert(rc == 0);
    PARSEC_DEBUG_VERBOSE(4, parsec_debug_output,
        "CP %p [%p] [type %p] lda %d", cp, cp->device_private, cp->dtt, info->lda);
    return info->lda;
}
#define LDA(ddc, FLOW_NAME)\
    LDA_internal(ddc, _f_##FLOW_NAME)

/* Obtain the appropriate parsec_arena_datatype_t for a given location and shape in the datacollection
 */
static inline const parsec_arena_datatype_t*
ADTT_DC(const dplasma_data_collection_t *ddc, int loc, int target_shape, int target_layout)
{
    int rc;
    lapack_info_t info;
    const parsec_arena_datatype_t *adt;
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
static inline const parsec_arena_datatype_t*
ADTT_CP(parsec_data_copy_t *cp, const dplasma_data_collection_t *ddc, int target_loc, int target_shape)
{
    const parsec_arena_datatype_t *adt;
    const lapack_info_t* cp_info;
    int rc;
    (void)rc;

    if(cp == NULL) {
      /* this flow is not actually originating an output dep */
      /* this case happens because the JDF generated code obtains the ADT before
       * evaluating the guards for the successors during iterate successors.
       */
      return &adt_null_dc;
    }

    /* obtain the location & layout of this dc->dtt */
    rc = dplasma_get_info_from_datatype(ddc, cp->dtt, &cp_info, &adt);
    assert(rc == 0);

    if(( cp_info->shape == target_shape )||(target_shape == SAME)){
      PARSEC_DEBUG_VERBOSE(8, parsec_debug_output,
                           "CP %p [type %p] -> target_shape %d target_loc %d dtt %p",
                            cp, cp->dtt, target_shape, target_loc, adt->opaque_dtt);
      return adt;
    }

    /* make a copy of the info, don't alter the one currently in the hash table */
    lapack_info_t info = *cp_info;
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
    const lapack_info_t* info;
    const parsec_arena_datatype_t *adt;
    (void)rc;
    assert(cp != NULL);
    /* obtain the location & layout of this dc->dtt */
    rc = dplasma_get_info_from_datatype(ddc, cp->dtt, &info, &adt);
    assert(rc == 0);
    *lda = info->lda;
    *rows = info->rows;
    *cols = info->cols;
}

#define ADTT_INFO(ddc, FLOW_NAME, lda, rows, cols)\
    ADTT_INFO_internal(_f_##FLOW_NAME, ddc, lda, rows, cols)


#if defined(DPLASMA_HAVE_CUDA)
/* Use cudaMemcpy2DAsync or loop with cudaMemcpyAsync for data transfers to device */
#define CUDA_COPY_2D

/* Functions to transfer data in and out of the GPU.
 * Assuming a full tiled has been allocated on the GPU (mb*nb*size(elem))
 */
static int
stage_in_lapack(parsec_gpu_task_t *gtask,
                uint32_t flow_mask,
                parsec_gpu_exec_stream_t *gpu_stream)
{
    cudaError_t ret;
    parsec_data_copy_t * copy_in;
    parsec_data_copy_t * copy_out;
    parsec_device_cuda_module_t *in_elem_dev;
    parsec_cuda_exec_stream_t *cuda_stream = (parsec_cuda_exec_stream_t *)gpu_stream;
    dplasma_data_collection_t * ddc;
    parsec_task_t *task = gtask->ec;
    int elem_sz;
    int i;
    for(i = 0; i < task->task_class->nb_flows; i++){
        if(flow_mask & (1U << i)){
            copy_in = task->data[i].data_in;
            copy_out = task->data[i].data_out;
            ddc = (dplasma_data_collection_t*)gtask->flow_dc[i];
            assert(ddc != NULL);
            elem_sz = parsec_datadist_getsizeoftype(ddc->dc_original->mtype);
            in_elem_dev = (parsec_device_cuda_module_t*)parsec_mca_device_get( copy_in->device_index);
            if( (in_elem_dev->super.super.type & PARSEC_DEV_CUDA) || (ddc->dc_original->storage != PARSEC_MATRIX_LAPACK)){
                ret = (cudaError_t)cudaMemcpyAsync( copy_out->device_private,
                                                    copy_in->device_private,
                                                    gtask->flow_nb_elts[i],
                                                    (in_elem_dev->super.super.type & PARSEC_DEV_CUDA)?
                                                            cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice,
                                                    cuda_stream->cuda_stream);
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return PARSEC_ERROR; } );
            }else{

#ifdef CUDA_COPY_2D
                int ldd, nrows, ncols;
                ADTT_INFO_internal(copy_in, ddc, &ldd, &nrows, &ncols);
                size_t dpitch = ddc->dc_original->mb * elem_sz;
                size_t spitch = ldd * elem_sz;
                size_t width  = nrows * elem_sz;
                size_t height = ncols;
                /* copy width bytes heigth times, skipping pitch - width bytes every time */
                ret = (cudaError_t)cudaMemcpy2DAsync( copy_out->device_private,
                                                      dpitch, /*dst pitch bytes*/
                                                      copy_in->device_private,
                                                      spitch, /*src pitch bytes*/
                                                      width, height,
                                                      cudaMemcpyHostToDevice,
                                                      cuda_stream->cuda_stream );
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpy2DAsync ", ret, { return PARSEC_ERROR; } );


#else

                int ldd, nrows, ncols;
                ADTT_INFO_internal(copy_in, ddc, &ldd, &nrows, &ncols);

                int j;
                for(j=0; j<ncols; j++) {
                    char*src = ((char*)copy_in->device_private) + j * ldd * elem_sz;
                    char*dst = ((char*)copy_out->device_private) + j * ddc->dc_original->mb * elem_sz;
                    ret = cudaMemcpyAsync(dst,
                                          src,
                                          nrows * elem_sz,
                                          cudaMemcpyHostToDevice,
                                          gpu_stream->cuda_stream );
                    PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return PARSEC_ERROR; } );

                }
#endif


            }
        }
    }
    return PARSEC_SUCCESS;
}

static int
stage_out_lapack(parsec_gpu_task_t *gtask,
                 uint32_t flow_mask,
                 parsec_gpu_exec_stream_t *gpu_stream)
{
    cudaError_t ret;
    parsec_cuda_exec_stream_t *cuda_stream = (parsec_cuda_exec_stream_t*)gpu_stream;
    parsec_data_copy_t * copy_in;
    parsec_data_copy_t * copy_out;
    parsec_device_cuda_module_t *out_elem_dev;
    parsec_task_t *task = gtask->ec;
    dplasma_data_collection_t * ddc;
    int elem_sz;
    int i;
    for(i = 0; i < task->task_class->nb_flows; i++){
        if(flow_mask & (1U << i)){
            copy_in = task->data[i].data_out;
            copy_out = copy_in->original->device_copies[0];
            ddc = (dplasma_data_collection_t*)gtask->flow_dc[i];
            assert(ddc != NULL);
            elem_sz = parsec_datadist_getsizeoftype(ddc->dc_original->mtype);
            out_elem_dev = (parsec_device_cuda_module_t*)parsec_mca_device_get( copy_out->device_index);

            if( (out_elem_dev->super.super.type & PARSEC_DEV_CUDA) || (ddc->dc_original->storage != PARSEC_MATRIX_LAPACK)){
                ret = (cudaError_t)cudaMemcpyAsync( copy_out->device_private,
                                                    copy_in->device_private,
                                                    gtask->flow_nb_elts[i],
                                                    out_elem_dev->super.super.type & PARSEC_DEV_CUDA ?
                                                            cudaMemcpyDeviceToDevice : cudaMemcpyDeviceToHost,
                                                    cuda_stream->cuda_stream);
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return PARSEC_ERROR; } );
            }else{

#ifdef CUDA_COPY_2D
                int ldd, nrows, ncols;
                ADTT_INFO_internal(copy_out, ddc, &ldd, &nrows, &ncols);
                size_t dpitch = ldd * elem_sz;
                size_t spitch = ddc->dc_original->mb * elem_sz;
                size_t width  = nrows * elem_sz;
                size_t height = ncols;
                /* copy width bytes heigth times, skipping pitch - width bytes every time */
                ret = (cudaError_t)cudaMemcpy2DAsync( copy_out->device_private,
                                                      dpitch, /*dst pitch bytes*/
                                                      copy_in->device_private,
                                                      spitch, /*src pitch bytes*/
                                                      width, height,
                                                      cudaMemcpyDeviceToHost,
                                                      cuda_stream->cuda_stream);
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpy2DAsync ", ret, { return PARSEC_ERROR; } );
#else
                int ldd, nrows, ncols;
                ADTT_INFO_internal(copy_out, ddc, &ldd, &nrows, &ncols);
                int j;
                for(j=0; j<ncols; j++) {
                    char*src = ((char*)copy_in->device_private) + j * ddc->dc_original->mb * elem_sz;
                    char*dst = ((char*)copy_out->device_private) + j * ldd * elem_sz;
                    ret = cudaMemcpyAsync(dst,
                                          src,
                                          nrows * elem_sz,
                                          cudaMemcpyDeviceToHost,
                                          cuda_stream->cuda_stream);
                    PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return PARSEC_ERROR; } );
                }
#endif
            }
        }
    }
    return PARSEC_SUCCESS;
}
#endif /* defined(DPLASMA_HAVE_CUDA) */

#endif /* _DPLASMAJDF_LAPACK_DTT_H_ */

