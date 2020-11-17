#ifndef _DPLASMAJDF_LAPACK_DTT_H_
#define _DPLASMAJDF_LAPACK_DTT_H_

#include "dplasma/types.h"
#include "dplasma/types_lapack.h"

/* DON'T CHANGE SHAPE */
#define SAME      -1

/* Obtain LDA of the current datacopy.
 */
static inline
int LDA_internal(const parsec_tiled_matrix_dc_t *dc, parsec_data_copy_t *cp)
{
    int lda, rows, cols, loc, shape, layout, rc;
    parsec_arena_datatype_t *adt;
    /* obtain the lda of this dc->dtt */
    rc = dplasma_get_datatype_info(dc, cp->dtt, &lda, &rows, &cols, &loc, &shape, &layout, &adt);
    assert(rc == 0);
    PARSEC_DEBUG_VERBOSE(4, parsec_debug_output, "CP %p [%p] [type %p] lda %d", cp, cp->device_private, cp->dtt, lda);
    return lda;
}
#define LDA(dc, FLOW_NAME)\
    LDA_internal(dc->dc_original, _f_##FLOW_NAME)

/* Obtain the appropriate parsec_arena_datatype_t for a given location and shape in the datacollection
 */
static inline
parsec_arena_datatype_t* ADTT_DC(const dplasma_data_collection_t *ddc, int loc, int target_shape, int target_layout)
{
    (void) loc; (void) target_shape;

    int lda, rc;
    parsec_arena_datatype_t *adt;

    /* obtain the dtt for the location & layout with the target_shape */
    rc = dplasma_get_datatype(ddc->dc_original, loc, target_shape, target_layout, &lda, &adt);
    assert(rc == 0);

    PARSEC_DEBUG_VERBOSE(8, parsec_debug_output, "LOC %d target_shape %d target_layout %d-> dtt %p ", loc, target_shape, target_layout, adt->opaque_dtt);
    return adt;
}

/* Obtain the appropriate parsec_arena_datatype_t that represents the type of the datacopy with the given shape
 */
static parsec_arena_datatype_t adt_null_dc;
static inline
parsec_arena_datatype_t* ADTT_CP(parsec_data_copy_t *cp, const dplasma_data_collection_t *ddc, int target_loc, int target_shape)
{
    (void) cp; (void) target_shape;

    if(cp == NULL){
      /* this flow is not actually originating an output dep */
      /* this case happens because the JDF generated code obtains the ADT before
       * evaluating the guards for the successors during iterate successors.
       */
      return &adt_null_dc;
    }
    int lda, rows, cols, loc, shape, layout, rc;
    parsec_arena_datatype_t *adt;
    /* obtain the location & layout of this dc->dtt */
    rc = dplasma_get_datatype_info(ddc->dc_original, cp->dtt, &lda, &rows, &cols, &loc, &shape, &layout, &adt);
    assert(rc == 0);

    if(( shape == target_shape )||(target_shape == SAME)){
      PARSEC_DEBUG_VERBOSE(8, parsec_debug_output, "CP %p [type %p] -> target_shape %d target_loc %d dtt %p", cp, cp->dtt, target_shape, target_loc, adt->opaque_dtt);
      return adt;
    }

    /* obtain the equivalent dtt for the same location with the target_shape */
    rc = dplasma_get_datatype(ddc->dc_original, target_loc, target_shape, layout, &lda, &adt);
    assert(rc == 0);
    PARSEC_DEBUG_VERBOSE(8, parsec_debug_output, "CP %p [type %p] loc %d layout %d -> dtt %p target_shape %d", cp, cp->dtt, target_loc, layout, adt->opaque_dtt, target_shape);
    return adt;
}

/* Obtain number of rows give datacopy and datacollection */
static inline
void ADTT_INFO_internal(parsec_data_copy_t *cp, const parsec_tiled_matrix_dc_t *dc, int *lda, int *rows, int *cols)
{
    assert(cp != NULL);
    int loc, shape, layout, rc;
    parsec_arena_datatype_t *adt;
    /* obtain the location & layout of this dc->dtt */
    rc = dplasma_get_datatype_info(dc, cp->dtt, lda, rows, cols, &loc, &shape, &layout, &adt);
    assert(rc == 0);
}

#define ADTT_INFO(dc, FLOW_NAME, lda, rows, cols)\
    ADTT_INFO_internal(_f_##FLOW_NAME, dc->dc_original, lda, rows, cols)


#if defined(DPLASMA_HAVE_CUDA)
/* Use cudaMemcpy2DAsync or loop with cudaMemcpyAsync for data transfers to device */
#define CUDA_COPY_2D

/* Functions to transfer data in and out of the GPU.
 * Assuming a full tiled has been allocated on the GPU (mb*nb*size(elem))
 */
static cudaError_t
stage_in_lapack(parsec_gpu_task_t *gtask,
                uint32_t flow_mask,
                parsec_gpu_exec_stream_t *gpu_stream)
{
    cudaError_t ret;
    parsec_data_copy_t * copy_in;
    parsec_data_copy_t * copy_out;
    parsec_device_cuda_module_t *in_elem_dev;
    parsec_tiled_matrix_dc_t * dc;
    parsec_task_t *task = gtask->ec;
    int elem_sz;
    int i;
    for(i = 0; i < task->task_class->nb_flows; i++){
        if(flow_mask & (1U << i)){
            copy_in = task->data[i].data_in;
            copy_out = task->data[i].data_out;
            dc = (parsec_tiled_matrix_dc_t*)gtask->flow_dc[i];
            assert(dc != NULL);
            elem_sz = parsec_datadist_getsizeoftype(dc->mtype);
            in_elem_dev = (parsec_device_cuda_module_t*)parsec_mca_device_get( copy_in->device_index);
            if( (in_elem_dev->super.type == PARSEC_DEV_CUDA) || (dc->storage != matrix_Lapack)){
                ret = (cudaError_t)cudaMemcpyAsync( copy_out->device_private,
                                                    copy_in->device_private,
                                                    gtask->flow_nb_elts[i],
                                                    (in_elem_dev->super.type != PARSEC_DEV_CUDA)?
                                                            cudaMemcpyHostToDevice : cudaMemcpyDeviceToDevice,
                                                    gpu_stream->cuda_stream);
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return ret; } );
            }else{

#ifdef CUDA_COPY_2D
                int ldd, nrows, ncols;
                ADTT_INFO_internal(copy_in, dc, &ldd, &nrows, &ncols);
                size_t dpitch = dc->mb * elem_sz;
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
                                                      gpu_stream->cuda_stream );
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpy2DAsync ", ret, { return ret; } );


#else

                int ldd, nrows, ncols;
                ADTT_INFO_internal(copy_in, dc, &ldd, &nrows, &ncols);

                int j;
                for(j=0; j<ncols; j++) {
                    char*src = ((char*)copy_in->device_private) + j * ldd * elem_sz;
                    char*dst = ((char*)copy_out->device_private) + j * dc->mb * elem_sz;
                    ret = cudaMemcpyAsync(dst,
                                          src,
                                          nrows * elem_sz,
                                          cudaMemcpyHostToDevice,
                                          gpu_stream->cuda_stream );
                    PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return ret; } );

                }
#endif


            }
        }
    }
    return ret;
}

static cudaError_t
stage_out_lapack(parsec_gpu_task_t *gtask,
                 uint32_t flow_mask,
                 parsec_gpu_exec_stream_t *gpu_stream)
{
    cudaError_t ret;
    parsec_data_copy_t * copy_in;
    parsec_data_copy_t * copy_out;
    parsec_device_cuda_module_t *out_elem_dev;
    parsec_task_t *task = gtask->ec;
    parsec_tiled_matrix_dc_t * dc;
    int elem_sz;
    int i;
    for(i = 0; i < task->task_class->nb_flows; i++){
        if(flow_mask & (1U << i)){
            copy_in = task->data[i].data_out;
            copy_out = copy_in->original->device_copies[0];
            dc = (parsec_tiled_matrix_dc_t*)gtask->flow_dc[i];
            assert(dc != NULL);
            elem_sz = parsec_datadist_getsizeoftype(dc->mtype);
            out_elem_dev = (parsec_device_cuda_module_t*)parsec_mca_device_get( copy_out->device_index);

            if( (out_elem_dev->super.type == PARSEC_DEV_CUDA) || (dc->storage != matrix_Lapack)){
                ret = (cudaError_t)cudaMemcpyAsync( copy_out->device_private,
                                                    copy_in->device_private,
                                                    gtask->flow_nb_elts[i],
                                                    out_elem_dev->super.type != PARSEC_DEV_CUDA ?
                                                            cudaMemcpyDeviceToHost : cudaMemcpyDeviceToDevice,
                                                    gpu_stream->cuda_stream);
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return ret; } );
            }else{

#ifdef CUDA_COPY_2D
                int ldd, nrows, ncols;
                ADTT_INFO_internal(copy_out, dc, &ldd, &nrows, &ncols);
                size_t dpitch = ldd * elem_sz;
                size_t spitch = dc->mb * elem_sz;
                size_t width  = nrows * elem_sz;
                size_t height = ncols;
                /* copy width bytes heigth times, skipping pitch - width bytes every time */
                ret = (cudaError_t)cudaMemcpy2DAsync( copy_out->device_private,
                                                      dpitch, /*dst pitch bytes*/
                                                      copy_in->device_private,
                                                      spitch, /*src pitch bytes*/
                                                      width, height,
                                                      cudaMemcpyDeviceToHost,
                                                      gpu_stream->cuda_stream);
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpy2DAsync ", ret, { return ret; } );
#else
                int ldd, nrows, ncols;
                ADTT_INFO_internal(copy_out, dc, &ldd, &nrows, &ncols);
                int j;
                for(j=0; j<ncols; j++) {
                    char*src = ((char*)copy_in->device_private) + j * dc->mb * elem_sz;
                    char*dst = ((char*)copy_out->device_private) + j * ldd * elem_sz;
                    ret = cudaMemcpyAsync(dst,
                                          src,
                                          nrows * elem_sz,
                                          cudaMemcpyDeviceToHost,
                                          gpu_stream->cuda_stream);
                    PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return ret; } );
                }
#endif
            }
        }
    }
    return ret;
}
#endif /* defined(DPLASMA_HAVE_CUDA) */

#endif /* _DPLASMAJDF_LAPACK_DTT_H_ */

