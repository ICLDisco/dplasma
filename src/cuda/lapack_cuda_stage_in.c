/*
 * Copyright (c) 2020-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation. All rights
 *                         reserved.
 *
 */

#include "dplasma.h"
#include "dplasmajdf_lapack_dtt.h"

#if defined(DPLASMA_HAVE_CUDA)
#include <cuda.h>
#include <parsec/mca/device/cuda/device_cuda.h>

/* Use cudaMemcpy2DAsync or loop with cudaMemcpyAsync for data transfers to device */
#define USE_COPY_2D

int
dplasma_cuda_lapack_stage_in(parsec_gpu_task_t *gtask,
                uint32_t flow_mask,
                parsec_gpu_exec_stream_t *gpu_stream)
{
    cudaError_t ret;
    parsec_data_copy_t * copy_in;
    parsec_data_copy_t * copy_out;
    parsec_device_gpu_module_t *in_elem_dev;
    parsec_cuda_exec_stream_t *cuda_stream = (parsec_cuda_exec_stream_t*)gpu_stream;
    dplasma_data_collection_t * ddc;
    parsec_task_t *task = gtask->ec;
    int elem_sz;
    int i;
    for(i = 0; i < task->task_class->nb_flows; i++){
        if(flow_mask & (1U << i)){
            copy_in = task->data[i].data_in;
            copy_out = task->data[i].data_out;
            ddc = (dplasma_data_collection_t*)gtask->flow_info[i].flow_dc;
            assert(ddc != NULL);
            elem_sz = parsec_datadist_getsizeoftype(ddc->dc_original->mtype);
            in_elem_dev = (parsec_device_gpu_module_t*)parsec_mca_device_get( copy_in->device_index);
            if( (in_elem_dev->super.type == PARSEC_DEV_CUDA) || (ddc->dc_original->storage != PARSEC_MATRIX_LAPACK)){
                ret = (cudaError_t)cudaMemcpyAsync( copy_out->device_private,
                                                    copy_in->device_private,
                                                    gtask->flow_info[i].flow_span,
                                                    (in_elem_dev->super.type != PARSEC_DEV_CUDA)?
                                                            cudaMemcpyHostToDevice : cudaMemcpyDeviceToDevice,
                                                    cuda_stream->cuda_stream);
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return PARSEC_ERROR; } );
            }else{

#ifdef USE_COPY_2D
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
                                          cuda_stream->cuda_stream );
                    PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return PARSEC_ERROR; } );

                }
#endif


            }
        }
    }
    return PARSEC_SUCCESS;
}

int
dplasma_cuda_lapack_stage_out(parsec_gpu_task_t *gtask,
                 uint32_t flow_mask,
                 parsec_gpu_exec_stream_t *gpu_stream)
{
    cudaError_t ret;
    parsec_data_copy_t * copy_in;
    parsec_data_copy_t * copy_out;
    parsec_device_gpu_module_t *out_elem_dev;
    parsec_cuda_exec_stream_t *cuda_stream = (parsec_cuda_exec_stream_t*)gpu_stream;
    parsec_task_t *task = gtask->ec;
    dplasma_data_collection_t * ddc;
    int elem_sz;
    int i;
    for(i = 0; i < task->task_class->nb_flows; i++){
        if(flow_mask & (1U << i)){
            copy_in = task->data[i].data_out;
            copy_out = copy_in->original->device_copies[0];
            ddc = (dplasma_data_collection_t*)gtask->flow_info[i].flow_dc;
            assert(ddc != NULL);
            elem_sz = parsec_datadist_getsizeoftype(ddc->dc_original->mtype);
            out_elem_dev = (parsec_device_gpu_module_t*)parsec_mca_device_get( copy_out->device_index);

            if( (out_elem_dev->super.type == PARSEC_DEV_CUDA) || (ddc->dc_original->storage != PARSEC_MATRIX_LAPACK)){
                ret = (cudaError_t)cudaMemcpyAsync( copy_out->device_private,
                                                    copy_in->device_private,
                                                    gtask->flow_info[i].flow_span,
                                                    out_elem_dev->super.type != PARSEC_DEV_CUDA ?
                                                            cudaMemcpyDeviceToHost : cudaMemcpyDeviceToDevice,
                                                    cuda_stream->cuda_stream);
                PARSEC_CUDA_CHECK_ERROR( "cudaMemcpyAsync ", ret, { return PARSEC_ERROR; } );
            }else{

#ifdef USE_COPY_2D
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
