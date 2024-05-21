/*
 * Copyright (c) 2023-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "dplasma/config.h"

#if defined(DPLASMA_HAVE_CUDA)
#include <cublas_v2.h>
#endif  /* defined(DPLASMA_HAVE_CUDA) */

#include "dplasma_z_dtd.h"

int
parsec_core_ztrsm(parsec_execution_stream_t *es, parsec_task_t *this_task)
{
    (void)es;
    int side, uplo, trans, diag;
    int  m, n, lda, ldc;
    dplasma_complex64_t alpha;
    dplasma_complex64_t *A, *C;

    parsec_dtd_unpack_args(this_task, &side, &uplo, &trans, &diag, &m, &n,
                           &alpha, &A, &lda, &C, &ldc);

    CORE_ztrsm(side, uplo, trans, diag,
               m, n, alpha,
               A, lda,
               C, ldc);

    return PARSEC_HOOK_RETURN_DONE;
}

#if defined(DPLASMA_HAVE_CUDA)

int
parsec_core_ztrsm_cuda(parsec_device_gpu_module_t* gpu_device,
                      parsec_gpu_task_t*          gpu_task,
                      parsec_gpu_exec_stream_t*   gpu_stream)
{
    int side, uplo, trans, diag;
    int  m, n, lda, ldc;
    dplasma_complex64_t alpha;
    dplasma_complex64_t *A, *Ag;
    dplasma_complex64_t *C, *Cg;
    parsec_task_t* this_task = gpu_task->ec;
    cublasStatus_t status;
    dplasma_cuda_handles_t* handles;


    parsec_dtd_unpack_args(this_task, &side, &uplo, &trans, &diag, &m, &n,
                           &alpha, &A, &lda, &C, &ldc);

    Ag = parsec_dtd_get_dev_ptr(this_task, 0);
    Cg = parsec_dtd_get_dev_ptr(this_task, 1);

    dplasma_cublas_side(side);
    dplasma_cublas_fill(uplo);
    dplasma_cublas_op(trans);
    dplasma_cublas_diag(diag);

    handles = parsec_info_get(&gpu_stream->infos, dplasma_dtd_cuda_infoid);

#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex alphag = make_cuDoubleComplex( creal(alpha), cimag(alpha));
#else
    double alphag = alpha;
#endif

#if defined(PARSEC_DEBUG_NOISIER)
    {
        char tmp[MAX_TASK_STRLEN];
        PARSEC_DEBUG_VERBOSE(10, parsec_gpu_output_stream, "GPU[%1d]:\tEnqueue on device %s priority %d",
                             gpu_device->super.device_index, parsec_task_snprintf(tmp, MAX_TASK_STRLEN,
                             (parsec_task_t *) this_task), this_task->priority);
    }
#endif /* defined(PARSEC_DEBUG_NOISIER) */

    parsec_cuda_exec_stream_t* cuda_stream = (parsec_cuda_exec_stream_t*)gpu_stream;
    cublasSetStream( handles->cublas_handle, cuda_stream->cuda_stream );
    status = cublasZtrsm(handles->cublas_handle,
                          side, uplo, trans, diag,
                          m, n, &alphag,
                          (cuDoubleComplex*)Ag, lda,
                          (cuDoubleComplex*)Cg, ldc);

    DPLASMA_CUBLAS_CHECK_STATUS( "cublasZtrsm ", status,
                                 {return PARSEC_HOOK_RETURN_ERROR;} );

    (void)gpu_device;
    return PARSEC_HOOK_RETURN_DONE;
}

#endif /* DPLASMA_HAVE_CUDA */


parsec_task_class_t*
parsec_dtd_create_ztrsm_task_class( parsec_taskpool_t* dtd_tp, int tile_full, int devices)
{
    parsec_task_class_t* ztrsm_tc =  parsec_dtd_create_task_class( dtd_tp, "ztrsm",
                /* side    */   sizeof(int), PARSEC_VALUE,
                /* uplo    */   sizeof(int), PARSEC_VALUE,
                /* trans   */   sizeof(int), PARSEC_VALUE,
                /* diag    */   sizeof(int), PARSEC_VALUE,
                /* tempmm  */   sizeof(int), PARSEC_VALUE,
                /* nb      */   sizeof(int), PARSEC_VALUE,
                /* alpha   */   sizeof(dplasma_complex64_t), PARSEC_VALUE,
                /* A(k, k) */   PASSED_BY_REF, PARSEC_INPUT | tile_full,
                /* ldak    */   sizeof(int), PARSEC_VALUE,
                /* A(m, k) */   PASSED_BY_REF, PARSEC_INOUT | tile_full | PARSEC_AFFINITY,
                /* ldam    */   sizeof(int), PARSEC_VALUE,
                PARSEC_DTD_ARG_END );

#if defined(DPLASMA_HAVE_CUDA)
    if( devices & PARSEC_DEV_CUDA )
        parsec_dtd_task_class_add_chore(dtd_tp, ztrsm_tc, PARSEC_DEV_CUDA, parsec_core_ztrsm_cuda);
#endif

    if( devices & PARSEC_DEV_CPU )
        parsec_dtd_task_class_add_chore(dtd_tp, ztrsm_tc, PARSEC_DEV_CPU, parsec_core_ztrsm);

    return ztrsm_tc;
}
