/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
 *
 * @precisions normal z -> s d c
 *
 */

#include "dplasma.h"
#include "dplasma/types.h"
#include "dplasma/types_lapack.h"
#include "dplasmaaux.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/mca/device/cuda/device_cuda.h"
#include "utils/dplasma_info.h"

#include "zgemm_NN.h"
#include "zgemm_NT.h"
#include "zgemm_TN.h"
#include "zgemm_TT.h"

#include "zgemm_NN_summa.h"
#include "zgemm_NT_summa.h"
#include "zgemm_TN_summa.h"
#include "zgemm_TT_summa.h"

#define MAX_SHAPES 3

#include "zgemm_NN_gpu.h"

#include "parsec/utils/mca_param.h"

static parsec_taskpool_t *
dplasma_Zgemm_New_summa(dplasma_enum_t transA, dplasma_enum_t transB,
                        dplasma_complex64_t alpha, const parsec_tiled_matrix_dc_t* A, const parsec_tiled_matrix_dc_t* B,
                        dplasma_complex64_t beta,  parsec_tiled_matrix_dc_t* C,
                        dplasma_info_t opt,
                        parsec_arena_datatype_t** adt)
{
    int P, Q, IP, JQ, m, n;
    parsec_taskpool_t *zgemm_tp;
    two_dim_block_cyclic_t *Cdist;

    P = ((two_dim_block_cyclic_t*)C)->grid.rows;
    Q = ((two_dim_block_cyclic_t*)C)->grid.cols;
    IP = ((two_dim_block_cyclic_t*)C)->grid.ip;
    JQ = ((two_dim_block_cyclic_t*)C)->grid.jq;

    dplasma_data_collection_t * ddc_A = dplasma_wrap_data_collection((parsec_tiled_matrix_dc_t*)A);
    dplasma_data_collection_t * ddc_B = dplasma_wrap_data_collection((parsec_tiled_matrix_dc_t*)B);
    dplasma_data_collection_t * ddc_C = dplasma_wrap_data_collection(C);

    m = dplasma_imax(C->mt, P);
    n = dplasma_imax(C->nt, Q);

    /* Create a copy of the C matrix to be used as a data distribution metric.
     * As it is used as a NULL value we must have a data_copy and a data associated
     * with it, so we can create them here.
     * Create the task distribution */
    Cdist = (two_dim_block_cyclic_t*)malloc(sizeof(two_dim_block_cyclic_t));

    two_dim_block_cyclic_init(
            Cdist, matrix_RealDouble, matrix_Tile,
            C->super.myrank,
            1, 1, /* Dimensions of the tiles              */
            m, n, /* Dimensions of the matrix             */
            0, 0, /* Starting points (not important here) */
            m, n, /* Dimensions of the submatrix          */
            P, Q, 1, 1, IP, JQ);
    Cdist->super.super.data_of = NULL;
    Cdist->super.super.data_of_key = NULL;

    if( dplasmaNoTrans == transA ) {
        if( dplasmaNoTrans == transB ) {
            PARSEC_DEBUG_VERBOSE(3, parsec_debug_output, "zgemm_NN_summa\n");
            parsec_zgemm_NN_summa_taskpool_t* tp;
            tp = parsec_zgemm_NN_summa_new(transA, transB, alpha, beta,
                                           ddc_A, ddc_B, ddc_C, (parsec_data_collection_t*)Cdist);
            *adt = &tp->arenas_datatypes[PARSEC_zgemm_TN_DEFAULT_ADT_IDX];
            zgemm_tp = (parsec_taskpool_t*)tp;
        } else {
            PARSEC_DEBUG_VERBOSE(3, parsec_debug_output, "zgemm_NT_summa\n");
            parsec_zgemm_NT_summa_taskpool_t* tp;
            tp = parsec_zgemm_NT_summa_new(transA, transB, alpha, beta,
                                           ddc_A, ddc_B, ddc_C, (parsec_data_collection_t*)Cdist);
            *adt = &tp->arenas_datatypes[PARSEC_zgemm_TN_DEFAULT_ADT_IDX];
            zgemm_tp = (parsec_taskpool_t*)tp;
        }
    } else {
        if( dplasmaNoTrans == transB ) {
            PARSEC_DEBUG_VERBOSE(3, parsec_debug_output, "zgemm_TN_summa\n");
            parsec_zgemm_TN_summa_taskpool_t* tp;
            tp = parsec_zgemm_TN_summa_new(transA, transB, alpha, beta,
                                           ddc_A, ddc_B, ddc_C, (parsec_data_collection_t*)Cdist);
            *adt = &tp->arenas_datatypes[PARSEC_zgemm_TN_DEFAULT_ADT_IDX];
            zgemm_tp = (parsec_taskpool_t*)tp;
        } else {
            PARSEC_DEBUG_VERBOSE(3, parsec_debug_output, "zgemm_TT_summa\n");
            parsec_zgemm_TT_summa_taskpool_t* tp;
            tp = parsec_zgemm_TT_summa_new(transA, transB, alpha, beta,
                                           ddc_A, ddc_B, ddc_C,
                                           (parsec_data_collection_t*)Cdist);
            *adt = &tp->arenas_datatypes[PARSEC_zgemm_TN_DEFAULT_ADT_IDX];
            zgemm_tp = (parsec_taskpool_t*)tp;
        }
    }

    int shape = 0;
    dplasma_setup_adtt_all_loc( ddc_A,
                                parsec_datatype_double_complex_t,
                                matrix_UpperLower/*uplo*/, 1/*diag:for matrix_Upper or matrix_Lower types*/,
                                &shape);

    dplasma_setup_adtt_all_loc( ddc_B,
                                parsec_datatype_double_complex_t,
                                matrix_UpperLower/*uplo*/, 1/*diag:for matrix_Upper or matrix_Lower types*/,
                                &shape);

    dplasma_setup_adtt_all_loc( ddc_C,
                                parsec_datatype_double_complex_t,
                                matrix_UpperLower/*uplo*/, 1/*diag:for matrix_Upper or matrix_Lower types*/,
                                &shape);

    assert(shape == MAX_SHAPES);

    (void)opt; //No user-defined options for this algorithm
    return zgemm_tp;
}

static parsec_taskpool_t *
dplasma_Zgemm_New_default(dplasma_enum_t transA, dplasma_enum_t transB,
                          dplasma_complex64_t alpha, const parsec_tiled_matrix_dc_t* A, const parsec_tiled_matrix_dc_t* B,
                          dplasma_complex64_t beta,  parsec_tiled_matrix_dc_t* C,
                          dplasma_info_t opt,
                          parsec_arena_datatype_t** adt)
{
    parsec_taskpool_t* zgemm_tp;

    dplasma_data_collection_t * ddc_A = dplasma_wrap_data_collection((parsec_tiled_matrix_dc_t*)A);
    dplasma_data_collection_t * ddc_B = dplasma_wrap_data_collection((parsec_tiled_matrix_dc_t*)B);
    dplasma_data_collection_t * ddc_C = dplasma_wrap_data_collection(C);

    if( dplasmaNoTrans == transA ) {
        if( dplasmaNoTrans == transB ) {
            parsec_zgemm_NN_taskpool_t* tp;
            tp = parsec_zgemm_NN_new(transA, transB, alpha, beta,
                                     ddc_A, ddc_B, ddc_C);
            *adt = &tp->arenas_datatypes[PARSEC_zgemm_NN_DEFAULT_ADT_IDX];
            zgemm_tp = (parsec_taskpool_t*)tp;
        } else {
            parsec_zgemm_NT_taskpool_t* tp;
            tp = parsec_zgemm_NT_new(transA, transB, alpha, beta,
                                     ddc_A, ddc_B, ddc_C);
            *adt = &tp->arenas_datatypes[PARSEC_zgemm_NT_DEFAULT_ADT_IDX];
            zgemm_tp = (parsec_taskpool_t*)tp;
        }
    } else {
        if( dplasmaNoTrans == transB ) {
            parsec_zgemm_TN_taskpool_t* tp;
            tp = parsec_zgemm_TN_new(transA, transB, alpha, beta,
                                     ddc_A, ddc_B, ddc_C);
            *adt = &tp->arenas_datatypes[PARSEC_zgemm_TN_DEFAULT_ADT_IDX];
            zgemm_tp = (parsec_taskpool_t*)tp;
        }
        else {
            parsec_zgemm_TT_taskpool_t* tp;
            tp = parsec_zgemm_TT_new(transA, transB, alpha, beta,
                                     ddc_A, ddc_B, ddc_C);
            *adt = &tp->arenas_datatypes[PARSEC_zgemm_TT_DEFAULT_ADT_IDX];
            zgemm_tp = (parsec_taskpool_t*)tp;
        }
    }

    int shape = 0;
    dplasma_setup_adtt_all_loc( ddc_A,
                                parsec_datatype_double_complex_t,
                                matrix_UpperLower/*uplo*/, 1/*diag:for matrix_Upper or matrix_Lower types*/,
                                &shape);

    dplasma_setup_adtt_all_loc( ddc_B,
                                parsec_datatype_double_complex_t,
                                matrix_UpperLower/*uplo*/, 1/*diag:for matrix_Upper or matrix_Lower types*/,
                                &shape);

    dplasma_setup_adtt_all_loc( ddc_C,
                                parsec_datatype_double_complex_t,
                                matrix_UpperLower/*uplo*/, 1/*diag:for matrix_Upper or matrix_Lower types*/,
                                &shape);

    assert(shape == MAX_SHAPES);

    (void)opt; //No user-defined options for this algorithm
    return zgemm_tp;
}

static parsec_taskpool_t*
dplasma_Zgemm_New_gpu( dplasma_enum_t transA, dplasma_enum_t transB,
                       dplasma_complex64_t alpha, const parsec_tiled_matrix_dc_t* A, const parsec_tiled_matrix_dc_t* B,
                       dplasma_complex64_t beta,  parsec_tiled_matrix_dc_t* C,
                       dplasma_info_t opt,
                       parsec_arena_datatype_t** adt)
{

    parsec_taskpool_t* zgemm_tp = NULL;
    int64_t gpu_mem_nb_blocks = -1;
    size_t gpu_mem_block_size = 0;
    size_t tile_size;
    int64_t nb_block_per_tile, nb_tile_per_gpu;
    int mt, nt, kt;
    int info_found;
    char info_value[DPLASMA_MAX_INFO_VAL];
    double vd;

    int *dev_index, nbgpu, dev;
    int u, v;
    int M, Mbound, Mlim;
    int N, Nbound, Nlim;
    int K;

    int b, c, d, p, q, look_ahead;

    if( dplasmaNoTrans != transA || dplasmaNoTrans != transB ) {
        dplasma_error("dplasma_Zgemm_New_gpu", "NoTrans for A or B not implemented yet in JDF for GPUs");
        return NULL;
    }

    nbgpu = 0;
    for(dev = 0; dev < (int)parsec_nb_devices; dev++) {
        parsec_device_module_t *device = parsec_mca_device_get(dev);
        if( PARSEC_DEV_CUDA == device->type ) {
            parsec_device_cuda_module_t *cuda_device = (parsec_device_cuda_module_t*)device;
            nbgpu++;
            if( 0 == gpu_mem_block_size )
                gpu_mem_block_size = cuda_device->mem_block_size;
            if( -1 == gpu_mem_nb_blocks || cuda_device->mem_nb_blocks < gpu_mem_nb_blocks )
                gpu_mem_nb_blocks = cuda_device->mem_nb_blocks;
        }
    }
    if(nbgpu == 0) {
        dplasma_error("dplasma_Zgemm_gpu_New", "Trying to instantiate JDF for GPUs on machine without GPUs");
        return NULL;
    }
    dev_index = (int*)malloc(nbgpu * sizeof(int));
    nbgpu= 0;
    for(dev = 0; dev < (int)parsec_nb_devices; dev++) {
        parsec_device_module_t *device = parsec_mca_device_get(dev);
        if( PARSEC_DEV_CUDA == device->type ) {
            dev_index[nbgpu++] = device->device_index;
        }
    }

    p = ((two_dim_block_cyclic_t*)C)->grid.rows;
    q = ((two_dim_block_cyclic_t*)C)->grid.cols;

    vd = 1.0; // Default percentage of available memory dedicated to this GEMM
    dplasma_info_get(opt, "DPLASMA:GEMM:GPU:mem_ratio", DPLASMA_MAX_INFO_VAL, info_value, &info_found);
    if( info_found ) {
        vd = strtod(info_value, NULL);
        if(vd <= 0.0 || vd > 1.0) {
            dplasma_error("dplasma_Zgemm_New_gpu",
                          "Invalid value for DPLASMA:GEMM:GPU:mem_ratio. Mem ratio must be real in ]0, 1]");
            goto cleanup;
        }
    }
    tile_size = A->mb*A->nb*sizeof(dplasma_complex64_t);
    nb_block_per_tile = (tile_size + gpu_mem_block_size -1 ) / gpu_mem_block_size;
    gpu_mem_nb_blocks = vd * gpu_mem_nb_blocks;
    nb_tile_per_gpu = gpu_mem_nb_blocks / nb_block_per_tile;

    mt = A->mt;
    nt = B->nt;
    kt = A->nt;

    // We find (b, c) such that b*c tiles of C fill at most 75% of the GPU memory
    // and b*p divides MT
    // and c*q divides NT
    vd = 0.75; // By default it's up to 75% of the memory to host C
    dplasma_info_get(opt, "DPLASMA:GEMM:GPU:c_ratio", DPLASMA_MAX_INFO_VAL, info_value, &info_found);
    if( info_found ) {
        vd = strtod(info_value, NULL);
        if(vd <= 0.0 || vd >= 1.0) {
            dplasma_error("dplasma_Zgemm_New_gpu",
                          "Invalid value for DPLASMA:GEMM:GPU:c_ratio. Ratio of memory dedicated to hosting tiles of C must be real in ]0, 1[");
            goto cleanup;
        }
    }
    int fact = 1;
    while( fact < mt && fact < nt && ((mt/fact) * (nt/fact)) / (p * q * nbgpu) > nb_tile_per_gpu * vd ) fact++;
    b = mt/(p*fact);
    c = nt/(q*fact);

    // Usually, look ahead is detrimental to performance when fact=1
    // and critical to performance when fact > 1
    look_ahead = 1 + (fact > 1);
    dplasma_info_get(opt, "DPLASMA:GEMM:GPU:look_ahead", DPLASMA_MAX_INFO_VAL, info_value, &info_found);
    if( info_found ) {
        look_ahead = atoi(info_value);
        if(look_ahead <= 0) {
            dplasma_error("dplasma_Zgemm_New_gpu",
                          "Invalid value for DPLASMA:GEMM:GPU:look_ahead. Look ahead must be 1 or more");
            goto cleanup;
        }
    }

    // OK, now we fill up each GPU with data from A and B
    int c_per_gpu = c / nbgpu;
    int maxd = (nb_tile_per_gpu - b*c_per_gpu)/(b+c_per_gpu) - 1;
    d = maxd < kt ? maxd : kt;

    // Now we let the user overwrite the b, c and d parameters
    dplasma_info_get(opt, "DPLASMA:GEMM:GPU:b", DPLASMA_MAX_INFO_VAL, info_value, &info_found);
    if( info_found ) {
        b = atoi(info_value);
        if(b <= 0 || b*p > A->mt) {
            dplasma_error("dplasma_Zgemm_New_gpu",
                          "Invalid value for DPLASMA:GEMM:GPU:b. b must be > 0 and b*P less or equal to A.mt");
            goto cleanup;
        }
    }
    dplasma_info_get(opt, "DPLASMA:GEMM:GPU:c", DPLASMA_MAX_INFO_VAL, info_value, &info_found);
    if( info_found ) {
        c = atoi(info_value);
        if(c <= 0 || c*q > A->nt) {
            dplasma_error("dplasma_Zgemm_New_gpu",
                          "Invalid value for DPLASMA:GEMM:GPU:c. c must be > 0 and c*Q less or equal to A.nt");
            goto cleanup;
        }
    }
    dplasma_info_get(opt, "DPLASMA:GEMM:GPU:d", DPLASMA_MAX_INFO_VAL, info_value, &info_found);
    if( info_found ) {
        d = atoi(info_value);
        if(d <= 0 || d > B->mt) {
            dplasma_error("dplasma_Zgemm_New_gpu",
                          "Invalid value for DPLASMA:GEMM:GPU:d. d must be > 0 and less or equal to B.mt");
            goto cleanup;
        }
    }

    assert(d <= B->mt);
    assert( b*p <= A->mt );
    assert( c*q <= C->nt );

    {
        parsec_zgemm_NN_gpu_taskpool_t *tp;
        tp = parsec_zgemm_NN_gpu_new(transA, transB, alpha, beta,
                                     A, B, C, b, c, d, p, q, look_ahead,
                                     nbgpu, dev_index);
        *adt = &tp->arenas_datatypes[PARSEC_zgemm_NN_gpu_DEFAULT_ADT_IDX];

        u = C->super.myrank / q;
        v = C->super.myrank % q;

        M = A->mt;
        Mbound = M / (p * b);
        Mlim = p * b * Mbound + u;
        tp->_g_xMax = Mbound + (Mlim < M) - 1;

        N = C->nt;
        Nbound = N / (c * q);
        Nlim = c * q * Nbound + v;
        tp->_g_yMax = Nbound + (Nlim < N) - 1;

        K = B->mt;
        tp->_g_zMax = (K + d - 1) / d - 1;

        *adt = &tp->arenas_datatypes[PARSEC_zgemm_TT_DEFAULT_ADT_IDX];
        zgemm_tp = (parsec_taskpool_t *) tp;

        return zgemm_tp;
    }

  cleanup:
    if(NULL != dev_index)
        free(dev_index);
    return NULL;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zgemm_New - Generates the taskpool that performs one of the following
 *  matrix-matrix operations. WARNING: The computations are not done by this call.
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = dplasmaNoTrans:   A is not transposed;
 *          = dplasmaTrans:     A is transposed;
 *          = dplasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = dplasmaNoTrans:   B is not transposed;
 *          = dplasmaTrans:     B is transposed;
 *          = dplasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A.
 *
 * @param[in] B
 *          Descriptor of the distributed matrix B.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          Descriptor of the distributed matrix C.
 *          On exit, the data described by C are overwritten by the matrix (
 *          alpha*op( A )*op( B ) + beta*C )
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_zgemm_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_zgemm
 * @sa dplasma_zgemm_Destruct
 * @sa dplasma_cgemm_New
 * @sa dplasma_dgemm_New
 * @sa dplasma_sgemm_New
 *
 ******************************************************************************/
parsec_taskpool_t*
dplasma_zgemm_New_ex( dplasma_enum_t transA, dplasma_enum_t transB,
                      dplasma_complex64_t alpha, const parsec_tiled_matrix_dc_t* A, const parsec_tiled_matrix_dc_t* B,
                      dplasma_complex64_t beta,  parsec_tiled_matrix_dc_t* C, dplasma_info_t opt)
{
    parsec_arena_datatype_t* adt = NULL;
    parsec_taskpool_t* zgemm_tp = NULL;

    /* Check input arguments */
    if ((transA != dplasmaNoTrans) && (transA != dplasmaTrans) && (transA != dplasmaConjTrans)) {
        dplasma_error("dplasma_zgemm_New", "illegal value of transA");
        return NULL /*-1*/;
    }
    if ((transB != dplasmaNoTrans) && (transB != dplasmaTrans) && (transB != dplasmaConjTrans)) {
        dplasma_error("dplasma_zgemm_New", "illegal value of transB");
        return NULL /*-2*/;
    }

    if ( C->dtype & two_dim_block_cyclic_type ) {
        int nb_gpu_devices = 0, devid;
		int p = ((two_dim_block_cyclic_t*)C)->grid.rows;
	    int q = ((two_dim_block_cyclic_t*)C)->grid.cols;
		int64_t gpu_mem_block_size = 0;
		int64_t gpu_mem_nb_blocks = -1;
        for(devid = 0; devid < (int)parsec_nb_devices; devid++) {
            parsec_device_module_t *device = parsec_mca_device_get(devid);
            if( PARSEC_DEV_CUDA == device->type ) {
				parsec_device_cuda_module_t *cuda_device = (parsec_device_cuda_module_t*)device;
                nb_gpu_devices++;
				if( 0 == gpu_mem_block_size )
				    gpu_mem_block_size = cuda_device->mem_block_size;
				if( -1 == gpu_mem_nb_blocks || cuda_device->mem_nb_blocks < gpu_mem_nb_blocks )
				    gpu_mem_nb_blocks = cuda_device->mem_nb_blocks;
            }
        }
        if(0 < nb_gpu_devices) {
            int64_t tile_size = A->mb*A->nb*sizeof(dplasma_complex64_t);
            int64_t nb_block_per_tile = (tile_size + gpu_mem_block_size -1 ) / gpu_mem_block_size;
            int64_t nb_tile_per_gpu = gpu_mem_nb_blocks / nb_block_per_tile;
            int64_t nb_active_tiles_per_gpu = C->mt * C->nt / (p*q) + dplasma_aux_getGEMMLookahead(C) * A->mt / p + dplasma_aux_getGEMMLookahead(C) * B->nt / q;
            if( (A->dtype & two_dim_block_cyclic_type) &&
                (B->dtype & two_dim_block_cyclic_type) &&
                transA == dplasmaNoTrans &&
                transB == dplasmaNoTrans &&
                (nb_active_tiles_per_gpu > 0.95* nb_tile_per_gpu) ) {
                zgemm_tp = dplasma_Zgemm_New_gpu(transA, transB, alpha, A, B, beta, C, opt, &adt);
                return zgemm_tp;
            } 
            zgemm_tp = dplasma_Zgemm_New_summa(transA, transB, alpha, A, B, beta, C, opt, &adt);
            return zgemm_tp;
        }
    }
    zgemm_tp = dplasma_Zgemm_New_default(transA, transB, alpha, A, B, beta, C, opt, &adt);
    return zgemm_tp;
}

parsec_taskpool_t*
dplasma_zgemm_New( dplasma_enum_t transA, dplasma_enum_t transB,
                   dplasma_complex64_t alpha, const parsec_tiled_matrix_dc_t* A, const parsec_tiled_matrix_dc_t* B,
                   dplasma_complex64_t beta,  parsec_tiled_matrix_dc_t* C)
{
    parsec_taskpool_t *tp;
    dplasma_info_t opt;
    dplasma_info_create(&opt);
    tp = dplasma_zgemm_New_ex(transA, transB, alpha, A, B, beta, C, opt);
    dplasma_info_free(&opt);
    return tp;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zgemm_Destruct - Free the data structure associated to an taskpool
 *  created with dplasma_zgemm_New().
 *
 *******************************************************************************
 *
 * @param[in,out] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 *******************************************************************************
 *
 * @sa dplasma_zgemm_New
 * @sa dplasma_zgemm
 *
 ******************************************************************************/
void
dplasma_zgemm_Destruct( parsec_taskpool_t *tp )
{
    parsec_zgemm_NN_taskpool_t *zgemm_tp = (parsec_zgemm_NN_taskpool_t *)tp;
    dplasma_data_collection_t *ddc_A = NULL, *ddc_B = NULL, *ddc_C = NULL;

    if( zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_NN_SUMMA ||
        zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_NT_SUMMA ||
        zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_TN_SUMMA ||
        zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_TT_SUMMA) {
        parsec_zgemm_NN_summa_taskpool_t *zgemm_summa_tp = (parsec_zgemm_NN_summa_taskpool_t *)tp;
        parsec_tiled_matrix_dc_t* Cdist = (parsec_tiled_matrix_dc_t*)zgemm_summa_tp->_g_Cdist;
        if ( NULL != Cdist ) {
            parsec_tiled_matrix_dc_destroy( Cdist );
            free( Cdist );
        }
        dplasma_clean_adtt_all_loc(zgemm_summa_tp->_g_ddescA, MAX_SHAPES);
        dplasma_clean_adtt_all_loc(zgemm_summa_tp->_g_ddescB, MAX_SHAPES);
        dplasma_clean_adtt_all_loc(zgemm_summa_tp->_g_ddescC, MAX_SHAPES);

        ddc_A = zgemm_summa_tp->_g_ddescA;
        ddc_B = zgemm_summa_tp->_g_ddescB;
        ddc_C = zgemm_summa_tp->_g_ddescC;
    } else if( zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_NN ||
               zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_NT ||
               zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_TN ||
               zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_TT) {
        dplasma_clean_adtt_all_loc(zgemm_tp->_g_ddescA, MAX_SHAPES);
        dplasma_clean_adtt_all_loc(zgemm_tp->_g_ddescB, MAX_SHAPES);
        dplasma_clean_adtt_all_loc(zgemm_tp->_g_ddescC, MAX_SHAPES);

        ddc_A = zgemm_tp->_g_ddescA;
        ddc_B = zgemm_tp->_g_ddescB;
        ddc_C = zgemm_tp->_g_ddescC;
    } else if( zgemm_tp->_g_gemm_type == DPLASMA_ZGEMM_NN_GPU ) {
        parsec_zgemm_NN_gpu_taskpool_t *zgemm_gpu_tp = (parsec_zgemm_NN_gpu_taskpool_t *)tp;
        dplasma_matrix_del2arena( &zgemm_gpu_tp->arenas_datatypes[PARSEC_zgemm_NN_gpu_DEFAULT_ADT_IDX] );
    }

    parsec_taskpool_free(tp);

    /* free the dplasma_data_collection_t, after the tp stops referring to them */
    if(NULL != ddc_A)
        dplasma_unwrap_data_collection(ddc_A);
    if(NULL != ddc_B)
        dplasma_unwrap_data_collection(ddc_B);
    if(NULL != ddc_C)
        dplasma_unwrap_data_collection(ddc_C);
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 *  dplasma_zgemm - Performs one of the following matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = dplasmaNoTrans:   A is not transposed;
 *          = dplasmaTrans:     A is transposed;
 *          = dplasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = dplasmaNoTrans:   B is not transposed;
 *          = dplasmaTrans:     B is transposed;
 *          = dplasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A.
 *
 * @param[in] B
 *          Descriptor of the distributed matrix B.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          Descriptor of the distributed matrix C.
 *          On exit, the data described by C are overwritten by the matrix (
 *          alpha*op( A )*op( B ) + beta*C )
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *
 *******************************************************************************
 *
 * @sa dplasma_zgemm_New
 * @sa dplasma_zgemm_Destruct
 * @sa dplasma_cgemm
 * @sa dplasma_dgemm
 * @sa dplasma_sgemm
 *
 ******************************************************************************/
int
dplasma_zgemm( parsec_context_t *parsec,
               dplasma_enum_t transA, dplasma_enum_t transB,
               dplasma_complex64_t alpha, const parsec_tiled_matrix_dc_t *A,
                                        const parsec_tiled_matrix_dc_t *B,
               dplasma_complex64_t beta,        parsec_tiled_matrix_dc_t *C)
{
    parsec_taskpool_t *parsec_zgemm = NULL;
    int M, N, K;
    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;

    /* Check input arguments */
    if ((transA != dplasmaNoTrans) && (transA != dplasmaTrans) && (transA != dplasmaConjTrans)) {
        dplasma_error("dplasma_zgemm", "illegal value of transA");
        return -1;
    }
    if ((transB != dplasmaNoTrans) && (transB != dplasmaTrans) && (transB != dplasmaConjTrans)) {
        dplasma_error("dplasma_zgemm", "illegal value of transB");
        return -2;
    }

    if ( transA == dplasmaNoTrans ) {
        Am  = A->m;
        An  = A->n;
        Amb = A->mb;
        Anb = A->nb;
        Ai  = A->i;
        Aj  = A->j;
    } else {
        Am  = A->n;
        An  = A->m;
        Amb = A->nb;
        Anb = A->mb;
        Ai  = A->j;
        Aj  = A->i;
    }

    if ( transB == dplasmaNoTrans ) {
        Bm  = B->m;
        Bn  = B->n;
        Bmb = B->mb;
        Bnb = B->nb;
        Bi  = B->i;
        Bj  = B->j;
    } else {
        Bm  = B->n;
        Bn  = B->m;
        Bmb = B->nb;
        Bnb = B->mb;
        Bi  = B->j;
        Bj  = B->i;
    }

    if ( (Amb != C->mb) || (Anb != Bmb) || (Bnb != C->nb) ) {
        dplasma_error("dplasma_zgemm", "tile sizes have to match");
        return -101;
    }
    if ( (Am != C->m) || (An != Bm) || (Bn != C->n) ) {
        dplasma_error("dplasma_zgemm", "sizes of matrices have to match");
        return -101;
    }
    if ( (Ai != C->i) || (Aj != Bi) || (Bj != C->j) ) {
        dplasma_error("dplasma_zgemm", "start indexes have to match");
        return -101;
    }

    M = C->m;
    N = C->n;
    K = An;

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (dplasma_complex64_t)0.0 || K == 0) && beta == (dplasma_complex64_t)1.0))
        return 0;

    parsec_zgemm = dplasma_zgemm_New(transA, transB,
                                    alpha, A, B,
                                    beta, C);

    if ( parsec_zgemm != NULL )
    {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zgemm);
        dplasma_wait_until_completion(parsec);
        dplasma_zgemm_Destruct( parsec_zgemm );
        return 0;
    }
    return -101;
}
