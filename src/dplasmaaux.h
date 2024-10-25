/*
 * Copyright (c) 2011-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 *
 */

#ifndef _DPLASMAAUX_H_INCLUDED
#define _DPLASMAAUX_H_INCLUDED

/**
 * Returns the priority limit of a specific function for a specific precision.
 *
 * @details
 *   This auxiliary helper function uses the process environment to
 *   find what priority limit the user wants to set for a given
 *   kernel.
 *
 *   The priority limit defines how priorities should pursue the
 *   critical path (see e.g. zpotrf_wrapper.c). By convention,
 *   the environment variable name is a capital S,D,C, or Z to
 *   specify the precision concatenated with the kernel name
 *   in capital (e.g. SPOTRF to set the priority limit of potrf
 *   in simple precision).
 *
 *    @param[IN] function: the base function name used to compose
 *                         the environment variable name
 *    @param[IN] dc:       a data collection used to find the precision
 *    @return the priority limit that the user wants to define for
 *            this function and precision.
 */
int dplasma_aux_get_priority_limit( char* function, const parsec_tiled_matrix_t* dc );

/**
 * Returns the lookahead to use for GEMM
 *
 * @details
 *   This auxiliary helper function apply internal heuristics
 *   to determine the look ahead used in some SUMMA GEMM algorithms
 *   used in dplasma.
 *
 *  @param[IN] the data collection pointing to the A matrix of the
 *             GEMM operation.
 *  @return depending on the number of nodes and the matrix size,
 *          the value to use for the look ahead in SUMMA.
 */
int dplasma_aux_getGEMMLookahead( parsec_tiled_matrix_t *A );

/**
 *  Create a dplasma-specific communicator
 *
 *  @details
 *    This function allocates a dplasma-specific communicator when
 *    compiling and running in distributed using MPI. It requires
 *    PaRSEC to be compiled with MPI.
 *
 *    Allocating a dplasma-specific communicator ensures that no
 *    dplasma communication in Wrappers will interfere with an
 *    application communication. It is recommended to use in
 *    MPI setups.
 *
 *    To ensure portability, the communicator is passed as a pointer,
 *    but it should point to an actual communicator that will be dupplicated.
 *
 *    If this function is not called, dplasma will use MPI_COMM_WORLD
 *    by default, creating risks of deadlocks and errors if some
 *    dplasma communication overlaps with other application communications
 *    on the same communicator.
 *
 *    dplasma_aux_free_comm must be called before MPI_Finalize() if
 *    dplasma_aux_dup_comm is called.
 *
 *    It is incorrect to call this function twice without calling
 *    dplasma_aux_free_comm between each call.
 *
 *    @param[IN] _psrc: a pointer to a valid MPI communicator
 *    @return the error code of MPI_Comm_dup
 */
int dplasma_aux_dup_comm(void *_psrc);

/**
 *  Free the dplasma-specific communicator
 *
 *  @details
 *     This function frees the communicator allocated via dplasma_aux_dup_comm.
 *
 *     @return the error code of MPI_Comm_free
 */
int dplasma_aux_free_comm(void);

/**
 * Globally visible pointer to the dplasma-specific communicator
 */
extern void *dplasma_pcomm;

#define dplasma_wait_until_completion( object ) \
    do {                                        \
        parsec_context_start( object );         \
        parsec_context_wait( object );          \
    } while (0)


#if defined(DPLASMA_DEBUG)
#include <stdio.h>
#define dplasma_error(__func, __msg) do { fprintf(stderr, "%s: %s\n", (__func), (__msg)); *((volatile int*)0) = 42; } while(0)
#else
#define dplasma_error(__func, __msg) do { fprintf(stderr, "%s: %s\n", (__func), (__msg)); } while(0)
#endif /* defined(DPLASMA_DEBUG) */

#if defined(DPLASMA_HAVE_CUDA)
#include "dplasmaaux_cuda.h"
#endif
#if defined(DPLASMA_HAVE_HIP)
#include "dplasmaaux_hip.h"
#endif

#if defined(DPLASMA_HAVE_CUDA) || defined(DPLASMA_HAVE_HIP)
/* Advise data on device arguments */
typedef struct dplasma_advise_data_on_device_s {
    int nb_gpu_devices;
    int *gpu_device_index;
    int gpu_rows;
    int gpu_cols;
    int grid_rows;
    int grid_cols;
} dplasma_advise_data_on_device_t;

/* Find all GPUs
 * Size of dev_index: at least parsec_nb_devices
 * */
void dplasma_find_nb_gpus(int *dev_index, int *nb);

/* Get the most suitable process/gpu grid */
int dplasma_grid_calculation( int nb_process );

/* Operator 2D */
int dplasma_advise_data_on_device_ops_2D(parsec_execution_stream_t *es,
                        const parsec_tiled_matrix_t *descA,
                        void *_A, parsec_matrix_uplo_t uplo,
                        int m, int n, void *args);

/* Set advise data on device
 *
 * If op_args == NULL, use dplasma_advise_data_on_device_t by default 
 */
int dplasma_advise_data_on_device( parsec_context_t *parsec,
        parsec_matrix_uplo_t uplo,
        parsec_tiled_matrix_t *A,
        parsec_tiled_matrix_unary_op_t operation,
        void *op_args );
#endif

#endif /* _DPLASMAAUX_H_INCLUDED */
