/*
 * Copyright (c) 2011-2024 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
 *
 */

#include "dplasma.h"
#include "parsec/vpmap.h"
#include <math.h>
#include <alloca.h>
#include <string.h>
#include "dplasmaaux.h"
#include "parsec/utils/show_help.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

#if defined(PARSEC_HAVE_MPI)
/*
 * dplasma falls back to MPI_COMM_WORLD by default.
 * This is sub-optimal, as this does not provide an opportunity
 * to insulate dplasma communications from potential others, but
 * this allows to maintain the behavior that dplasma does not
 * need initialization / finalization.
 *
 * The dplasma API provides two functions to provide and free
 * a dplasma-specific communicator if needed (these should be called
 * before any other dplasma call and after any dplasma call, respectively)
 */

static MPI_Comm dplasma_comm = MPI_COMM_WORLD;
void *dplasma_pcomm = &dplasma_comm;

int dplasma_aux_dup_comm(void *_psrc)
{
    MPI_Comm *src = (MPI_Comm*)_psrc;
    return MPI_Comm_dup(*src, &dplasma_comm);
}

int dplasma_aux_free_comm(void)
{
    return MPI_Comm_free(&dplasma_comm);
}
#else
void *dplasma_pcomm = NULL;
int dplasma_aux_dup_comm(void *comm)
{
    return -1;
}

int dplasma_aux_free_comm(void)
{
    return -1;
}
#endif


int
dplasma_aux_get_priority_limit( char* function, const parsec_tiled_matrix_t* dc )
{
    char *v;
    char *keyword;

    if( NULL == function || NULL == dc )
        return 0;

    keyword = alloca( strlen(function)+2 );

    switch( dc->mtype ) {
    case PARSEC_MATRIX_FLOAT:
        sprintf(keyword, "S%s", function);
        break;
    case PARSEC_MATRIX_DOUBLE:
        sprintf(keyword, "D%s", function);
        break;
    case PARSEC_MATRIX_COMPLEX_FLOAT:
        sprintf(keyword, "C%s", function);
        break;
    case PARSEC_MATRIX_COMPLEX_DOUBLE:
        sprintf(keyword, "Z%s", function);
        break;
    default:
        return 0;
    }

    if( (v = getenv(keyword)) != NULL ) {
        return atoi(v);
    }
    return 0;
}

int
dplasma_aux_getGEMMLookahead( parsec_tiled_matrix_t *A )
{
    /**
     * Assume that the number of threads per node is constant, and compute the
     * look ahead based on the global information to get the same one on all
     * nodes.
     */
    int nbunits = vpmap_get_nb_total_threads() * A->super.nodes;
    double alpha =  3. * (double)nbunits / ( A->mt * A->nt );

    if ( A->super.nodes == 1 ) {
        /* No look ahaead */
        return dplasma_imax( A->mt, A->nt );
    }
    else {
        /* Look ahaed of at least 2, and that provides 3 tiles per computational units */
        return dplasma_imax( ceil( alpha ), 2 );
    }
}

#if defined(DPLASMA_HAVE_CUDA) || defined(DPLASMA_HAVE_HIP)

/** Find all GPUs
 * Size of dev_index: at least parsec_nb_devices
 */
void dplasma_find_nb_gpus(int *dev_index, int *nb) {
    *nb = 0;
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( PARSEC_DEV_CUDA & device->type || PARSEC_DEV_HIP & device->type ) {
            dev_index[(*nb)++] = device->device_index;
        }
    }

#if defined(DPLASMA_DEBUG)
    if((*nb) == 0) {
        char hostname[256];
        gethostname(hostname, 256);
        parsec_warning(stderr, "No CUDA device found on rank %d on %s\n",
                parsec->my_rank, hostname);
    }
#endif
}

/** Get the most suitable process/gpu grid */
int dplasma_grid_calculation( int nb_process ) {
    int P;
    for( P = (int)(sqrt(nb_process + 1.0)); P > 0; P-- ) {
        if( 0 == nb_process % P ) break;
    }
    return P;
}

/* Operator 2D */
int dplasma_advise_data_on_device_ops_2D(parsec_execution_stream_t *es,
                        const parsec_tiled_matrix_t *A,
                        void *_A, parsec_matrix_uplo_t uplo,
                        int m, int n, void *op_args) {
	dplasma_advise_data_on_device_t *args = (dplasma_advise_data_on_device_t *)op_args;

    if( args->nb_gpu_devices > 0 ) {
        /* Nested 2D grid on GPU */
        int g = (m / args->grid_rows % args->gpu_rows) * args->gpu_cols + n / args->grid_cols % args->gpu_cols;
        parsec_advise_data_on_device(A->super.data_of((parsec_data_collection_t*)A, m, n), 
                                    args->gpu_device_index[g],   
                                    PARSEC_DEV_DATA_ADVICE_PREFERRED_DEVICE );
    }

    (void)es; (void)uplo;
    return 0;
}

/* Set advise data on device
 *
 * If op_args == NULL, use dplasma_advise_data_on_device_t by default 
 */
int dplasma_advise_data_on_device(parsec_context_t *parsec,
		parsec_matrix_uplo_t uplo,
		parsec_tiled_matrix_t *A,
		parsec_tiled_matrix_unary_op_t operation,
		void *op_args) {

	if(NULL != op_args) {
		parsec_apply(parsec, uplo, A, operation, op_args);
	} else {
		/* Find the number of GPUs */
		dplasma_advise_data_on_device_t *args = (dplasma_advise_data_on_device_t *)malloc(sizeof(dplasma_advise_data_on_device_t));
        args->gpu_device_index = (int *)malloc(parsec_nb_devices * sizeof(int));
		dplasma_find_nb_gpus(args->gpu_device_index, &args->nb_gpu_devices);

        /* Calculate the nested grid for the multiple GPUs on one process
         * gpu_rows >= gpu_cols and as square as possible */
        if(dplasmaUpper == uplo) {
            args->gpu_rows = dplasma_grid_calculation(args->nb_gpu_devices);
            args->gpu_cols = args->nb_gpu_devices/args->gpu_rows;
        } else {
            args->gpu_cols = dplasma_grid_calculation(args->nb_gpu_devices);
            args->gpu_rows = args->nb_gpu_devices/args->gpu_cols;
        }

        if(dplasmaUpper == uplo || dplasmaLower == uplo) {
            args->grid_rows = ((parsec_matrix_sym_block_cyclic_t *)A)->grid.rows;
            args->grid_cols = ((parsec_matrix_sym_block_cyclic_t *)A)->grid.cols;
        } else if(dplasmaUpperLower == uplo) {
            args->grid_rows = ((parsec_matrix_block_cyclic_t *)A)->grid.rows;
            args->grid_cols = ((parsec_matrix_block_cyclic_t *)A)->grid.cols;
        } else {
            dplasma_error("dplasma_advise_data_on_device", "illegal value of uplo");
        }

#if defined(DPLASMA_DEBUG)
        parsec_warning("nb_gpu_devices %d gpu_rows %d gpu_cols %d grid_rows %d grid_cols %d\n",
                args->nb_gpu_devices, args->gpu_rows, args->gpu_cols, args->grid_rows, args->grid_cols);
#endif

        parsec_apply(parsec, uplo, A, operation, (void *)args);

        free(args->gpu_device_index);
        free(args);
    }

    return 0;
}

#endif
