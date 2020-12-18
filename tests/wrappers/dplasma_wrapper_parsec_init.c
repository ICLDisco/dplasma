#include "common.h"


#define APP_NAME "SCALAPACK_WRAPPED_CALL"

static void parsec_init_wrapper_internal(){
  if( parsec_ctx == NULL ){

    char *var_REDIS_BLOCKSZ = getenv("REDIS_BLOCKSZ");
    if(var_REDIS_BLOCKSZ!=NULL){
        REDIS_BLOCKSZ = atoi(var_REDIS_BLOCKSZ);
    }

    char *var_REDIS_BLOCKSZ_MIN = getenv("REDIS_BLOCKSZ_MIN");
    if(var_REDIS_BLOCKSZ_MIN!=NULL){
        REDIS_BLOCKSZ_MIN = atoi(var_REDIS_BLOCKSZ_MIN);
    }

    char *var_ncores = getenv("PARSEC_WRAPPER_CORES");
    int ncores=-1;
    if(var_ncores!=NULL){
        ncores = atoi(var_ncores);
    }
    int parsec_argc = 1;

    char** parsec_argv = (char**)calloc(parsec_argc+1, sizeof(char*));
    int i;
    for(i=0; i<parsec_argc; i++){
        parsec_argv[i] = (char*)malloc(sizeof(char)*80);
    }
    parsec_argv[parsec_argc] = NULL;
    sprintf(parsec_argv[0], APP_NAME);

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef MEASURE_INTERNAL_TIMES
    SYNC_TIME_START();
#endif

    parsec_ctx = parsec_init(ncores, &parsec_argc, &parsec_argv);

#ifdef MEASURE_INTERNAL_TIMES
    SYNC_TIME_PRINT(rank, ("PaRSEC initialized\n"));
#endif

    for(i=0; i<parsec_argc; i++){
        free(parsec_argv[i]);
    }
    free(parsec_argv);

  }
}

static void parsec_fini_wrapper_internal(){

#ifdef COUNT_WRAPPED_CALLS
    PARSEC_DEBUG_VERBOSE(1, parsec_debug_output,
      " WRAPPED CALLS "
      "count_PDGEMM %d count_PDPOTRF %d count_PDTRMM %d count_PDTRSM %d  count_PDGEQRF %d count_PDGETRF_1D %d count_PDGETRF_NOPIV %d"
      "REDIS [ IN %d, OUT %d]",
       count_PDGEMM,   count_PDPOTRF,   count_PDTRMM,   count_PDTRSM,    count_PDGEQRF,   count_PDGETRF_1D,   count_PDGETRF_NOPIV,
       count_REDIS_IN, count_REDIS_OUT);
#endif

	if( parsec_ctx != NULL ){
	    parsec_fini(&parsec_ctx);
	}
}

#ifdef APPLEVEL
#define ALWAYSGENERATE_F77_BINDINGS(upper_case,\
                              lower_case,\
                              single_underscore,\
                              double_underscore,\
                              wrapper_function,\
                              signature,\
                              params)\
            void single_underscore signature { wrapper_function params; }


ALWAYSGENERATE_F77_BINDINGS (PARSEC_INIT_WRAPPER,
                           parsec_init_wrapper,
                           parsec_init_wrapper_,
                           parsec_init_wrapper__,
                           parsec_init_wrapper_internal,
                           (void), ())

ALWAYSGENERATE_F77_BINDINGS (PARSEC_FINI_WRAPPER,
                           parsec_fini_wrapper,
                           parsec_fini_wrapper_,
                           parsec_fini_wrapper__,
                           parsec_fini_wrapper_internal,
                           (void), ())
#endif

GENERATE_F77_BINDINGS (PARSEC_INIT_WRAPPER,
                           parsec_init_wrapper,
                           parsec_init_wrapper_,
                           parsec_init_wrapper__,
                           parsec_init_wrapper_internal,
                           (void), ())

GENERATE_F77_BINDINGS (PARSEC_FINI_WRAPPER,
                           parsec_fini_wrapper,
                           parsec_fini_wrapper_,
                           parsec_fini_wrapper__,
                           parsec_fini_wrapper_internal,
                           (void), ())

void parsec_init_wrapped_call(MPI_Comm comm){
    parsec_remote_dep_set_ctx(parsec_ctx, (intptr_t)comm);
}

void parsec_wrapper_devices_release_memory_(void){
    parsec_devices_release_memory();
}

void parsec_wrapper_devices_reset_load_(void){
	if( parsec_ctx != NULL ){
        parsec_devices_reset_load(parsec_ctx);
	}
}

