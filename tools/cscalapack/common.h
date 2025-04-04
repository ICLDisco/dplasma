/*
 * Copyright (c) 2009-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2010      University of Denver, Colorado.
 */

#ifndef _SCALAPACK_COMMON_H_
#define _SCALAPACK_COMMON_H_

#include "../../src/flops.h"
#include "myscalapack.h"

#ifdef DPLASMA_WRAPPER_ON

#ifdef SCALAPACK_SUP_UNDERSCORE
    #define parsec_init_wrapper_ parsec_init_wrapper
    #define parsec_fini_wrapper_ parsec_fini_wrapper
#endif

    extern void parsec_init_wrapper_();
    extern void parsec_fini_wrapper_();
    extern void parsec_wrapper_devices_release_memory_(void);
    extern void parsec_wrapper_devices_reset_load_(void);
#endif

static int i0=0, i1=1;
static double m1=-1e0, p0=0e0, p1=1e0;

#ifndef imax
#define imax(_a, _b) ( (_a) < (_b) ? (_b) : (_a) )
#define imin(_a, _b) ( (_a) > (_b) ? (_b) : (_a) )
#endif

typedef enum {
    PARAM_BLACS_CTX,
    PARAM_RANK,
    PARAM_M,
    PARAM_N,
    PARAM_K,
    PARAM_NB,
    PARAM_MB,
    PARAM_SEED,
    PARAM_VALIDATE,
    PARAM_NRHS,
    PARAM_NRUNS,
    PARAM_THREAD_MT,
    PARAMS_SIZE
} params_enum_t;

void setup_params( int params[], int argc, char* argv[] );

void scalapack_pdplrnt( double *A,
                        int m, int n,
                        int mb, int nb,
                        int myrow, int mycol,
                        int nprow, int npcol,
                        int mloc,
                        int seed );

void scalapack_pdplghe( double *A,
                        int m, int n,
                        int mb, int nb,
                        int myrow, int mycol,
                        int nprow, int npcol,
                        int mloc,
                        int seed );

#endif /* _SCALAPACK_COMMON_H_ */
