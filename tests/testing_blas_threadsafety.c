/*
 * Copyright (c) 2026      NVIDIA Corporation.  All rights reserved.
 */

/*
 * Does the BLAS underneath dplasma survive being called from several threads
 * at once?
 *
 * dplasma runs one kernel per parsec worker thread, so every level-3 call is
 * concurrent with as many others as there are threads. Not every BLAS build
 * allows that. A sequential OpenBLAS packs its panels through a central
 * buffer pool in memory.c, and the allocation out of that pool is only
 * guarded when the library was built with USE_LOCKING -- which its own
 * Makefile.rule leaves off by default and turns on automatically only for
 * the threaded builds:
 *
 *   # If you want to build a single-threaded OpenBLAS, but expect to call
 *   # this from several concurrent threads in some other program, comment
 *   # this in for thread safety. (This is done automatically for
 *   # USE_THREAD=1 , and should not be necessary when USE_OPENMP=1)
 *   # USE_LOCKING = 1
 *
 * Two threads are then handed the same scratch and each returns a result
 * mixed with the other's panel. Nothing fails, nothing is reported: the
 * operands are correct, the arguments are correct, and the answer is wrong.
 * Downstream that surfaces as a residual that is occasionally too large, on
 * a run that is otherwise indistinguishable from a good one, which is a
 * miserable thing to track down from the far end.
 *
 * Every call below computes the same product from the same read-only
 * operands, so every call owes the same bits. A difference is the library
 * failing under concurrency and nothing else -- rounding cannot vary when
 * the inputs and the arguments do not.
 */

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cores/core_blas.h"

/* Big enough that the level-3 kernels pack rather than take a small-case
 * path, small enough that thousands of calls stay quick. */
static int N = 240;
static int rounds = 400;
static int nthreads = 8;

static double *A, *B, *Cin, *reference;

typedef struct {
    int id;
    int differed;       /* rounds whose result was not the reference    */
    int first_round;    /* the earliest one, or -1                      */
    int first_element;  /* where that result first left the reference   */
    double got, want;
} worker_t;

static void *hammer( void *arg )
{
    worker_t *w = (worker_t*)arg;
    double *C = malloc((size_t)N * N * sizeof(double));
    int r, i;

    if( NULL == C ) return NULL;

    for( r = 0; r < rounds; r++ ) {
        memcpy(C, Cin, (size_t)N * N * sizeof(double));
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    N, N, N, 1.0, A, N, B, N, 1.0, C, N);

        if( 0 == memcmp(C, reference, (size_t)N * N * sizeof(double)) )
            continue;

        w->differed++;
        if( w->first_round >= 0 ) continue;
        for( i = 0; i < N * N; i++ ) {
            if( C[i] == reference[i] ) continue;
            w->first_round = r;
            w->first_element = i;
            w->got = C[i];
            w->want = reference[i];
            break;
        }
    }
    free(C);
    return NULL;
}

int main( int argc, char *argv[] )
{
    pthread_t *threads;
    worker_t *workers;
    int i, failed = 0, total = 0;

    if( argc > 1 ) nthreads = atoi(argv[1]);
    if( argc > 2 ) N        = atoi(argv[2]);
    if( argc > 3 ) rounds   = atoi(argv[3]);
    if( nthreads < 2 ) {
        printf("this check needs at least two threads to mean anything\n");
        return 77;  /* ctest's conventional 'skipped' */
    }

    A         = malloc((size_t)N * N * sizeof(double));
    B         = malloc((size_t)N * N * sizeof(double));
    Cin       = malloc((size_t)N * N * sizeof(double));
    reference = malloc((size_t)N * N * sizeof(double));
    threads   = malloc((size_t)nthreads * sizeof(pthread_t));
    workers   = calloc((size_t)nthreads, sizeof(worker_t));
    if( NULL == A || NULL == B || NULL == Cin || NULL == reference ||
        NULL == threads || NULL == workers ) {
        fprintf(stderr, "out of memory\n");
        return 2;
    }

    srand(51);
    for( i = 0; i < N * N; i++ ) {
        A[i]   = (double)rand() / RAND_MAX - 0.5;
        B[i]   = (double)rand() / RAND_MAX - 0.5;
        Cin[i] = (double)rand() / RAND_MAX - 0.5;
    }

    /* The answer every thread owes, taken while nothing else is running. */
    memcpy(reference, Cin, (size_t)N * N * sizeof(double));
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                N, N, N, 1.0, A, N, B, N, 1.0, reference, N);

    printf("%d threads x %d calls of dgemm(%d,%d,%d), all computing the same "
           "product\n", nthreads, rounds, N, N, N);

    for( i = 0; i < nthreads; i++ ) {
        workers[i].id = i;
        workers[i].first_round = -1;
        pthread_create(&threads[i], NULL, hammer, &workers[i]);
    }
    for( i = 0; i < nthreads; i++ )
        pthread_join(threads[i], NULL);

    for( i = 0; i < nthreads; i++ ) {
        total += workers[i].differed;
        if( 0 == workers[i].differed ) continue;
        failed++;
        printf("  thread %d: %d of %d results differed, first in round %d at "
               "element %d, %.17g instead of %.17g\n",
               i, workers[i].differed, rounds, workers[i].first_round,
               workers[i].first_element, workers[i].got, workers[i].want);
    }

    free(A); free(B); free(Cin); free(reference); free(threads); free(workers);

    if( 0 == failed ) {
        printf("PASSED: every call returned the same bits\n");
        return 0;
    }
    printf("FAILED: %d of %d calls across %d threads returned something other\n"
           "than the single-threaded answer, from operands none of them wrote.\n"
           "This BLAS is not safe to call concurrently, and dplasma calls it\n"
           "that way from every parsec worker thread. Results will be wrong\n"
           "occasionally and silently.\n"
           "If this is a sequential OpenBLAS, it needs USE_LOCKING=1 at build\n"
           "time; the packaged pthread or OpenMP builds already have it, and\n"
           "can be held to one thread with OPENBLAS_NUM_THREADS=1.\n",
           total, nthreads * rounds, nthreads);
    return 1;
}
