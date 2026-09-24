/*
 * Copyright (c) 2010-2022 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "dplasma.h"
#include <math.h>
#include <lapacke.h>
#include <inttypes.h>
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/utils/debug.h"

/* The factorization and the solve have been shown to produce bit-identical
 * matrices on runs the checker disagrees about, so the checker's own five
 * parallel steps are what diverge. Hash its workspace after each one. */
static void check_fingerprint( const char *what, int rank, void *mat, size_t len )
{
    const uint8_t *p = (const uint8_t*)mat;
    uint64_t h = 14695981039346656037ULL;
    size_t i;
    for( i = 0; i < len; i++ ) {
        h ^= (uint64_t)p[i];
        h *= 1099511628211ULL;
    }
    printf("CHECKPRINT rank %d %s %zu bytes %016"PRIx64"\n", rank, what, len, h);
    fflush(stdout);
}

#define CHECK_FINGERPRINT(WHAT, DC) \
    check_fingerprint(WHAT, (DC).grid.rank, (DC).mat, \
                      (size_t)(DC).super.nb_local_tiles * (size_t)(DC).super.bsiz * \
                      (size_t)parsec_datadist_getsizeoftype((DC).super.mtype))

/* One hash per matrix says an operation went wrong; one hash per tile says
 * which task produced it. The race is timing sensitive enough that printing
 * as we go suppresses it, so collect quietly and report once at the end. */
static uint64_t check_hash( const void *buf, size_t len )
{
    const uint8_t *p = (const uint8_t*)buf;
    uint64_t h = 14695981039346656037ULL;
    size_t i;
    for( i = 0; i < len; i++ ) {
        h ^= (uint64_t)p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

/* Hash a tile the way ztrmm_RLT.jdf's TILEHASH records do: the live elements
 * only, at the tile's leading dimension. Matching conventions is the point --
 * it makes the value the data collection holds directly comparable to the
 * value the last task claims to have written into it, which is the one place
 * the in-task checksums cannot see. */
static uint64_t check_hash_tile( const void *ptr, int m, int n, int ld, size_t es )
{
    const char *base = (const char*)ptr;
    uint64_t h = 14695981039346656037ULL;
    size_t col = (size_t)m * es;
    int j;

    for( j = 0; j < n; j++ ) {
        const char *p = base + (size_t)j * ld * es;
        size_t i = 0;
        for( ; i + sizeof(uint64_t) <= col; i += sizeof(uint64_t) ) {
            uint64_t w;
            memcpy(&w, p + i, sizeof(uint64_t));
            h = (h ^ w) * 1099511628211ULL;
        }
        for( ; i < col; i++ )
            h = (h ^ (uint64_t)(unsigned char)p[i]) * 1099511628211ULL;
    }
    return h;
}

static void check_hash_tiles( parsec_matrix_block_cyclic_t *dc, const char *what,
                              uint64_t *out )
{
    parsec_data_collection_t *o = (parsec_data_collection_t*)dc;
    parsec_tiled_matrix_t *t = &dc->super;
    size_t es = (size_t)parsec_datadist_getsizeoftype(t->mtype);
    int m, n;

    for( m = 0; m < t->mt; m++ ) {
        int mm = (m == t->mt-1) ? t->m - m*t->mb : t->mb;
        for( n = 0; n < t->nt; n++ ) {
            int nn = (n == t->nt-1) ? t->n - n*t->nb : t->nb;
            uint64_t h;
            void *p;

            if( o->myrank != o->rank_of(o, m, n) ) continue;
            p = parsec_data_copy_get_ptr(parsec_data_get_copy(o->data_of(o, m, n), 0));
            h = check_hash_tile(p, mm, nn, t->mb, es);
            out[m * t->nt + n] = h;
            parsec_debug_history_add("TILEHASH r%d descB(%d,%d) %s %016"PRIx64" @%p\n",
                                     o->myrank, m, n, what, h, p);
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64_check
 *
 * check_zpotrf - Check the correctness of the Cholesky factorization computed
 * Cholesky functions with the following criteria:
 *
 *    \f[ ||L'L-A||_oo/(||A||_oo.N.eps) < 60. \f]
 *
 *  or
 *
 *    \f[ ||UU'-A||_oo/(||A||_oo.N.eps) < 60. \f]
 *
 *  where A is the original matrix, and L, or U, the result of the Cholesky
 *  factorization.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] loud
 *          The level of verbosity required.
 *
 * @param[in] uplo
 *          = dplasmaUpper: Upper triangle of A and A0 are referenced;
 *          = dplasmaLower: Lower triangle of A and A0 are referenced.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A result of the Cholesky
 *          factorization. Holds L or U. If uplo == dplasmaUpper, the only the
 *          upper part is referenced, otherwise if uplo == dplasmaLower, the
 *          lower part is referenced.
 *
 * @param[in] A0
 *          Descriptor of the original distributed matrix A before
 *          factorization. If uplo == dplasmaUpper, the only the upper part is
 *          referenced, otherwise if uplo == dplasmaLower, the lower part is
 *          referenced.
 *
 *******************************************************************************
 *
 * @return
 *          \retval 1, if the result is incorrect
 *          \retval 0, if the result is correct
 *
 ******************************************************************************/
int check_zpotrf( parsec_context_t *parsec, int loud,
                  dplasma_enum_t uplo,
                  parsec_tiled_matrix_t *A,
                  parsec_tiled_matrix_t *A0 )
{
    parsec_matrix_block_cyclic_t *twodA = (parsec_matrix_block_cyclic_t *)A0;
    parsec_matrix_block_cyclic_t LLt;
    int info_factorization;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double result = 0.0;
    int M = A->m;
    int N = A->n;
    double eps = LAPACKE_dlamch_work('e');
    dplasma_enum_t side;

    parsec_matrix_block_cyclic_init(&LLt, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                              twodA->grid.rank,
                              A->mb, A->nb, M, N, 0, 0,
                              M, N, twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols, twodA->grid.ip, twodA->grid.jq);

    LLt.mat = parsec_data_allocate((size_t)LLt.super.nb_local_tiles *
                                  (size_t)LLt.super.bsiz *
                                  (size_t)parsec_datadist_getsizeoftype(LLt.super.mtype));

    /* The trmm below has been caught producing a different answer than a
     * passing run from bit-identical A and LLt, on one rank, while the
     * standalone trmm tester never fails. Repeating it in-process says
     * whether it is racy every time or only on its first invocation after
     * the factorization taskpool tore down. */
    {
        const char *env = getenv("DPLASMA_CHECK_REPEAT");
        int repeat = (NULL == env) ? 1 : atoi(env);
        int mt = LLt.super.mt, nt = LLt.super.nt;
        uint64_t *in = NULL, *out = NULL;
        int it, m, n, caught = 0;

        if( repeat < 1 ) repeat = 1;
        side = (uplo == dplasmaUpper ) ? dplasmaLeft : dplasmaRight;

        if( repeat > 1 ) {
            in  = calloc((size_t)repeat * mt * nt, sizeof(uint64_t));
            out = calloc((size_t)repeat * mt * nt, sizeof(uint64_t));
        }

        /* Every iteration is kept, not just the one that goes wrong. A run is
         * only interesting when some iteration disagrees with the others, and
         * the useful comparison is then task by task against an iteration that
         * agreed -- which means both have to be in the same dump. Recording is
         * in-memory, so unlike printing it does not smother the race. */
        if( repeat > 1 ) parsec_debug_history_purge();

        for( it = 0; it < repeat; it++ ) {
            dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0.,(parsec_tiled_matrix_t *)&LLt );
            dplasma_zlacpy( parsec, uplo, A, (parsec_tiled_matrix_t *)&LLt );

            if( repeat > 1 ) {
                parsec_debug_history_add("=== ztrmm iteration %d on rank %d begins\n",
                                         it, LLt.grid.rank);
                check_hash_tiles(&LLt, "lacpy", in + (size_t)it * mt * nt);
            } else {
                CHECK_FINGERPRINT("LLt-after-lacpy", LLt);
            }

            /* Compute LL' or U'U  */
            dplasma_ztrmm( parsec, side, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0,
                           A, (parsec_tiled_matrix_t*)&LLt);

            if( repeat > 1 ) check_hash_tiles(&LLt, "final", out + (size_t)it * mt * nt);
            else             CHECK_FINGERPRINT("LLt-after-trmm", LLt);

            if( repeat > 1 && it > 0 && !caught ) {
                int o;
                for( o = 0; o < mt * nt; o++ ) {
                    if( out[(size_t)it * mt * nt + o] == out[o] ) continue;
                    printf("CHECKCAUGHT rank %d iteration %d tile(%d,%d) "
                           "%016"PRIx64" != %016"PRIx64"\n",
                           LLt.grid.rank, it, o / nt, o % nt,
                           out[(size_t)it * mt * nt + o], out[o]);
                    fflush(stdout);
                    caught = 1;
                    break;
                }
            }
        }

        /* Dumped by every rank, not just the ones holding a wrong tile: a bad
         * value can be produced anywhere and only surface where it lands. */
        if( repeat > 1 && NULL != getenv("DPLASMA_CHECK_TRACE") )
            parsec_debug_history_dump();

        /* Report once, after every iteration is done, so the printing cannot
         * perturb the race it is trying to observe. */
        for( m = 0; NULL != out && m < mt; m++ ) {
            for( n = 0; n < nt; n++ ) {
                size_t o = (size_t)m * nt + n;
                int differs = 0;

                for( it = 1; it < repeat; it++ )
                    if( out[(size_t)it * mt * nt + o] != out[o] ) differs = 1;
                if( !differs ) continue;

                printf("CHECKVARY rank %d tile(%d,%d) trmm output varies:", LLt.grid.rank, m, n);
                for( it = 0; it < repeat; it++ )
                    printf(" %016"PRIx64, out[(size_t)it * mt * nt + o]);
                printf("\n            input was:");
                for( it = 0; it < repeat; it++ )
                    printf(" %016"PRIx64, in[(size_t)it * mt * nt + o]);
                printf("\n");
            }
        }
        if( NULL != out ) fflush(stdout);
        free(in); free(out);
    }

    /* compute LL' - A or U'U - A */
    dplasma_ztradd( parsec, uplo, dplasmaNoTrans,
                    -1.0, A0, 1., (parsec_tiled_matrix_t*)&LLt);
    CHECK_FINGERPRINT("LLt-after-tradd", LLt);

    Anorm = dplasma_zlanhe(parsec, dplasmaInfNorm, uplo, A0);
    Rnorm = dplasma_zlanhe(parsec, dplasmaInfNorm, uplo,
                           (parsec_tiled_matrix_t*)&LLt);

    result = Rnorm / ( Anorm * N * eps ) ;

    if ( loud > 2 ) {
        printf("============\n");
        printf("Checking the Cholesky factorization \n");

        if ( loud > 3 )
            printf( "-- ||A||_oo = %e, ||L'L-A||_oo = %e\n", Anorm, Rnorm );

        printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n", result);
    }

    if ( isnan(Rnorm)  || isinf(Rnorm)  ||
         isnan(result) || isinf(result) ||
         (result > 60.0) )
    {
        if( loud ) printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else
    {
        if( loud ) printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    parsec_data_free(LLt.mat); LLt.mat = NULL;
    parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&LLt);

    return info_factorization;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64_check
 *
 * check_zaxmb - Returns the result of the following test
 *
 *    \f[ (|| A x - b ||_oo / ((||A||_oo * ||x||_oo + ||b||_oo) * N * eps) ) < 60. \f]
 *
 *  where A is the original matrix, b the original right hand side, and x the
 *  solution computed through any factorization.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] loud
 *          The level of verbosity required.
 *
 * @param[in] uplo
 *          = dplasmaUpper: Upper triangle of A is referenced;
 *          = dplasmaLower: Lower triangle of A is referenced.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A result of the Cholesky
 *          factorization. Holds L or U. If uplo == dplasmaUpper, the only the
 *          upper part is referenced, otherwise if uplo == dplasmaLower, the
 *          lower part is referenced.
 *
 * @param[in,out] b
 *          Descriptor of the original distributed right hand side b.
 *          On exit, b is overwritten by (b - A * x).
 *
 * @param[in] x
 *          Descriptor of the solution to the problem, x.
 *
 *******************************************************************************
 *
 * @return
 *          \retval 1, if the result is incorrect
 *          \retval 0, if the result is correct
 *
 ******************************************************************************/
int check_zaxmb( parsec_context_t *parsec, int loud,
                 dplasma_enum_t uplo,
                 parsec_tiled_matrix_t *A,
                 parsec_tiled_matrix_t *b,
                 parsec_tiled_matrix_t *x )
{
    int info_solution;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double Bnorm = 0.0;
    double Xnorm, result;
    int N = b->m;
    double eps = LAPACKE_dlamch_work('e');

    Anorm = dplasma_zlanhe(parsec, dplasmaInfNorm, uplo, A);
    Bnorm = dplasma_zlange(parsec, dplasmaInfNorm, b);
    Xnorm = dplasma_zlange(parsec, dplasmaInfNorm, x);

    /* Compute b - A*x */
    dplasma_zhemm( parsec, dplasmaLeft, uplo, -1.0, A, x, 1.0, b);

    Rnorm = dplasma_zlange(parsec, dplasmaInfNorm, b);

    result = Rnorm / ( ( Anorm * Xnorm + Bnorm ) * N * eps ) ;

    if ( loud > 2 ) {
        printf("============\n");
        printf("Checking the Residual of the solution \n");
        if ( loud > 3 )
            printf( "-- ||A||_oo = %e, ||X||_oo = %e, ||B||_oo= %e, ||A X - B||_oo = %e\n",
                    Anorm, Xnorm, Bnorm, Rnorm );

        printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);
    }

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        if( loud ) printf("-- Solution is suspicious ! \n");
        info_solution = 1;
    }
    else{
        if( loud ) printf("-- Solution is CORRECT ! \n");
        info_solution = 0;
    }

    return info_solution;
}


/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64_check
 *
 * check_zpoinv - Returns the result of the following test
 *
 *    \f[ (|| I - A * A^(-1) ||_one / (||A||_one * ||A^(-1)||_one * N * eps) ) < 10. \f]
 *
 *  where A is the original matrix, and Ainv the result of a cholesky inversion.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] loud
 *          The level of verbosity required.
 *
 * @param[in] uplo
 *          = dplasmaUpper: Upper triangle of A is referenced;
 *          = dplasmaLower: Lower triangle of A is referenced.
 *
 * @param[in] A
 *          Descriptor of the distributed original matrix A.
 *          A must be parsec_matrix_block_cyclic and fully generated.
 *
 * @param[in] Ainv
 *          Descriptor of the computed distributed A inverse.
 *
 *******************************************************************************
 *
 * @return
 *          \retval 1, if the result is incorrect
 *          \retval 0, if the result is correct
 *
 ******************************************************************************/
int check_zpoinv( parsec_context_t *parsec, int loud,
                  dplasma_enum_t uplo,
                  parsec_tiled_matrix_t *A,
                  parsec_tiled_matrix_t *Ainv )
{
    parsec_matrix_block_cyclic_t *twodA = (parsec_matrix_block_cyclic_t *)A;
    parsec_matrix_block_cyclic_t Id;
    int info_solution;
    double Anorm, Ainvnorm, Rnorm;
    double eps, result;

    eps = LAPACKE_dlamch_work('e');

    parsec_matrix_block_cyclic_init(&Id, PARSEC_MATRIX_COMPLEX_DOUBLE, PARSEC_MATRIX_TILE,
                               twodA->grid.rank,
                               A->mb, A->nb, A->n, A->n, 0, 0,
                               A->n, A->n, twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols, twodA->grid.ip, twodA->grid.jq);

    Id.mat = parsec_data_allocate((size_t)Id.super.nb_local_tiles *
                                  (size_t)Id.super.bsiz *
                                  (size_t)parsec_datadist_getsizeoftype(Id.super.mtype));

    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_t *)&Id);

    /* Id - A^-1 * A */
    dplasma_zhemm(parsec, dplasmaLeft, uplo,
                  -1., Ainv, A,
                  1., (parsec_tiled_matrix_t *)&Id );

    Anorm    = dplasma_zlanhe( parsec, dplasmaOneNorm, uplo, A );
    Ainvnorm = dplasma_zlanhe( parsec, dplasmaOneNorm, uplo, Ainv );
    Rnorm    = dplasma_zlange( parsec, dplasmaOneNorm, (parsec_tiled_matrix_t*)&Id );

    result = Rnorm / ( (Anorm*Ainvnorm)*A->n*eps );
    if ( loud > 2 ) {
        printf("  ||A||_one = %e, ||A^(-1)||_one = %e, ||I - A * A^(-1)||_one = %e, result = %e\n",
               Anorm, Ainvnorm, Rnorm, result);
    }

    if ( isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    parsec_data_free(Id.mat);
    parsec_tiled_matrix_destroy((parsec_tiled_matrix_t*)&Id);

    return info_solution;
}
