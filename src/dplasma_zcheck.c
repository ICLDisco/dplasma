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
#include <math.h>
#include <lapacke.h>
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"

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
                  parsec_tiled_matrix_dc_t *A,
                  parsec_tiled_matrix_dc_t *A0 )
{
    two_dim_block_cyclic_t *twodA = (two_dim_block_cyclic_t *)A0;
    two_dim_block_cyclic_t LLt;
    int info_factorization;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double result = 0.0;
    int M = A->m;
    int N = A->n;
    double eps = LAPACKE_dlamch_work('e');
    dplasma_enum_t side;

    two_dim_block_cyclic_init(&LLt, matrix_ComplexDouble, matrix_Tile,
                              twodA->grid.rank,
                              A->mb, A->nb, M, N, 0, 0,
                              M, N, twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols, twodA->grid.ip, twodA->grid.jq);

    LLt.mat = parsec_data_allocate((size_t)LLt.super.nb_local_tiles *
                                  (size_t)LLt.super.bsiz *
                                  (size_t)parsec_datadist_getsizeoftype(LLt.super.mtype));

    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 0.,(parsec_tiled_matrix_dc_t *)&LLt );
    dplasma_zlacpy( parsec, uplo, A, (parsec_tiled_matrix_dc_t *)&LLt );

    /* Compute LL' or U'U  */
    side = (uplo == dplasmaUpper ) ? dplasmaLeft : dplasmaRight;
    dplasma_ztrmm( parsec, side, uplo, dplasmaConjTrans, dplasmaNonUnit, 1.0,
                   A, (parsec_tiled_matrix_dc_t*)&LLt);

    /* compute LL' - A or U'U - A */
    dplasma_ztradd( parsec, uplo, dplasmaNoTrans,
                    -1.0, A0, 1., (parsec_tiled_matrix_dc_t*)&LLt);

    Anorm = dplasma_zlanhe(parsec, dplasmaInfNorm, uplo, A0);
    Rnorm = dplasma_zlanhe(parsec, dplasmaInfNorm, uplo,
                           (parsec_tiled_matrix_dc_t*)&LLt);

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
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&LLt);

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
                 parsec_tiled_matrix_dc_t *A,
                 parsec_tiled_matrix_dc_t *b,
                 parsec_tiled_matrix_dc_t *x )
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
 *          A must be two_dim_block_cyclic and fully generated.
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
                  parsec_tiled_matrix_dc_t *A,
                  parsec_tiled_matrix_dc_t *Ainv )
{
    two_dim_block_cyclic_t *twodA = (two_dim_block_cyclic_t *)A;
    two_dim_block_cyclic_t Id;
    int info_solution;
    double Anorm, Ainvnorm, Rnorm;
    double eps, result;

    eps = LAPACKE_dlamch_work('e');

    two_dim_block_cyclic_init(&Id, matrix_ComplexDouble, matrix_Tile,
                               twodA->grid.rank,
                               A->mb, A->nb, A->n, A->n, 0, 0,
                               A->n, A->n, twodA->grid.rows, twodA->grid.cols, twodA->grid.krows, twodA->grid.kcols, twodA->grid.ip, twodA->grid.jq);

    Id.mat = parsec_data_allocate((size_t)Id.super.nb_local_tiles *
                                  (size_t)Id.super.bsiz *
                                  (size_t)parsec_datadist_getsizeoftype(Id.super.mtype));

    dplasma_zlaset( parsec, dplasmaUpperLower, 0., 1., (parsec_tiled_matrix_dc_t *)&Id);

    /* Id - A^-1 * A */
    dplasma_zhemm(parsec, dplasmaLeft, uplo,
                  -1., Ainv, A,
                  1., (parsec_tiled_matrix_dc_t *)&Id );

    Anorm    = dplasma_zlanhe( parsec, dplasmaOneNorm, uplo, A );
    Ainvnorm = dplasma_zlanhe( parsec, dplasmaOneNorm, uplo, Ainv );
    Rnorm    = dplasma_zlange( parsec, dplasmaOneNorm, (parsec_tiled_matrix_dc_t*)&Id );

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
    parsec_tiled_matrix_dc_destroy((parsec_tiled_matrix_dc_t*)&Id);

    return info_solution;
}
