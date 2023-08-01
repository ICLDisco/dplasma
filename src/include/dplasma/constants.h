/*
 * Copyright (c) 2020      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * DPLASMA API constants and enumerated types
 *
 * imported from:
 *
 * @file plasmatypes.h
 *
 *  PLASMA types header
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * Contains all common types to libcoreblas and libplasma.
 *
 **/
#ifndef _DPLASMA_CONSTANTS_H_
#define _DPLASMA_CONSTANTS_H_


/** ****************************************************************************
 *  PLASMA constants - success & error codes
 **/
#define DPLASMA_SUCCESS                 0
#define DPLASMA_ERR_NOT_INITIALIZED  -101
#define DPLASMA_ERR_REINITIALIZED    -102
#define DPLASMA_ERR_NOT_SUPPORTED    -103
#define DPLASMA_ERR_ILLEGAL_VALUE    -104
#define DPLASMA_ERR_NOT_FOUND        -105
#define DPLASMA_ERR_OUT_OF_RESOURCES -106
#define DPLASMA_ERR_INTERNAL_LIMIT   -107
#define DPLASMA_ERR_UNALLOCATED      -108
#define DPLASMA_ERR_FILESYSTEM       -109
#define DPLASMA_ERR_UNEXPECTED       -110
#define DPLASMA_ERR_SEQUENCE_FLUSHED -111

/** ****************************************************************************
 *  PLASMA types
 **/
typedef int  dplasma_enum_t;
typedef int dplasma_bool_t;
typedef long dplasma_index_t;
typedef long dplasma_size_t;

/** ****************************************************************************
 *  PLASMA constants - precisions
 **/
#define dplasmaByte          0
#define dplasmaInteger       1
#define dplasmaRealFloat     2
#define dplasmaRealDouble    3
#define dplasmaComplexFloat  4
#define dplasmaComplexDouble 5

/** ****************************************************************************
 *
 *  PLASMA constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 **/
#define dplasmaRM            101
#define dplasmaCM            102
#define dplasmaCCRB          103
#define dplasmaCRRB          104
#define dplasmaRCRB          105
#define dplasmaRRRB          106

#define dplasmaNoTrans       111
#define dplasmaTrans         112
#define dplasmaConjTrans     113

#define dplasmaUpper         121
#define dplasmaLower         122
#define dplasmaUpperLower    123

#define dplasmaNonUnit       131
#define dplasmaUnit          132

#define dplasmaLeft          141
#define dplasmaRight         142

#define dplasmaOneNorm       171
#define dplasmaRealOneNorm   172
#define dplasmaTwoNorm       173
#define dplasmaFrobeniusNorm 174
#define dplasmaInfNorm       175
#define dplasmaRealInfNorm   176
#define dplasmaMaxNorm       177
#define dplasmaRealMaxNorm   178

#define dplasmaIncreasingOrder 181
#define dplasmaDecreasingOrder 182

#define dplasmaDistUniform   201
#define dplasmaDistSymmetric 202
#define dplasmaDistNormal    203

/* PLASMA_symmetry_type */
#define dplasmaGeneral                 231
#define dplasmaSymmetric               232
#define dplasmaHermitian               233
#define dplasmaTriangular              234
#define dplasmaLowerTriangular         235
#define dplasmaUpperTriangular         236
#define dplasmaLowerSymmetric          237
#define dplasmaUpperSymmetric          238
#define dplasmaLowerHermitian          239
#define dplasmaUpperHermitian          240
#define dplasmaSymetricBandLowerStored 241
#define dplasmaSymetricBandUpperStored 242
#define dplasmaBand                    243
#define dplasmaUpperHessenberg         244

/* Not found in CLapack - used by zlatms in libtmg */
#define dplasmaHermGeev      246
#define dplasmaHermPoev      247
#define dplasmaNonsymPosv    248
#define dplasmaSymPosv       249

#define dplasmaNoPacking     291
#define dplasmaPackSubdiag   292
#define dplasmaPackSupdiag   293
#define dplasmaPackColumn    294
#define dplasmaPackRow       295
#define dplasmaPackLowerBand 296
#define dplasmaPackUpeprBand 297
#define dplasmaPackAll       298

#define dplasmaNoVec         301
#define dplasmaVec           302
#define dplasmaIvec          303
#define dplasmaAllVec        304

#define dplasmaForward       391
#define dplasmaBackward      392

#define dplasmaColumnwise    401
#define dplasmaRowwise       402

#define dplasmaW             501
#define dplasmaA2            502

#define dplasmaLaed3Update1   0x01
#define dplasmaLaed3Update2   0x10
#define dplasmaLaed3UpdateAll 0x11


/**
 * Type of matrices that can be generated with PLASMA_zplrntx family
 * functions. (See PLASMA_zplrntx() for more details)
 * The list is coming from the LAWN 263.
 */
enum dplasma_matrix_type_e {
    dplasmaMatrixRandom    = 0,
    dplasmaMatrixHadamard  = 1,
    dplasmaMatrixHouse     = 2,
    dplasmaMatrixParter    = 3,
    dplasmaMatrixRis       = 4,
    dplasmaMatrixKms       = 5,
    dplasmaMatrixToeppen   = 6,   /* Unavailable */
    dplasmaMatrixCondex    = 7,
    dplasmaMatrixMoler     = 8,
    dplasmaMatrixCircul    = 9,
    dplasmaMatrixRandcorr  = 10,  /* Unavailable */
    dplasmaMatrixPoisson   = 11,  /* Unavailable */
    dplasmaMatrixHankel    = 12,
    dplasmaMatrixJordbloc  = 13,  /* Unavailable */
    dplasmaMatrixCompan    = 14,
    dplasmaMatrixPei       = 15,  /* Unavailable */
    dplasmaMatrixRandcolu  = 16,  /* Unavailable */
    dplasmaMatrixSprandn   = 17,  /* Unavailable */
    dplasmaMatrixRiemann   = 18,
    dplasmaMatrixCompar    = 19,  /* Unavailable */
    dplasmaMatrixTridiag   = 20,  /* Unavailable */
    dplasmaMatrixChebspec  = 21,  /* Unavailable */
    dplasmaMatrixLehmer    = 22,
    dplasmaMatrixToeppd    = 23,
    dplasmaMatrixMinij     = 24,
    dplasmaMatrixRandsvd   = 25,  /* Unavailable */
    dplasmaMatrixForsythe  = 26,  /* Unavailable */
    dplasmaMatrixFiedler   = 27,
    dplasmaMatrixDorr      = 28,
    dplasmaMatrixDemmel    = 29,
    dplasmaMatrixChebvand  = 30,
    dplasmaMatrixInvhess   = 31,
    dplasmaMatrixProlate   = 32,  /* Unavailable */
    dplasmaMatrixFrank     = 33,  /* Unavailable */
    dplasmaMatrixCauchy    = 34,
    dplasmaMatrixHilb      = 35,
    dplasmaMatrixLotkin    = 36,
    dplasmaMatrixKahan     = 37,  /* Unavailable */
    dplasmaMatrixOrthog    = 38,
    dplasmaMatrixWilkinson = 39,
    dplasmaMatrixFoster    = 40,
    dplasmaMatrixWright    = 41,
    dplasmaMatrixLangou    = 42
};

extern char *dplasma_lapack_const_strings[];
#define dplasma_lapack_const(plasma_const) (dplasma_lapack_const_strings[plasma_const][0])

#if defined(DPLASMA_HAVE_CUDA)
#include <cublas.h>

#define dplasma_cublas_side(side)                                         \
    assert( (side == dplasmaRight) || (side == dplasmaLeft) );            \
    side = (side == dplasmaRight) ? CUBLAS_SIDE_RIGHT : CUBLAS_SIDE_LEFT;


#define dplasma_cublas_diag(diag)                                              \
    assert( (diag == dplasmaNonUnit) || (diag == dplasmaUnit) );               \
    diag = (diag == dplasmaNonUnit) ? CUBLAS_DIAG_NON_UNIT : CUBLAS_DIAG_UNIT;

#define dplasma_cublas_fill(fill)                                                    \
    assert( (fill == dplasmaLower) || (fill == dplasmaUpper) );                      \
    fill = (fill == dplasmaLower) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;

#if defined(PRECISION_z) || defined(PRECISION_c)
#define dplasma_cublas_op(trans)                 \
    assert( (trans == dplasmaNoTrans) || (trans == dplasmaTrans) || (trans == dplasmaConjTrans) ); \
    switch(trans){                               \
        case dplasmaNoTrans:                     \
            trans = CUBLAS_OP_N;                 \
            break;                               \
        case dplasmaTrans:                       \
            trans = CUBLAS_OP_T;                 \
            break;                               \
        case dplasmaConjTrans:                   \
            trans = CUBLAS_OP_C;                 \
            break;                               \
        default:                                 \
            trans = CUBLAS_OP_N;                 \
            break;                               \
    }
#else
#define dplasma_cublas_op(trans)                                    \
    assert( (trans == dplasmaNoTrans) || (trans == dplasmaTrans) ); \
    trans = (trans == dplasmaNoTrans) ? CUBLAS_OP_N : CUBLAS_OP_T;
#endif /* PRECISION_z || PRECISION_c */

#endif /* DPLASMA_HAVE_CUDA */

#endif /* _DPLASMA_CONSTANTS_H_ */
