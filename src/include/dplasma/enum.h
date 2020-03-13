/*
 * Copyright (c) 2020      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 */
#ifndef _DPLASMA_ENUM_H_
#define _DPLASMA_ENUM_H_




/* these numbers must fit the numbers declared in cores/plasmatypes.h as
 * imported from PLASMA
 */

typedef enum dplasma_enum_e {
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
    dplasmaMatrixLangou    = 42,

    dplasmaNoTrans = 111,
    dplasmaTrans = 112,
    dplasmaConjTrans = 113,

    dplasmaUpper = 121,
    dplasmaLower = 122,
    dplasmaUpperLower = 123,

    dplasmaNonUnit = 131,
    dplasmaUnit = 132,

    dplasmaLeft = 141,
    dplasmaRight = 142,

    dplasmaOneNorm = 171,
    dplasmaRealOneNorm = 172,
    dplasmaTwoNorm = 173,
    dplasmaFrobeniusNorm = 174,
    dplasmaInfNorm = 175,
    dplasmaRealInfNorm = 176,
    dplasmaMaxNorm = 177,
    dplasmaRealMaxNorm = 178,

    dplasmaGeneral = 231,
    dplasmaSymmetric = 232,
    dplasmaHermitian = 233,
    dplasmaTriangular = 234,
    dplasmaLowerTriangular = 235,
    dplasmaUpperTriangular = 236,
    dplasmaLowerSymetric = 237,
    dplasmaUpperSymetric = 238,
    dplasmaLowerHermitian = 239,
    dplasmaUpperHermitian = 240,
    dplasmaSymetricBandLowerStored = 241,
    dplasmaSymetricBandUpperStored = 242,
    dplasmaBand = 243,
    dplasmaUpperHessenberg = 244,


    dplasmaNoVec = 301,
    dplasmaVec = 302,
    dplasmaIvec = 303,
    dplasmaAllVec = 304,

    dplasmaForward = 391,
    dplasmaBackward = 392,

    dplasmaColumnwise = 401,
    dplasmaRowwise = 402
} dplasma_enum_t;

extern char *dplasma_lapack_constants[];
#define dplasma_lapack_const(plasma_const) dplasma_lapack_constants[plasma_const][0]

#endif /* _DPLASMA_ENUM_H_ */
