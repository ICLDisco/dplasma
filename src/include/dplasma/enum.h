/*
 * Copyright (c) 2020      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 */
#ifndef _DPLASMA_ENUM_H_
#define _DPLASMA_ENUM_H_

typedef enum dplasma_enum_e {
    NoTrans = 111,
    Trans = 112,
    ConjTrans = 113,

    Upper = 121,
    Lower = 122,
    UpperLower = 123,

    NonUnit = 131,
    Unit = 132,

    Left = 141,
    Right = 142,

    OneNorm = 171,
    RealOneNorm = 172,
    TwoNorm = 173,
    FrobeniusNorm = 174,
    InfNorm = 175,
    RealInfNorm = 176,
    MaxNorm = 177,
    RealMaxNorm = 178,

    General = 231,
    Symmetric = 232,
    Hermitian = 233,
    Triangular = 234,
    LowerTriangular = 235,
    UpperTriangular = 236,
    LowerSymetric = 237,
    UpperSymetric = 238,
    LowerHermitian = 239,
    UpperHermitian = 240,
    SymetricBandLowerStored = 241,
    SymetricBandUpperStored = 242,
    Band = 243,
    UpperHessenberg = 244
} dplasma_enum_t;


#endif /* _DPLASMA_ENUM_H_ */
