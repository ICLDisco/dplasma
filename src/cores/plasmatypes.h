/**
 * Copyright (c) 2020      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * Type correspondance between PLASMA Cores and DPLASMA
 * The original PLASMA plasmatypes.h is split between the DPLASMA API
 * include/dplasma/constants.h and include/dplasma/complex.h
 *
 * When updating to a new PLASMA version, import the plasmatypes.h in these API 
 * files. Then apply the substitutions documented below to produce this file.
 *
 **/
#ifndef _PLASMATYPES_H_
#define _PLASMATYPES_H_

/** ****************************************************************************
 *  PLASMA Release number
 **/
#define PLASMA_VERSION_MAJOR 2
#define PLASMA_VERSION_MINOR 8
#define PLASMA_VERSION_MICRO 0


#include "dplasma/constants.h"

/* converted from plasmatypes.h with:'<,'>s/PLASMA_\([a-zA-Z_]\+\).*$/PLASMA_\1 DPLASMA_\1/ */
#define PLASMA_SUCCESS DPLASMA_SUCCESS
#define PLASMA_ERR_NOT_INITIALIZED DPLASMA_ERR_NOT_INITIALIZED
#define PLASMA_ERR_REINITIALIZED DPLASMA_ERR_REINITIALIZED
#define PLASMA_ERR_NOT_SUPPORTED DPLASMA_ERR_NOT_SUPPORTED
#define PLASMA_ERR_ILLEGAL_VALUE DPLASMA_ERR_ILLEGAL_VALUE
#define PLASMA_ERR_NOT_FOUND DPLASMA_ERR_NOT_FOUND
#define PLASMA_ERR_OUT_OF_RESOURCES DPLASMA_ERR_OUT_OF_RESOURCES
#define PLASMA_ERR_INTERNAL_LIMIT DPLASMA_ERR_INTERNAL_LIMIT
#define PLASMA_ERR_UNALLOCATED DPLASMA_ERR_UNALLOCATED
#define PLASMA_ERR_FILESYSTEM DPLASMA_ERR_FILESYSTEM
#define PLASMA_ERR_UNEXPECTED DPLASMA_ERR_UNEXPECTED
#define PLASMA_ERR_SEQUENCE_FLUSHED DPLASMA_ERR_SEQUENCE_FLUSHED

/** ****************************************************************************
 *  PLASMA types
 **/
typedef dplasma_enum_t  PLASMA_enum;
typedef dplasma_bool_t  PLASMA_bool;
typedef dplasma_index_t PLASMA_index;
typedef dplasma_size_t  PLASMA_size;

/** ****************************************************************************
 *  PLASMA constants - precisions
 **/
/* converted from plasmatypes.h with:'<,'>s/Plasma\([a-zA-Z0-9]\+\).*$/Plasma\1 dplasma\1/ */
#define PlasmaByte dplasmaByte
#define PlasmaInteger dplasmaInteger
#define PlasmaRealFloat dplasmaRealFloat
#define PlasmaRealDouble dplasmaRealDouble
#define PlasmaComplexFloat dplasmaComplexFloat
#define PlasmaComplexDouble dplasmaComplexDouble

/** ****************************************************************************
 *
 *  PLASMA constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 **/
#define PlasmaRM dplasmaRM
#define PlasmaCM dplasmaCM
#define PlasmaCCRB dplasmaCCRB
#define PlasmaCRRB dplasmaCRRB
#define PlasmaRCRB dplasmaRCRB
#define PlasmaRRRB dplasmaRRRB

#define PlasmaNoTrans dplasmaNoTrans
#define PlasmaTrans dplasmaTrans
#define PlasmaConjTrans dplasmaConjTrans

#define PlasmaUpper dplasmaUpper
#define PlasmaLower dplasmaLower
#define PlasmaUpperLower dplasmaUpperLower

#define PlasmaNonUnit dplasmaNonUnit
#define PlasmaUnit dplasmaUnit

#define PlasmaLeft dplasmaLeft
#define PlasmaRight dplasmaRight

#define PlasmaOneNorm dplasmaOneNorm
#define PlasmaRealOneNorm dplasmaRealOneNorm
#define PlasmaTwoNorm dplasmaTwoNorm
#define PlasmaFrobeniusNorm dplasmaFrobeniusNorm
#define PlasmaInfNorm dplasmaInfNorm
#define PlasmaRealInfNorm dplasmaRealInfNorm
#define PlasmaMaxNorm dplasmaMaxNorm
#define PlasmaRealMaxNorm dplasmaRealMaxNorm

#define PlasmaIncreasingOrder dplasmaIncreasingOrder
#define PlasmaDecreasingOrder dplasmaDecreasingOrder

#define PlasmaDistUniform dplasmaDistUniform
#define PlasmaDistSymmetric dplasmaDistSymmetric
#define PlasmaDistNormal dplasmaDistNormal

/* PLASMA_symmetry_type */
#define PlasmaGeneral dplasmaGeneral
#define PlasmaSymmetric dplasmaSymmetric
#define PlasmaHermitian dplasmaHermitian
#define PlasmaTriangular dplasmaTriangular
#define PlasmaLowerTriangular dplasmaLowerTriangular
#define PlasmaUpperTriangular dplasmaUpperTriangular
#define PlasmaLowerSymmetric dplasmaLowerSymmetric
#define PlasmaUpperSymmetric dplasmaUpperSymmetric
#define PlasmaLowerHermitian dplasmaLowerHermitian
#define PlasmaUpperHermitian dplasmaUpperHermitian
#define PlasmaSymetricBandLowerStored dplasmaSymetricBandLowerStored
#define PlasmaSymetricBandUpperStored dplasmaSymetricBandUpperStored
#define PlasmaBand dplasmaBand
#define PlasmaUpperHessenberg dplasmaUpperHessenberg

/* Not found in CLapack - used by zlatms in libtmg */
#define PlasmaHermGeev dplasmaHermGeev
#define PlasmaHermPoev dplasmaHermPoev
#define PlasmaNonsymPosv dplasmaNonsymPosv
#define PlasmaSymPosv dplasmaSymPosv

#define PlasmaNoPacking dplasmaNoPacking
#define PlasmaPackSubdiag dplasmaPackSubdiag
#define PlasmaPackSupdiag dplasmaPackSupdiag
#define PlasmaPackColumn dplasmaPackColumn
#define PlasmaPackRow dplasmaPackRow
#define PlasmaPackLowerBand dplasmaPackLowerBand
#define PlasmaPackUpeprBand dplasmaPackUpeprBand
#define PlasmaPackAll dplasmaPackAll

#define PlasmaNoVec dplasmaNoVec
#define PlasmaVec dplasmaVec
#define PlasmaIvec dplasmaIvec
#define PlasmaAllVec dplasmaAllVec

#define PlasmaForward dplasmaForward
#define PlasmaBackward dplasmaBackward

#define PlasmaColumnwise dplasmaColumnwise
#define PlasmaRowwise dplasmaRowwise

#define PlasmaW dplasmaW
#define PlasmaA2 dplasmaA2

#define PlasmaLaed3Update1 dplasmaLaed3Update1
#define PlasmaLaed3Update2 dplasmaLaed3Update2
#define PlasmaLaed3UpdateAll dplasmaLaed3UpdatAll


/**
 * Type of matrices that can be generated with PLASMA_zplrntx family
 * functions. (See PLASMA_zplrntx() for more details)
 * The list is coming from the LAWN 263.
 */
enum plasma_matrix_type_e {
    PlasmaMatrixRandom = dplasmaMatrixRandom,
    PlasmaMatrixHadamard = dplasmaMatrixHadamard,
    PlasmaMatrixHouse = dplasmaMatrixHouse,
    PlasmaMatrixParter = dplasmaMatrixParter,
    PlasmaMatrixRis = dplasmaMatrixRis,
    PlasmaMatrixKms = dplasmaMatrixKms,
    PlasmaMatrixToeppen = dplasmaMatrixToeppen,
    PlasmaMatrixCondex = dplasmaMatrixCondex,
    PlasmaMatrixMoler = dplasmaMatrixMoler,
    PlasmaMatrixCircul = dplasmaMatrixCircul,
    PlasmaMatrixRandcorr = dplasmaMatrixRandcorr,
    PlasmaMatrixPoisson = dplasmaMatrixPoisson,
    PlasmaMatrixHankel = dplasmaMatrixHankel,
    PlasmaMatrixJordbloc = dplasmaMatrixJordbloc,
    PlasmaMatrixCompan = dplasmaMatrixCompan,
    PlasmaMatrixPei = dplasmaMatrixPei,
    PlasmaMatrixRandcolu = dplasmaMatrixRandcolu,
    PlasmaMatrixSprandn = dplasmaMatrixSprandn,
    PlasmaMatrixRiemann = dplasmaMatrixRiemann,
    PlasmaMatrixCompar = dplasmaMatrixCompar,
    PlasmaMatrixTridiag = dplasmaMatrixTridiag,
    PlasmaMatrixChebspec = dplasmaMatrixChebspec,
    PlasmaMatrixLehmer = dplasmaMatrixLehmer,
    PlasmaMatrixToeppd = dplasmaMatrixToeppd,
    PlasmaMatrixMinij = dplasmaMatrixMinij,
    PlasmaMatrixRandsvd = dplasmaMatrixRandsvd,
    PlasmaMatrixForsythe = dplasmaMatrixForsythe,
    PlasmaMatrixFiedler = dplasmaMatrixFiedler,
    PlasmaMatrixDorr = dplasmaMatrixDorr,
    PlasmaMatrixDemmel = dplasmaMatrixDemmel,
    PlasmaMatrixChebvand = dplasmaMatrixChebvand,
    PlasmaMatrixInvhess = dplasmaMatrixInvhess,
    PlasmaMatrixProlate = dplasmaMatrixProlate,
    PlasmaMatrixFrank = dplasmaMatrixFrank,
    PlasmaMatrixCauchy = dplasmaMatrixCauchy,
    PlasmaMatrixHilb = dplasmaMatrixHilb,
    PlasmaMatrixLotkin = dplasmaMatrixLotkin,
    PlasmaMatrixKahan = dplasmaMatrixKahan,
    PlasmaMatrixOrthog = dplasmaMatrixOrthog,
    PlasmaMatrixWilkinson = dplasmaMatrixWilkinson,
    PlasmaMatrixFoster = dplasmaMatrixFoster,
    PlasmaMatrixWright = dplasmaMatrixWright,
    PlasmaMatrixLangou = dplasmaMatrixLangou,
};

/* rename the global symbol */
#define plasma_lapack_constants dplasma_lapack_constants
#define lapack_const(i) dplasma_lapack_const(i)

#include "dplasma/complex.h"

typedef dplasma_complex32_t PLASMA_Complex32_t;
typedef dplasma_complex64_t PLASMA_Complex64_t;

#endif /* _PLASMATYPES_H_ */
