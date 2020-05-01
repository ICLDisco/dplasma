/*
 * Copyright (c) 2010-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * expanded from:
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
 */

#ifndef _DPLASMA_COMPLEX_H_
#define _DPLASMA_COMPLEX_H_

/******************************************************************************
 * PaRSEC Complex numbers
 **/

#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
/* Windows and non-Intel compiler */
#include <float.h>
#include <complex>
typedef std::complex<float>  dplasma_complex32_t;
typedef std::complex<double> dplasma_complex64_t;
/* For LAPACKE lapacke.h force usage of Windows C++ Complex types */
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#define isnan _isnan
#define isinf !_finite
#else
typedef float  _Complex dplasma_complex32_t;
typedef double _Complex dplasma_complex64_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "dplasma/config.h"
#if !defined(__cplusplus) && defined(DPLASMA_HAVE_COMPLEX_H)
#include <complex.h>
#else

/* These declarations will not clash with what C++ provides because the names in C++ are name-mangled. */

extern double cabs     (dplasma_complex64_t z);
extern double carg     (dplasma_complex64_t z);
extern double creal    (dplasma_complex64_t z);
extern double cimag    (dplasma_complex64_t z);

extern float  cabsf    (dplasma_complex32_t z);
extern float  cargf    (dplasma_complex32_t z);
extern float  crealf   (dplasma_complex32_t z);
extern float  cimagf   (dplasma_complex32_t z);

extern dplasma_complex64_t conj  (dplasma_complex64_t z);
extern dplasma_complex64_t cproj (dplasma_complex64_t z);
extern dplasma_complex64_t csqrt (dplasma_complex64_t z);
extern dplasma_complex64_t cexp  (dplasma_complex64_t z);
extern dplasma_complex64_t clog  (dplasma_complex64_t z);
extern dplasma_complex64_t cpow  (dplasma_complex64_t z, dplasma_complex64_t w);

extern dplasma_complex32_t conjf (dplasma_complex32_t z);
extern dplasma_complex32_t cprojf(dplasma_complex32_t z);
extern dplasma_complex32_t csqrtf(dplasma_complex32_t z);
extern dplasma_complex32_t cexpf (dplasma_complex32_t z);
extern dplasma_complex32_t clogf (dplasma_complex32_t z);
extern dplasma_complex32_t cpowf (dplasma_complex32_t z, dplasma_complex32_t w);

#endif /* DPLASMA_HAVE_COMPLEX_H */

#ifdef __cplusplus
}
#endif

#endif /* _DPLASMA_COMPLEX_H_ */
