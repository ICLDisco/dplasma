/**
 * Copyright (c) 2020      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 **/
#ifndef _PLASMA_COMMON_H_
#define _PLASMA_COMMON_H_

/**
 * Just a placeholder to avoid changing imported plasma files. It contains some
 * of the macro definitions that are in control/xyz.h in PLASMA that we did not
 * import straight.
 */
#include "cores/dplasma_plasmatypes.h"
#include "cores/descriptor.h"
#include "core_blas.h"

#define ELTADDR(A, type, m, n)  (type *)plasma_geteltaddr(A, m, n)
#define ELTLDD(A, k) ( ( (((k)-1)/(A).mb) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )
#define BLKADDR(A, type, m, n)  (type *)plasma_getaddr(A, m, n)
#define BLKLDD(A, k) ( ( (k) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )

static inline int coreblas_imin(int a, int b) { return (a < b) ? a : b; };
#ifndef min
#define min(a, b) coreblas_imin(a, b)
#endif
static inline int coreblas_imax(int a, int b) { return (a > b) ? a : b; };
#ifndef max
#define max(a, b) coreblas_imax(a,b)
#endif

#endif /* _PLASMA_COMMON_H_ */
