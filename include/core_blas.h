/**
 *
 * @file core_blas.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_CORE_BLAS_H_
#define _PLASMA_CORE_BLAS_H_

#include <cblas.h>
#include <lapacke.h>
#include "plasmatypes.h"
#include "cores/descriptor.h"

#include "cores/core_zblas.h"
#include "cores/core_dblas.h"
#include "cores/core_cblas.h"
#include "cores/core_sblas.h"
#include "cores/core_zcblas.h"
#include "cores/core_dsblas.h"

#ifdef __cplusplus
extern "C" {
#endif

  static inline int coreblas_imin(int a, int b) { return (a < b) ? a : b; };
  static inline int coreblas_imax(int a, int b) { return (a > b) ? a : b; };
  /*
   * Coreblas Error
   */
#define coreblas_error(k, str) fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str);

 /** ****************************************************************************
  *  LAPACK Constants
  **/
extern char *plasma_lapack_constants[];
#define lapack_const(plasma_const) plasma_lapack_constants[plasma_const][0]

/* CBLAS requires for scalar arguments to be passed by address rather than by value */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif

 /** ****************************************************************************
  *  External interface of the GKK algorithm for InPlace Layout Translation
  **/
int  GKK_minloc(int n, int *T);
void GKK_BalanceLoad(int thrdnbr, int *Tp, int *leaders, int nleaders, int L);
int  GKK_getLeaderNbr(int me, int ne, int *nleaders, int **leaders);

void CORE_pivot_update(int m, int n, int *ipiv, int *indices,
                       int offset, int init);

#ifdef __cplusplus
}
#endif

#endif /* _PLASMA_CORE_BLAS_H_ */
