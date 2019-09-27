/**
 *
 * @file gkkleader.h
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/

#ifndef GKKLEADERS_H
#define GKKLEADERS_H

int  GKK_doublingtable(int x, int m, int emax, int *dt);
int  GKK_modpow(int *dt, int e, int m);
int  GKK_primroot(int p, int e, primedec_t *pr_p1, int t_p1);
int  GKK_multorder(int n, int p, int e, int pe, primedec_t *pr_p1, int t_p1);
void GKK_prepare(int q, int n, primedec_t *pr_q, int *t,
                 primedec_t **pr_pi1, int *t_pi1, int *gi,
                 int *Nif, int *Kif, int *dif);
void GKK_L(int t, primedec_t *pr_q, int *fi, int *Kif, int *dif,
           int *Li, int *diLi, int *cl, int *nl);
void GKK_precompute_terms(int q, primedec_t *pr_q, int t, int *gi,
                          int *diLi, int *rp, int *Mg, int nMg);
void GKK_output_tables(int m, int n, int q, primedec_t *pr_q, int t,
                       int *gi, int *Nif, int *Kif, int *dif);

#endif /* GKKLEADERS_H */
