/**
 * Copyright (c) 2019      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Imported from:
 *
 * @file core_zblas.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#ifndef _PLASMA_CORE_ZBLAS_H_
#define _PLASMA_CORE_ZBLAS_H_

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

struct CORE_zgetrf_data_s;
typedef struct CORE_zgetrf_data_s CORE_zgetrf_data_t;

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/
int CORE_zamax(PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
               const PLASMA_Complex64_t *A, int lda, double *work);
int CORE_zamax_tile( PLASMA_enum storev, PLASMA_enum uplo, const PLASMA_desc descA, double *work);
void CORE_dzasum(int storev, PLASMA_enum uplo, int M, int N,
                 const PLASMA_Complex64_t *A, int lda, double *work);
void CORE_zbrdalg1( PLASMA_enum uplo,
                    int n,
                    int nb,
                    PLASMA_Complex64_t *A,
                    int lda,
                    PLASMA_Complex64_t *VQ,
                    PLASMA_Complex64_t *TAUQ,
                    PLASMA_Complex64_t *VP,
                    PLASMA_Complex64_t *TAUP,
                    int Vblksiz, int wantz,
                    int i, int sweepid, int m, int grsiz,
                    PLASMA_Complex64_t *work);
int CORE_zgbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zgbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zgblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zgeadd(PLASMA_enum trans, int M, int N,
                      PLASMA_Complex64_t alpha,
                const PLASMA_Complex64_t *A, int LDA,
                      PLASMA_Complex64_t beta,
                      PLASMA_Complex64_t *B, int LDB);
int  CORE_zgelqt(int M, int N, int IB,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *TAU,
                 PLASMA_Complex64_t *WORK);
void CORE_zgemm(PLASMA_enum transA, PLASMA_enum transB,
                int M, int N, int K,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                          const PLASMA_Complex64_t *B, int LDB,
                PLASMA_Complex64_t beta,        PLASMA_Complex64_t *C, int LDC);
void CORE_zgemv(PLASMA_enum trans, int M, int N,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                          const PLASMA_Complex64_t *x, int incx,
                PLASMA_Complex64_t beta,        PLASMA_Complex64_t *y, int incy);
void CORE_zgeqp3_init( int n, int *jpvt );
void CORE_zgeqp3_larfg( PLASMA_desc A, int ii, int jj, int i, int j,
                        PLASMA_Complex64_t *tau, PLASMA_Complex64_t *beta );
void CORE_zgeqp3_norms( PLASMA_desc A, int ioff, int joff, double *norms1, double *norms2 );
void CORE_zgeqp3_pivot( PLASMA_desc A, PLASMA_Complex64_t *F, int ldf,
                        int jj, int k, int *jpvt,
                        double *norms1, double *norms2, int *info );
int  CORE_zgeqp3_tntpiv(int m, int n,
                        PLASMA_Complex64_t *A, int lda,
                        int *IPIV, PLASMA_Complex64_t *tau,
                        int *iwork);
void CORE_zgeqp3_update( const PLASMA_Complex64_t *Ajj, int lda1,
                         PLASMA_Complex64_t       *Ajk, int lda2,
                         const PLASMA_Complex64_t *Fk,  int ldf,
                         int joff, int k, int koff, int nb,
                         double *norms1, double *norms2,
                         int *info );
int  CORE_zgeqrt(int M, int N, int IB,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK);
int  CORE_zgessm(int M, int N, int K, int IB,
                 const int *IPIV,
                 const PLASMA_Complex64_t *L, int LDL,
                 PLASMA_Complex64_t *A, int LDA);
int  CORE_zgessq(int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *scale, double *sumsq);
int  CORE_zgetf2_nopiv(int m, int n,
                      PLASMA_Complex64_t *A, int lda);
int  CORE_zgetrf(int M, int N,
                 PLASMA_Complex64_t *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_zgetrf_incpiv(int M, int N, int IB,
                        PLASMA_Complex64_t *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_zgetrf_nopiv(int m, int n, int ib,
                      PLASMA_Complex64_t *A, int lda);
int  CORE_zgetrf_reclap(CORE_zgetrf_data_t *data, int M, int N,
                        PLASMA_Complex64_t *A, int LDA,
                        int *IPIV, int *info);
CORE_zgetrf_data_t *CORE_zgetrf_reclap_init(int nbthrd);
int  CORE_zgetrf_rectil(CORE_zgetrf_data_t *data, const PLASMA_desc A, int *IPIV, int *info);
CORE_zgetrf_data_t *CORE_zgetrf_rectil_init(int nbthrd);
void CORE_zgetrip(int m, int n, PLASMA_Complex64_t *A,
                  PLASMA_Complex64_t *work);
int CORE_zhbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zhblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zhbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                int st, int ed, int eltsize);
void CORE_zhbtype1cb(int N, int NB,
                     PLASMA_Complex64_t *A, int LDA,
                     PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     PLASMA_Complex64_t *WORK);
void CORE_zhbtype2cb(int N, int NB,
                     PLASMA_Complex64_t *A, int LDA,
                     PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     PLASMA_Complex64_t *WORK);
void CORE_zhbtype3cb(int N, int NB,
                     PLASMA_Complex64_t *A, int LDA,
                     const PLASMA_Complex64_t *V, const PLASMA_Complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     PLASMA_Complex64_t *WORK);
void CORE_zgbtype1cb(PLASMA_enum uplo, int N, int NB,
                PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t *VQ, PLASMA_Complex64_t *TAUQ,
                PLASMA_Complex64_t *VP, PLASMA_Complex64_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                PLASMA_Complex64_t *WORK);
void CORE_zgbtype2cb(PLASMA_enum uplo, int N, int NB,
                PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t *VQ, PLASMA_Complex64_t *TAUQ,
                PLASMA_Complex64_t *VP, PLASMA_Complex64_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                PLASMA_Complex64_t *WORK);
void CORE_zgbtype3cb(PLASMA_enum uplo, int N, int NB,
                PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t *VQ, PLASMA_Complex64_t *TAUQ,
                PLASMA_Complex64_t *VP, PLASMA_Complex64_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                PLASMA_Complex64_t *WORK);
void CORE_zhegst(int itype, PLASMA_enum uplo, int N,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *B, int LDB, int *INFO);
#ifdef COMPLEX
void CORE_zhemm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                          const PLASMA_Complex64_t *B, int LDB,
                PLASMA_Complex64_t beta,        PLASMA_Complex64_t *C, int LDC);
void CORE_zherk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                double alpha, const PLASMA_Complex64_t *A, int LDA,
                double beta,        PLASMA_Complex64_t *C, int LDC);
void CORE_zher2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                           const PLASMA_Complex64_t *B, int LDB,
                 double beta,                    PLASMA_Complex64_t *C, int LDC);
int  CORE_zhessq(PLASMA_enum uplo, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *scale, double *sumsq);
#endif
int  CORE_zherfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 const PLASMA_Complex64_t *A,    int LDA,
                 const PLASMA_Complex64_t *T,    int LDT,
                       PLASMA_Complex64_t *C,    int LDC,
                       PLASMA_Complex64_t *WORK, int LDWORK);
void CORE_zlacpy(PLASMA_enum uplo, int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                       PLASMA_Complex64_t *B, int LDB);
int CORE_zlacpy_pivot( const PLASMA_desc descA,
                       PLASMA_enum direct,
                       int k1, int k2, const int *ipiv,
                       int *rankin, int *rankout,
                       PLASMA_Complex64_t *A, int lda,
                       int init);
void CORE_zlange(int norm, int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *work, double *normA);
#ifdef COMPLEX
void CORE_zlanhe(int norm, PLASMA_enum uplo, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *work, double *normA);
#endif
void CORE_zlansy(int norm, PLASMA_enum uplo, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *work, double *normA);
void CORE_zlantr(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *work, double *normA);
int CORE_zlarfb_gemm(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                     int M, int N, int K,
                     const PLASMA_Complex64_t *V, int LDV,
                     const PLASMA_Complex64_t *T, int LDT,
                           PLASMA_Complex64_t *C, int LDC,
                           PLASMA_Complex64_t *WORK, int LDWORK);
int CORE_zlarfx2(PLASMA_enum side, int N,
                 PLASMA_Complex64_t V,
                 PLASMA_Complex64_t TAU,
                 PLASMA_Complex64_t *C1, int LDC1,
                 PLASMA_Complex64_t *C2, int LDC2);
int CORE_zlarfx2c(PLASMA_enum uplo,
                  PLASMA_Complex64_t V,
                  PLASMA_Complex64_t TAU,
                  PLASMA_Complex64_t *C1,
                  PLASMA_Complex64_t *C2,
                  PLASMA_Complex64_t *C3);
int CORE_zlarfx2ce(PLASMA_enum uplo,
                   PLASMA_Complex64_t *V,
                   PLASMA_Complex64_t *TAU,
                   PLASMA_Complex64_t *C1,
                   PLASMA_Complex64_t *C2,
                   PLASMA_Complex64_t *C3);
void CORE_zlarfy(int N,
                 PLASMA_Complex64_t *A, int LDA,
                 const PLASMA_Complex64_t *V,
                 const PLASMA_Complex64_t *TAU,
                 PLASMA_Complex64_t *WORK);
int  CORE_zlascal(PLASMA_enum uplo, int m, int n,
                  PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int lda);
void CORE_zlaset(PLASMA_enum uplo, int n1, int n2,
                 PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta,
                 PLASMA_Complex64_t *tileA, int ldtilea);
void CORE_zlaset2(PLASMA_enum uplo, int n1, int n2, PLASMA_Complex64_t alpha,
                  PLASMA_Complex64_t *tileA, int ldtilea);
void CORE_zlaswp(int N, PLASMA_Complex64_t *A, int LDA,
                 int I1,  int I2, const int *IPIV, int INC);
int  CORE_zlaswp_ontile( PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_zlaswpc_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_zlatro(PLASMA_enum uplo, PLASMA_enum trans,
                 int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                       PLASMA_Complex64_t *B, int LDB);
void CORE_zlauum(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA);
int CORE_zpamm(int op, PLASMA_enum side, PLASMA_enum storev,
               int M, int N, int K, int L,
               const PLASMA_Complex64_t *A1, int LDA1,
                     PLASMA_Complex64_t *A2, int LDA2,
               const PLASMA_Complex64_t *V, int LDV,
                     PLASMA_Complex64_t *W, int LDW);
int  CORE_zparfb(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                 int M1, int N1, int M2, int N2, int K, int L,
                       PLASMA_Complex64_t *A1, int LDA1,
                       PLASMA_Complex64_t *A2, int LDA2,
                 const PLASMA_Complex64_t *V, int LDV,
                 const PLASMA_Complex64_t *T, int LDT,
                       PLASMA_Complex64_t *WORK, int LDWORK);
int CORE_zpemv(PLASMA_enum trans, PLASMA_enum storev,
               int M, int N, int L,
               PLASMA_Complex64_t ALPHA,
               const PLASMA_Complex64_t *A, int LDA,
               const PLASMA_Complex64_t *X, int INCX,
               PLASMA_Complex64_t BETA,
               PLASMA_Complex64_t *Y, int INCY,
               PLASMA_Complex64_t *WORK);
void CORE_zplghe(double bump, int m, int n, PLASMA_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_zplgsy(PLASMA_Complex64_t bump, int m, int n, PLASMA_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_zplrnt(int m, int n, PLASMA_Complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
int  CORE_zpltmg(PLASMA_enum mtxtype, int m, int n, PLASMA_Complex64_t *A, int lda,
                  int gM, int gN, int m0, int n0, unsigned long long int seed );
int  CORE_zpltmg_chebvand( int M, int N, PLASMA_Complex64_t *A, int LDA,
                           int gN, int m0, int n0,
                           PLASMA_Complex64_t *W );
int  CORE_zpltmg_circul( int M, int N, PLASMA_Complex64_t *A, int LDA,
                         int gM, int m0, int n0,
                         const PLASMA_Complex64_t *V );
void CORE_zpltmg_condexq( int M, int N, PLASMA_Complex64_t *Q, int LDQ );
void CORE_zpltmg_fiedler(int m, int n,
                         const PLASMA_Complex64_t *X, int incX,
                         const PLASMA_Complex64_t *Y, int incY,
                         PLASMA_Complex64_t *A, int lda);
int  CORE_zpltmg_hankel( PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t *A, int LDA,
                         int m0, int n0, int nb,
                         const PLASMA_Complex64_t *V1,
                         const PLASMA_Complex64_t *V2 );
void CORE_zpltmg_toeppd1( int gM, int m0, int M, PLASMA_Complex64_t *W,
                          unsigned long long int seed );
void CORE_zpltmg_toeppd2( int M, int N, int K, int m0, int n0,
                          const PLASMA_Complex64_t *W,
                          PLASMA_Complex64_t *A, int LDA );
void CORE_zpotrf(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, int *INFO);
void CORE_zsetvar(const PLASMA_Complex64_t *alpha, PLASMA_Complex64_t *x);
void CORE_zshift(int s, int m, int n, int L,
                 PLASMA_Complex64_t *A);
void CORE_zshiftw(int s, int cl, int m, int n, int L,
                  PLASMA_Complex64_t *A, PLASMA_Complex64_t *W);
int  CORE_zssssm(int M1, int N1, int M2, int N2, int K, int IB,
                       PLASMA_Complex64_t *A1, int LDA1,
                       PLASMA_Complex64_t *A2, int LDA2,
                 const PLASMA_Complex64_t *L1, int LDL1,
                 const PLASMA_Complex64_t *L2, int LDL2,
                 const int *IPIV);
int CORE_zstedc(PLASMA_enum compz, int n,
                double *D, double *E,
                PLASMA_Complex64_t *Z, int LDZ,
                PLASMA_Complex64_t *WORK, int LWORK,
#ifdef COMPLEX
                double *RWORK, int LRWORK,
#endif
                int *IWORK, int LIWORK);
int CORE_zsteqr(PLASMA_enum compz, int n,
                double *D, double *E,
                PLASMA_Complex64_t *Z, int LDZ,
                double *WORK);
void CORE_zsymm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                          const PLASMA_Complex64_t *B, int LDB,
                PLASMA_Complex64_t beta,        PLASMA_Complex64_t *C, int LDC);
void CORE_zsyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t beta,        PLASMA_Complex64_t *C, int LDC);
void CORE_zsyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                           const PLASMA_Complex64_t *B, int LDB,
                 PLASMA_Complex64_t beta,        PLASMA_Complex64_t *C, int LDC);
int  CORE_zsyssq(PLASMA_enum uplo, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *scale, double *sumsq);
void CORE_zswpab(int i, int n1, int n2,
                 PLASMA_Complex64_t *A, PLASMA_Complex64_t *work);
int  CORE_zswptr_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc,
                        const PLASMA_Complex64_t *Akk, int ldak);
int CORE_ztradd(PLASMA_enum uplo, PLASMA_enum trans, int M, int N,
                      PLASMA_Complex64_t alpha,
                const PLASMA_Complex64_t *A, int LDA,
                      PLASMA_Complex64_t beta,
                      PLASMA_Complex64_t *B, int LDB);
void CORE_ztrasm(PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N, const PLASMA_Complex64_t *A, int lda, double *work);
void CORE_ztrdalg1(int n,
                        int nb,
                        PLASMA_Complex64_t *A,
                        int lda,
                        PLASMA_Complex64_t *V,
                        PLASMA_Complex64_t *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        PLASMA_Complex64_t *work);
void CORE_ztrmm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                                PLASMA_Complex64_t *B, int LDB);
void CORE_ztrsm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                                PLASMA_Complex64_t *B, int LDB);
int  CORE_ztrssq(PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *scale, double *sumsq);
void CORE_ztrtri(PLASMA_enum uplo, PLASMA_enum diag, int N,
                 PLASMA_Complex64_t *A, int LDA, int *info);
int  CORE_ztslqt(int M, int N, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK);
int  CORE_ztsmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 const PLASMA_Complex64_t *V, int LDV,
                 const PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *WORK, int LDWORK);
int CORE_ztsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        PLASMA_Complex64_t *A1, int lda1,
                        PLASMA_Complex64_t *A2, int lda2,
                        PLASMA_Complex64_t *A3, int lda3,
                        const PLASMA_Complex64_t *V, int ldv,
                        const PLASMA_Complex64_t *T, int ldt,
                        PLASMA_Complex64_t *WORK, int ldwork);
int CORE_ztsmlq_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        PLASMA_Complex64_t *A1, int lda1,
                        PLASMA_Complex64_t *A2, int lda2,
                        const PLASMA_Complex64_t *V, int ldv,
                        const PLASMA_Complex64_t *T, int ldt,
                        PLASMA_Complex64_t *WORK, int ldwork);
int  CORE_ztsmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 const PLASMA_Complex64_t *V, int LDV,
                 const PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *WORK, int LDWORK);
int CORE_ztsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        PLASMA_Complex64_t *A1, int lda1,
                        PLASMA_Complex64_t *A2, int lda2,
                        PLASMA_Complex64_t *A3, int lda3,
                        const PLASMA_Complex64_t *V, int ldv,
                        const PLASMA_Complex64_t *T, int ldt,
                        PLASMA_Complex64_t *WORK, int ldwork);
int CORE_ztsmqr_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        PLASMA_Complex64_t *A1, int lda1,
                        PLASMA_Complex64_t *A2, int lda2,
                        const PLASMA_Complex64_t *V, int ldv,
                        const PLASMA_Complex64_t *T, int ldt,
                        PLASMA_Complex64_t *WORK, int ldwork);
int  CORE_ztsqrt(int M, int N, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *TAU, PLASMA_Complex64_t *WORK);
int  CORE_ztstrf(int M, int N, int IB, int NB,
                 PLASMA_Complex64_t *U, int LDU,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *L, int LDL,
                 int *IPIV, PLASMA_Complex64_t *WORK,
                 int LDWORK, int *INFO);
int  CORE_zttmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 const PLASMA_Complex64_t *V, int LDV,
                 const PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *WORK, int LDWORK);
int  CORE_zttqrt(int M, int N, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *TAU,
                 PLASMA_Complex64_t *WORK);
int  CORE_zttmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 const PLASMA_Complex64_t *V, int LDV,
                 const PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *WORK, int LDWORK);
int  CORE_zttlqt(int M, int N, int IB,
                 PLASMA_Complex64_t *A1, int LDA1,
                 PLASMA_Complex64_t *A2, int LDA2,
                 PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *TAU,
                 PLASMA_Complex64_t *WORK);
int  CORE_zunmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const PLASMA_Complex64_t *V, int LDV,
                 const PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *C, int LDC,
                 PLASMA_Complex64_t *WORK, int LDWORK);
int  CORE_zunmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const PLASMA_Complex64_t *V, int LDV,
                 const PLASMA_Complex64_t *T, int LDT,
                 PLASMA_Complex64_t *C, int LDC,
                 PLASMA_Complex64_t *WORK, int LDWORK);

#ifndef COMPLEX
void CORE_dlaed2_computeK(int *K, int n, int n1,
                          double *beta, double *D, double *Q, int LDQ,
                          double *Z, double *DLAMBDA, double *W,
                          int *INDX, int *INDXC, int *INDXP, int *INDXQ,
                          int *COLTYP);
void CORE_dlaed2_compressq(int n, int n1, const int *INDX, const int *ctot,
                           const double *Q, int LDQ, double *Q2,
                           int start, int end);
void CORE_dlaed2_copydef(int n, int n1, int K, const int *ctot,
                         double *Q, int LDQ, const double *Q2,
                         int start, int end);
int CORE_dlaed4(int n, int K,
                double *D, double beta,
                double *Q, int LDQ,
                const double *D0, const double *Z,
                const int *INDX,
                int start, int end );
void CORE_dlaed3_computeW(int n, int K,
                          const double *Q, int LDQ,
                          const double *DLAMBDA, double *W,
                          const int *INDX,
                          int start, int end);
void CORE_dlaed3_reduceW(int n, int n1, int K, int l,
                         const double *Q, int LDQ,
                         const double *Wred, double *W);
void CORE_dlaed3_computevectors(int K, int il_nondef, int iu_nondef,
                                double *Q, int LDQ, double *W, double *S,
                                const int *INDXC,
                                int start, int end);
void CORE_dlaed3_merge( int n, int K, double *D, int *INDXQ );
void CORE_dlaed3_updatevectors(int op, int wsmode, int n, int n1, int K,
                   int il_nondef, int iu_nondef,
                               double *Q, int ldq, double *Q2,
                               const int *ctot, double *WORK, int start, int end);
#endif
void CORE_zswap(int m, int n, PLASMA_Complex64_t *Q, int ldq,
                const PLASMA_Complex64_t *work, const int *perm,
                int start, int end);
int CORE_zlascl(PLASMA_enum type, int kl, int ku, double cfrom, double cto,
                int m, int n, PLASMA_Complex64_t *A, int lda);
#ifdef COMPLEX
int CORE_dlag2z(int m, int n, const double *Q, int LDQ,
                 PLASMA_Complex64_t *Z, int LDZ);
#endif

#ifndef COMPLEX
void CORE_dlaed3_freebigwork(int oper, double **WORK);
void CORE_dlaed0_betaapprox(int subpbs, const int *subpbs_info,
                            double *D, const double *E);
int CORE_dlapst(PLASMA_enum type, int n,
                const double *D, int *INDX);
#endif

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
