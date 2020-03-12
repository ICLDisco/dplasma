/**
 * Copyright (c) 2019      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Imported from:
 *
 * @file core_cblas.h
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
 * @generated c Fri Apr  1 11:02:20 2016
 *
 **/
#ifndef _PLASMA_CORE_CBLAS_H_
#define _PLASMA_CORE_CBLAS_H_

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

struct CORE_cgetrf_data_s;
typedef struct CORE_cgetrf_data_s CORE_cgetrf_data_t;

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/
void CORE_scasum(int storev, PLASMA_enum uplo, int M, int N,
                 const PLASMA_Complex32_t *A, int lda, float *work);
void CORE_cbrdalg1( PLASMA_enum uplo,
                    int n,
                    int nb,
                    PLASMA_Complex32_t *A,
                    int lda,
                    PLASMA_Complex32_t *VQ,
                    PLASMA_Complex32_t *TAUQ,
                    PLASMA_Complex32_t *VP,
                    PLASMA_Complex32_t *TAUP,
                    int Vblksiz, int wantz,
                    int i, int sweepid, int m, int grsiz,
                    PLASMA_Complex32_t *work);
int CORE_cgbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_cgbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_cgblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_cgeadd(PLASMA_enum trans, int M, int N,
                      PLASMA_Complex32_t alpha,
                const PLASMA_Complex32_t *A, int LDA,
                      PLASMA_Complex32_t beta,
                      PLASMA_Complex32_t *B, int LDB);
int  CORE_cgelqt(int M, int N, int IB,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU,
                 PLASMA_Complex32_t *WORK);
void CORE_cgemm(PLASMA_enum transA, PLASMA_enum transB,
                int M, int N, int K,
                PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                                          const PLASMA_Complex32_t *B, int LDB,
                PLASMA_Complex32_t beta,        PLASMA_Complex32_t *C, int LDC);
void CORE_cgemv(PLASMA_enum trans, int M, int N,
                PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                                          const PLASMA_Complex32_t *x, int incx,
                PLASMA_Complex32_t beta,        PLASMA_Complex32_t *y, int incy);
void CORE_cgeqp3_init( int n, int *jpvt );
void CORE_cgeqp3_larfg( PLASMA_desc A, int ii, int jj, int i, int j,
                        PLASMA_Complex32_t *tau, PLASMA_Complex32_t *beta );
void CORE_cgeqp3_norms( PLASMA_desc A, int ioff, int joff, float *norms1, float *norms2 );
void CORE_cgeqp3_pivot( PLASMA_desc A, PLASMA_Complex32_t *F, int ldf,
                        int jj, int k, int *jpvt,
                        float *norms1, float *norms2, int *info );
int  CORE_cgeqp3_tntpiv(int m, int n,
                        PLASMA_Complex32_t *A, int lda,
                        int *IPIV, PLASMA_Complex32_t *tau,
                        int *iwork);
void CORE_cgeqp3_update( const PLASMA_Complex32_t *Ajj, int lda1,
                         PLASMA_Complex32_t       *Ajk, int lda2,
                         const PLASMA_Complex32_t *Fk,  int ldf,
                         int joff, int k, int koff, int nb,
                         float *norms1, float *norms2,
                         int *info );
int  CORE_cgeqrt(int M, int N, int IB,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU, PLASMA_Complex32_t *WORK);
int  CORE_cgessm(int M, int N, int K, int IB,
                 const int *IPIV,
                 const PLASMA_Complex32_t *L, int LDL,
                 PLASMA_Complex32_t *A, int LDA);
int  CORE_cgessq(int M, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                 float *scale, float *sumsq);
int  CORE_cgetf2_nopiv(int m, int n,
                      PLASMA_Complex32_t *A, int lda);
int  CORE_cgetrf(int M, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_cgetrf_incpiv(int M, int N, int IB,
                        PLASMA_Complex32_t *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_cgetrf_nopiv(int m, int n, int ib,
                      PLASMA_Complex32_t *A, int lda);
int  CORE_cgetrf_reclap(CORE_cgetrf_data_t *data, int M, int N,
                        PLASMA_Complex32_t *A, int LDA,
                        int *IPIV, int *info);
CORE_cgetrf_data_t *CORE_cgetrf_reclap_init(int nbthrd);
int  CORE_cgetrf_rectil(CORE_cgetrf_data_t *data, const PLASMA_desc A, int *IPIV, int *info);
CORE_cgetrf_data_t *CORE_cgetrf_rectil_init(int nbthrd);
void CORE_cgetrip(int m, int n, PLASMA_Complex32_t *A,
                  PLASMA_Complex32_t *work);
int CORE_chbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_chblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_chbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU,
                int st, int ed, int eltsize);
void CORE_chbtype1cb(int N, int NB,
                     PLASMA_Complex32_t *A, int LDA,
                     PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     PLASMA_Complex32_t *WORK);
void CORE_chbtype2cb(int N, int NB,
                     PLASMA_Complex32_t *A, int LDA,
                     PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     PLASMA_Complex32_t *WORK);
void CORE_chbtype3cb(int N, int NB,
                     PLASMA_Complex32_t *A, int LDA,
                     const PLASMA_Complex32_t *V, const PLASMA_Complex32_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     PLASMA_Complex32_t *WORK);
void CORE_cgbtype1cb(PLASMA_enum uplo, int N, int NB,
                PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *VQ, PLASMA_Complex32_t *TAUQ,
                PLASMA_Complex32_t *VP, PLASMA_Complex32_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                PLASMA_Complex32_t *WORK);
void CORE_cgbtype2cb(PLASMA_enum uplo, int N, int NB,
                PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *VQ, PLASMA_Complex32_t *TAUQ,
                PLASMA_Complex32_t *VP, PLASMA_Complex32_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                PLASMA_Complex32_t *WORK);
void CORE_cgbtype3cb(PLASMA_enum uplo, int N, int NB,
                PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *VQ, PLASMA_Complex32_t *TAUQ,
                PLASMA_Complex32_t *VP, PLASMA_Complex32_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                PLASMA_Complex32_t *WORK);
void CORE_chegst(int itype, PLASMA_enum uplo, int N,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *B, int LDB, int *INFO);
#ifdef COMPLEX
void CORE_chemm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                                          const PLASMA_Complex32_t *B, int LDB,
                PLASMA_Complex32_t beta,        PLASMA_Complex32_t *C, int LDC);
void CORE_cherk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                float alpha, const PLASMA_Complex32_t *A, int LDA,
                float beta,        PLASMA_Complex32_t *C, int LDC);
void CORE_cher2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                                           const PLASMA_Complex32_t *B, int LDB,
                 float beta,                    PLASMA_Complex32_t *C, int LDC);
int  CORE_chessq(PLASMA_enum uplo, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                 float *scale, float *sumsq);
#endif
int  CORE_cherfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 const PLASMA_Complex32_t *A,    int LDA,
                 const PLASMA_Complex32_t *T,    int LDT,
                       PLASMA_Complex32_t *C,    int LDC,
                       PLASMA_Complex32_t *WORK, int LDWORK);
void CORE_clacpy(PLASMA_enum uplo, int M, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                       PLASMA_Complex32_t *B, int LDB);
int CORE_clacpy_pivot( const PLASMA_desc descA,
                       PLASMA_enum direct,
                       int k1, int k2, const int *ipiv,
                       int *rankin, int *rankout,
                       PLASMA_Complex32_t *A, int lda,
                       int init);
void CORE_clange(int norm, int M, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                 float *work, float *normA);
#ifdef COMPLEX
void CORE_clanhe(int norm, PLASMA_enum uplo, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                 float *work, float *normA);
#endif
void CORE_clansy(int norm, PLASMA_enum uplo, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                 float *work, float *normA);
void CORE_clantr(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                 float *work, float *normA);
int CORE_clarfb_gemm(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                     int M, int N, int K,
                     const PLASMA_Complex32_t *V, int LDV,
                     const PLASMA_Complex32_t *T, int LDT,
                           PLASMA_Complex32_t *C, int LDC,
                           PLASMA_Complex32_t *WORK, int LDWORK);
int CORE_clarfx2(PLASMA_enum side, int N,
                 PLASMA_Complex32_t V,
                 PLASMA_Complex32_t TAU,
                 PLASMA_Complex32_t *C1, int LDC1,
                 PLASMA_Complex32_t *C2, int LDC2);
int CORE_clarfx2c(PLASMA_enum uplo,
                  PLASMA_Complex32_t V,
                  PLASMA_Complex32_t TAU,
                  PLASMA_Complex32_t *C1,
                  PLASMA_Complex32_t *C2,
                  PLASMA_Complex32_t *C3);
int CORE_clarfx2ce(PLASMA_enum uplo,
                   PLASMA_Complex32_t *V,
                   PLASMA_Complex32_t *TAU,
                   PLASMA_Complex32_t *C1,
                   PLASMA_Complex32_t *C2,
                   PLASMA_Complex32_t *C3);
void CORE_clarfy(int N,
                 PLASMA_Complex32_t *A, int LDA,
                 const PLASMA_Complex32_t *V,
                 const PLASMA_Complex32_t *TAU,
                 PLASMA_Complex32_t *WORK);
int  CORE_clascal(PLASMA_enum uplo, int m, int n,
                  PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda);
void CORE_claset(PLASMA_enum uplo, int n1, int n2,
                 PLASMA_Complex32_t alpha, PLASMA_Complex32_t beta,
                 PLASMA_Complex32_t *tileA, int ldtilea);
void CORE_claset2(PLASMA_enum uplo, int n1, int n2, PLASMA_Complex32_t alpha,
                  PLASMA_Complex32_t *tileA, int ldtilea);
void CORE_claswp(int N, PLASMA_Complex32_t *A, int LDA,
                 int I1,  int I2, const int *IPIV, int INC);
int  CORE_claswp_ontile( PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_claswpc_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_clatro(PLASMA_enum uplo, PLASMA_enum trans,
                 int M, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                       PLASMA_Complex32_t *B, int LDB);
void CORE_clauum(PLASMA_enum uplo, int N, PLASMA_Complex32_t *A, int LDA);
int CORE_cpamm(int op, PLASMA_enum side, PLASMA_enum storev,
               int M, int N, int K, int L,
               const PLASMA_Complex32_t *A1, int LDA1,
                     PLASMA_Complex32_t *A2, int LDA2,
               const PLASMA_Complex32_t *V, int LDV,
                     PLASMA_Complex32_t *W, int LDW);
int  CORE_cparfb(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                 int M1, int N1, int M2, int N2, int K, int L,
                       PLASMA_Complex32_t *A1, int LDA1,
                       PLASMA_Complex32_t *A2, int LDA2,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                       PLASMA_Complex32_t *WORK, int LDWORK);
int CORE_cpemv(PLASMA_enum trans, PLASMA_enum storev,
               int M, int N, int L,
               PLASMA_Complex32_t ALPHA,
               const PLASMA_Complex32_t *A, int LDA,
               const PLASMA_Complex32_t *X, int INCX,
               PLASMA_Complex32_t BETA,
               PLASMA_Complex32_t *Y, int INCY,
               PLASMA_Complex32_t *WORK);
void CORE_cplghe(float bump, int m, int n, PLASMA_Complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_cplgsy(PLASMA_Complex32_t bump, int m, int n, PLASMA_Complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_cplrnt(int m, int n, PLASMA_Complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
int  CORE_cpltmg(PLASMA_enum mtxtype, int m, int n, PLASMA_Complex32_t *A, int lda,
                  int gM, int gN, int m0, int n0, unsigned long long int seed );
int  CORE_cpltmg_chebvand( int M, int N, PLASMA_Complex32_t *A, int LDA,
                           int gN, int m0, int n0,
                           PLASMA_Complex32_t *W );
int  CORE_cpltmg_circul( int M, int N, PLASMA_Complex32_t *A, int LDA,
                         int gM, int m0, int n0,
                         const PLASMA_Complex32_t *V );
void CORE_cpltmg_condexq( int M, int N, PLASMA_Complex32_t *Q, int LDQ );
void CORE_cpltmg_fiedler(int m, int n,
                         const PLASMA_Complex32_t *X, int incX,
                         const PLASMA_Complex32_t *Y, int incY,
                         PLASMA_Complex32_t *A, int lda);
int  CORE_cpltmg_hankel( PLASMA_enum uplo, int M, int N, PLASMA_Complex32_t *A, int LDA,
                         int m0, int n0, int nb,
                         const PLASMA_Complex32_t *V1,
                         const PLASMA_Complex32_t *V2 );
void CORE_cpltmg_toeppd1( int gM, int m0, int M, PLASMA_Complex32_t *W,
                          unsigned long long int seed );
void CORE_cpltmg_toeppd2( int M, int N, int K, int m0, int n0,
                          const PLASMA_Complex32_t *W,
                          PLASMA_Complex32_t *A, int LDA );
void CORE_cpotrf(PLASMA_enum uplo, int N, PLASMA_Complex32_t *A, int LDA, int *INFO);
void CORE_csetvar(const PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *x);
void CORE_cshift(int s, int m, int n, int L,
                 PLASMA_Complex32_t *A);
void CORE_cshiftw(int s, int cl, int m, int n, int L,
                  PLASMA_Complex32_t *A, PLASMA_Complex32_t *W);
int  CORE_cssssm(int M1, int N1, int M2, int N2, int K, int IB,
                       PLASMA_Complex32_t *A1, int LDA1,
                       PLASMA_Complex32_t *A2, int LDA2,
                 const PLASMA_Complex32_t *L1, int LDL1,
                 const PLASMA_Complex32_t *L2, int LDL2,
                 const int *IPIV);
int CORE_cstedc(PLASMA_enum compz, int n,
                float *D, float *E,
                PLASMA_Complex32_t *Z, int LDZ,
                PLASMA_Complex32_t *WORK, int LWORK,
#ifdef COMPLEX
                float *RWORK, int LRWORK,
#endif
                int *IWORK, int LIWORK);
int CORE_csteqr(PLASMA_enum compz, int n,
                float *D, float *E,
                PLASMA_Complex32_t *Z, int LDZ,
                float *WORK);
void CORE_csymm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                                          const PLASMA_Complex32_t *B, int LDB,
                PLASMA_Complex32_t beta,        PLASMA_Complex32_t *C, int LDC);
void CORE_csyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t beta,        PLASMA_Complex32_t *C, int LDC);
void CORE_csyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                                           const PLASMA_Complex32_t *B, int LDB,
                 PLASMA_Complex32_t beta,        PLASMA_Complex32_t *C, int LDC);
int  CORE_csyssq(PLASMA_enum uplo, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                 float *scale, float *sumsq);
void CORE_cswpab(int i, int n1, int n2,
                 PLASMA_Complex32_t *A, PLASMA_Complex32_t *work);
int  CORE_cswptr_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc,
                        const PLASMA_Complex32_t *Akk, int ldak);
int CORE_ctradd(PLASMA_enum uplo, PLASMA_enum trans, int M, int N,
                      PLASMA_Complex32_t alpha,
                const PLASMA_Complex32_t *A, int LDA,
                      PLASMA_Complex32_t beta,
                      PLASMA_Complex32_t *B, int LDB);
void CORE_ctrasm(PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N, const PLASMA_Complex32_t *A, int lda, float *work);
void CORE_ctrdalg1(int n,
                        int nb,
                        PLASMA_Complex32_t *A,
                        int lda,
                        PLASMA_Complex32_t *V,
                        PLASMA_Complex32_t *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        PLASMA_Complex32_t *work);
void CORE_ctrmm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                                                PLASMA_Complex32_t *B, int LDB);
void CORE_ctrsm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int LDA,
                                                PLASMA_Complex32_t *B, int LDB);
int  CORE_ctrssq(PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                 const PLASMA_Complex32_t *A, int LDA,
                 float *scale, float *sumsq);
void CORE_ctrtri(PLASMA_enum uplo, PLASMA_enum diag, int N,
                 PLASMA_Complex32_t *A, int LDA, int *info);
int  CORE_ctslqt(int M, int N, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU, PLASMA_Complex32_t *WORK);
int  CORE_ctsmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int CORE_ctsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        PLASMA_Complex32_t *A3, int lda3,
                        const PLASMA_Complex32_t *V, int ldv,
                        const PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork);
int CORE_ctsmlq_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        const PLASMA_Complex32_t *V, int ldv,
                        const PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork);
int  CORE_ctsmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int CORE_ctsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        PLASMA_Complex32_t *A3, int lda3,
                        const PLASMA_Complex32_t *V, int ldv,
                        const PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork);
int CORE_ctsmqr_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        const PLASMA_Complex32_t *V, int ldv,
                        const PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork);
int  CORE_ctsqrt(int M, int N, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU, PLASMA_Complex32_t *WORK);
int  CORE_ctstrf(int M, int N, int IB, int NB,
                 PLASMA_Complex32_t *U, int LDU,
                 PLASMA_Complex32_t *A, int LDA,
                 PLASMA_Complex32_t *L, int LDL,
                 int *IPIV, PLASMA_Complex32_t *WORK,
                 int LDWORK, int *INFO);
int  CORE_cttmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_cttqrt(int M, int N, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU,
                 PLASMA_Complex32_t *WORK);
int  CORE_cttmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_cttlqt(int M, int N, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *TAU,
                 PLASMA_Complex32_t *WORK);
int  CORE_cunmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *C, int LDC,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_cunmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *C, int LDC,
                 PLASMA_Complex32_t *WORK, int LDWORK);

#ifndef COMPLEX
void CORE_slaed2_computeK(int *K, int n, int n1,
                          float *beta, float *D, float *Q, int LDQ,
                          float *Z, float *DLAMBDA, float *W,
                          int *INDX, int *INDXC, int *INDXP, int *INDXQ,
                          int *COLTYP);
void CORE_slaed2_compressq(int n, int n1, const int *INDX, const int *ctot,
                           const float *Q, int LDQ, float *Q2,
                           int start, int end);
void CORE_slaed2_copydef(int n, int n1, int K, const int *ctot,
                         float *Q, int LDQ, const float *Q2,
                         int start, int end);
int CORE_slaed4(int n, int K,
                float *D, float beta,
                float *Q, int LDQ,
                const float *D0, const float *Z,
                const int *INDX,
                int start, int end );
void CORE_slaed3_computeW(int n, int K,
                          const float *Q, int LDQ,
                          const float *DLAMBDA, float *W,
                          const int *INDX,
                          int start, int end);
void CORE_slaed3_reduceW(int n, int n1, int K, int l,
                         const float *Q, int LDQ,
                         const float *Wred, float *W);
void CORE_slaed3_computevectors(int K, int il_nondef, int iu_nondef,
                                float *Q, int LDQ, float *W, float *S,
                                const int *INDXC,
                                int start, int end);
void CORE_slaed3_merge( int n, int K, float *D, int *INDXQ );
void CORE_slaed3_updatevectors(int op, int wsmode, int n, int n1, int K,
                   int il_nondef, int iu_nondef,
                               float *Q, int ldq, float *Q2,
                               const int *ctot, float *WORK, int start, int end);
#endif
void CORE_cswap(int m, int n, PLASMA_Complex32_t *Q, int ldq,
                const PLASMA_Complex32_t *work, const int *perm,
                int start, int end);
int CORE_clascl(PLASMA_enum type, int kl, int ku, float cfrom, float cto,
                int m, int n, PLASMA_Complex32_t *A, int lda);
#ifdef COMPLEX
int CORE_slag2c(int m, int n, const float *Q, int LDQ,
                 PLASMA_Complex32_t *Z, int LDZ);
#endif

#ifndef COMPLEX
void CORE_slaed3_freebigwork(int oper, float **WORK);
void CORE_slaed0_betaapprox(int subpbs, const int *subpbs_info,
                            float *D, const float *E);
int CORE_slapst(PLASMA_enum type, int n,
                const float *D, int *INDX);
#endif

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
