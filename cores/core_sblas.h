/**
 * Copyright (c) 2019      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Imported from:
 *
 * @file core_sblas.h
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
 * @generated s Fri Apr  1 11:02:20 2016
 *
 **/
#ifndef _PLASMA_CORE_SBLAS_H_
#define _PLASMA_CORE_SBLAS_H_

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

struct CORE_sgetrf_data_s;
typedef struct CORE_sgetrf_data_s CORE_sgetrf_data_t;

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/
void CORE_sasum(int storev, PLASMA_enum uplo, int M, int N,
                 const float *A, int lda, float *work);
void CORE_sbrdalg1( PLASMA_enum uplo,
                    int n,
                    int nb,
                    float *A,
                    int lda,
                    float *VQ,
                    float *TAUQ,
                    float *VP,
                    float *TAUP,
                    int Vblksiz, int wantz,
                    int i, int sweepid, int m, int grsiz,
                    float *work);
int CORE_sgbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, float *V, float *TAU,
                int st, int ed, int eltsize);
int CORE_sgbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, float *V, float *TAU,
                int st, int ed, int eltsize);
int CORE_sgblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, float *V, float *TAU,
                int st, int ed, int eltsize);
int CORE_sgeadd(PLASMA_enum trans, int M, int N,
                      float alpha,
                const float *A, int LDA,
                      float beta,
                      float *B, int LDB);
int  CORE_sgelqt(int M, int N, int IB,
                 float *A, int LDA,
                 float *T, int LDT,
                 float *TAU,
                 float *WORK);
void CORE_sgemm(PLASMA_enum transA, PLASMA_enum transB,
                int M, int N, int K,
                float alpha, const float *A, int LDA,
                                          const float *B, int LDB,
                float beta,        float *C, int LDC);
void CORE_sgemv(PLASMA_enum trans, int M, int N,
                float alpha, const float *A, int LDA,
                                          const float *x, int incx,
                float beta,        float *y, int incy);
void CORE_sgeqp3_init( int n, int *jpvt );
void CORE_sgeqp3_larfg( PLASMA_desc A, int ii, int jj, int i, int j,
                        float *tau, float *beta );
void CORE_sgeqp3_norms( PLASMA_desc A, int ioff, int joff, float *norms1, float *norms2 );
void CORE_sgeqp3_pivot( PLASMA_desc A, float *F, int ldf,
                        int jj, int k, int *jpvt,
                        float *norms1, float *norms2, int *info );
int  CORE_sgeqp3_tntpiv(int m, int n,
                        float *A, int lda,
                        int *IPIV, float *tau,
                        int *iwork);
void CORE_sgeqp3_update( const float *Ajj, int lda1,
                         float       *Ajk, int lda2,
                         const float *Fk,  int ldf,
                         int joff, int k, int koff, int nb,
                         float *norms1, float *norms2,
                         int *info );
int  CORE_sgeqrt(int M, int N, int IB,
                 float *A, int LDA,
                 float *T, int LDT,
                 float *TAU, float *WORK);
int  CORE_sgessm(int M, int N, int K, int IB,
                 const int *IPIV,
                 const float *L, int LDL,
                 float *A, int LDA);
int  CORE_sgessq(int M, int N,
                 const float *A, int LDA,
                 float *scale, float *sumsq);
int  CORE_sgetf2_nopiv(int m, int n,
                      float *A, int lda);
int  CORE_sgetrf(int M, int N,
                 float *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_sgetrf_incpiv(int M, int N, int IB,
                        float *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_sgetrf_nopiv(int m, int n, int ib,
                      float *A, int lda);
int  CORE_sgetrf_reclap(CORE_sgetrf_data_t *data, int M, int N,
                        float *A, int LDA,
                        int *IPIV, int *info);
CORE_sgetrf_data_t *CORE_sgetrf_reclap_init(int nbthrd);
int  CORE_sgetrf_rectil(CORE_sgetrf_data_t *data, const PLASMA_desc A, int *IPIV, int *info);
CORE_sgetrf_data_t *CORE_sgetrf_rectil_init(int nbthrd);
void CORE_sgetrip(int m, int n, float *A,
                  float *work);
int CORE_shbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, float *V, float *TAU,
                int st, int ed, int eltsize);
int CORE_shblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, float *V, float *TAU,
                int st, int ed, int eltsize);
int CORE_shbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, float *V, float *TAU,
                int st, int ed, int eltsize);
void CORE_ssbtype1cb(int N, int NB,
                     float *A, int LDA,
                     float *V, float *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     float *WORK);
void CORE_ssbtype2cb(int N, int NB,
                     float *A, int LDA,
                     float *V, float *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     float *WORK);
void CORE_ssbtype3cb(int N, int NB,
                     float *A, int LDA,
                     const float *V, const float *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     float *WORK);
void CORE_sgbtype1cb(PLASMA_enum uplo, int N, int NB,
                float *A, int LDA,
                float *VQ, float *TAUQ,
                float *VP, float *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                float *WORK);
void CORE_sgbtype2cb(PLASMA_enum uplo, int N, int NB,
                float *A, int LDA,
                float *VQ, float *TAUQ,
                float *VP, float *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                float *WORK);
void CORE_sgbtype3cb(PLASMA_enum uplo, int N, int NB,
                float *A, int LDA,
                float *VQ, float *TAUQ,
                float *VP, float *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                float *WORK);
void CORE_ssygst(int itype, PLASMA_enum uplo, int N,
                 float *A, int LDA,
                 float *B, int LDB, int *INFO);
#ifdef COMPLEX
void CORE_ssymm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                float alpha, const float *A, int LDA,
                                          const float *B, int LDB,
                float beta,        float *C, int LDC);
void CORE_ssyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                float alpha, const float *A, int LDA,
                float beta,        float *C, int LDC);
void CORE_ssyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 float alpha, const float *A, int LDA,
                                           const float *B, int LDB,
                 float beta,                    float *C, int LDC);
int  CORE_shessq(PLASMA_enum uplo, int N,
                 const float *A, int LDA,
                 float *scale, float *sumsq);
#endif
int  CORE_ssyrfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 const float *A,    int LDA,
                 const float *T,    int LDT,
                       float *C,    int LDC,
                       float *WORK, int LDWORK);
void CORE_slacpy(PLASMA_enum uplo, int M, int N,
                 const float *A, int LDA,
                       float *B, int LDB);
int CORE_slacpy_pivot( const PLASMA_desc descA,
                       PLASMA_enum direct,
                       int k1, int k2, const int *ipiv,
                       int *rankin, int *rankout,
                       float *A, int lda,
                       int init);
void CORE_slange(int norm, int M, int N,
                 const float *A, int LDA,
                 float *work, float *normA);
#ifdef COMPLEX
void CORE_slansy(int norm, PLASMA_enum uplo, int N,
                 const float *A, int LDA,
                 float *work, float *normA);
#endif
void CORE_slansy(int norm, PLASMA_enum uplo, int N,
                 const float *A, int LDA,
                 float *work, float *normA);
void CORE_slantr(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N,
                 const float *A, int LDA,
                 float *work, float *normA);
int CORE_slarfb_gemm(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                     int M, int N, int K,
                     const float *V, int LDV,
                     const float *T, int LDT,
                           float *C, int LDC,
                           float *WORK, int LDWORK);
int CORE_slarfx2(PLASMA_enum side, int N,
                 float V,
                 float TAU,
                 float *C1, int LDC1,
                 float *C2, int LDC2);
int CORE_slarfx2c(PLASMA_enum uplo,
                  float V,
                  float TAU,
                  float *C1,
                  float *C2,
                  float *C3);
int CORE_slarfx2ce(PLASMA_enum uplo,
                   float *V,
                   float *TAU,
                   float *C1,
                   float *C2,
                   float *C3);
void CORE_slarfy(int N,
                 float *A, int LDA,
                 const float *V,
                 const float *TAU,
                 float *WORK);
int  CORE_slascal(PLASMA_enum uplo, int m, int n,
                  float alpha, float *A, int lda);
void CORE_slaset(PLASMA_enum uplo, int n1, int n2,
                 float alpha, float beta,
                 float *tileA, int ldtilea);
void CORE_slaset2(PLASMA_enum uplo, int n1, int n2, float alpha,
                  float *tileA, int ldtilea);
void CORE_slaswp(int N, float *A, int LDA,
                 int I1,  int I2, const int *IPIV, int INC);
int  CORE_slaswp_ontile( PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_slaswpc_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_slatro(PLASMA_enum uplo, PLASMA_enum trans,
                 int M, int N,
                 const float *A, int LDA,
                       float *B, int LDB);
void CORE_slauum(PLASMA_enum uplo, int N, float *A, int LDA);
int CORE_spamm(int op, PLASMA_enum side, PLASMA_enum storev,
               int M, int N, int K, int L,
               const float *A1, int LDA1,
                     float *A2, int LDA2,
               const float *V, int LDV,
                     float *W, int LDW);
int  CORE_sparfb(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                 int M1, int N1, int M2, int N2, int K, int L,
                       float *A1, int LDA1,
                       float *A2, int LDA2,
                 const float *V, int LDV,
                 const float *T, int LDT,
                       float *WORK, int LDWORK);
int CORE_spemv(PLASMA_enum trans, PLASMA_enum storev,
               int M, int N, int L,
               float ALPHA,
               const float *A, int LDA,
               const float *X, int INCX,
               float BETA,
               float *Y, int INCY,
               float *WORK);
void CORE_splgsy(float bump, int m, int n, float *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_splgsy(float bump, int m, int n, float *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_splrnt(int m, int n, float *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
int  CORE_spltmg(PLASMA_enum mtxtype, int m, int n, float *A, int lda,
                  int gM, int gN, int m0, int n0, unsigned long long int seed );
int  CORE_spltmg_chebvand( int M, int N, float *A, int LDA,
                           int gN, int m0, int n0,
                           float *W );
int  CORE_spltmg_circul( int M, int N, float *A, int LDA,
                         int gM, int m0, int n0,
                         const float *V );
void CORE_spltmg_condexq( int M, int N, float *Q, int LDQ );
void CORE_spltmg_fiedler(int m, int n,
                         const float *X, int incX,
                         const float *Y, int incY,
                         float *A, int lda);
int  CORE_spltmg_hankel( PLASMA_enum uplo, int M, int N, float *A, int LDA,
                         int m0, int n0, int nb,
                         const float *V1,
                         const float *V2 );
void CORE_spltmg_toeppd1( int gM, int m0, int M, float *W,
                          unsigned long long int seed );
void CORE_spltmg_toeppd2( int M, int N, int K, int m0, int n0,
                          const float *W,
                          float *A, int LDA );
void CORE_spotrf(PLASMA_enum uplo, int N, float *A, int LDA, int *INFO);
void CORE_ssetvar(const float *alpha, float *x);
void CORE_sshift(int s, int m, int n, int L,
                 float *A);
void CORE_sshiftw(int s, int cl, int m, int n, int L,
                  float *A, float *W);
int  CORE_sssssm(int M1, int N1, int M2, int N2, int K, int IB,
                       float *A1, int LDA1,
                       float *A2, int LDA2,
                 const float *L1, int LDL1,
                 const float *L2, int LDL2,
                 const int *IPIV);
int CORE_sstedc(PLASMA_enum compz, int n,
                float *D, float *E,
                float *Z, int LDZ,
                float *WORK, int LWORK,
#ifdef COMPLEX
                float *RWORK, int LRWORK,
#endif
                int *IWORK, int LIWORK);
int CORE_ssteqr(PLASMA_enum compz, int n,
                float *D, float *E,
                float *Z, int LDZ,
                float *WORK);
void CORE_ssymm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                float alpha, const float *A, int LDA,
                                          const float *B, int LDB,
                float beta,        float *C, int LDC);
void CORE_ssyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                float alpha, const float *A, int LDA,
                float beta,        float *C, int LDC);
void CORE_ssyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 float alpha, const float *A, int LDA,
                                           const float *B, int LDB,
                 float beta,        float *C, int LDC);
int  CORE_ssyssq(PLASMA_enum uplo, int N,
                 const float *A, int LDA,
                 float *scale, float *sumsq);
void CORE_sswpab(int i, int n1, int n2,
                 float *A, float *work);
int  CORE_sswptr_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc,
                        const float *Akk, int ldak);
int CORE_stradd(PLASMA_enum uplo, PLASMA_enum trans, int M, int N,
                      float alpha,
                const float *A, int LDA,
                      float beta,
                      float *B, int LDB);
void CORE_strasm(PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N, const float *A, int lda, float *work);
void CORE_strdalg1(int n,
                        int nb,
                        float *A,
                        int lda,
                        float *V,
                        float *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        float *work);
void CORE_strmm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                float alpha, const float *A, int LDA,
                                                float *B, int LDB);
void CORE_strsm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                float alpha, const float *A, int LDA,
                                                float *B, int LDB);
int  CORE_strssq(PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                 const float *A, int LDA,
                 float *scale, float *sumsq);
void CORE_strtri(PLASMA_enum uplo, PLASMA_enum diag, int N,
                 float *A, int LDA, int *info);
int  CORE_stslqt(int M, int N, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *T, int LDT,
                 float *TAU, float *WORK);
int  CORE_stsmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 const float *V, int LDV,
                 const float *T, int LDT,
                 float *WORK, int LDWORK);
int CORE_stsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        float *A1, int lda1,
                        float *A2, int lda2,
                        float *A3, int lda3,
                        const float *V, int ldv,
                        const float *T, int ldt,
                        float *WORK, int ldwork);
int CORE_stsmlq_sytra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        float *A1, int lda1,
                        float *A2, int lda2,
                        const float *V, int ldv,
                        const float *T, int ldt,
                        float *WORK, int ldwork);
int  CORE_stsmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 const float *V, int LDV,
                 const float *T, int LDT,
                 float *WORK, int LDWORK);
int CORE_stsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        float *A1, int lda1,
                        float *A2, int lda2,
                        float *A3, int lda3,
                        const float *V, int ldv,
                        const float *T, int ldt,
                        float *WORK, int ldwork);
int CORE_stsmqr_sytra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        float *A1, int lda1,
                        float *A2, int lda2,
                        const float *V, int ldv,
                        const float *T, int ldt,
                        float *WORK, int ldwork);
int  CORE_stsqrt(int M, int N, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *T, int LDT,
                 float *TAU, float *WORK);
int  CORE_ststrf(int M, int N, int IB, int NB,
                 float *U, int LDU,
                 float *A, int LDA,
                 float *L, int LDL,
                 int *IPIV, float *WORK,
                 int LDWORK, int *INFO);
int  CORE_sttmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 const float *V, int LDV,
                 const float *T, int LDT,
                 float *WORK, int LDWORK);
int  CORE_sttqrt(int M, int N, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *T, int LDT,
                 float *TAU,
                 float *WORK);
int  CORE_sttmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 const float *V, int LDV,
                 const float *T, int LDT,
                 float *WORK, int LDWORK);
int  CORE_sttlqt(int M, int N, int IB,
                 float *A1, int LDA1,
                 float *A2, int LDA2,
                 float *T, int LDT,
                 float *TAU,
                 float *WORK);
int  CORE_sormlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const float *V, int LDV,
                 const float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);
int  CORE_sormqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const float *V, int LDV,
                 const float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);

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
void CORE_sswap(int m, int n, float *Q, int ldq,
                const float *work, const int *perm,
                int start, int end);
int CORE_slascl(PLASMA_enum type, int kl, int ku, float cfrom, float cto,
                int m, int n, float *A, int lda);
#ifdef COMPLEX
int CORE_slag2c(int m, int n, const float *Q, int LDQ,
                 float *Z, int LDZ);
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
