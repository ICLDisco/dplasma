/**
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

#if defined(QUARK_H)
/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by PLASMA) - alphabetical order
 **/
void QUARK_CORE_sasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                       const float *A, int lda, int szeA,
                       float *work, int szeW);
void QUARK_CORE_sasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                          const float *A, int lda, int szeA,
                          float *work, int szeW,
                          float *fake, int szeF);
void QUARK_CORE_sgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum trans, int m, int n, int nb,
                       float  alpha,
                       const float *A, int lda,
                       float  beta,
                       float *B, int ldb);
void QUARK_CORE_sbrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo,
                        int n, int nb,
                        float *A,
                        int lda,
                        float *VQ,
                        float *TAUQ,
                        float *VP,
                        float *TAUP,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_sgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt);
void QUARK_CORE_sgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int k, int nb,
                      float alpha, const float *A, int lda,
                      const float *B, int ldb,
                      float beta, float *C, int ldc);
void QUARK_CORE_sgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum transA, PLASMA_enum transB,
                        int m, int n, int k, int nb,
                        float alpha, const float *A, int lda,
                        const float *B, int ldb,
                        float beta, float *C, int ldc);
void QUARK_CORE_sgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         float alpha, const float *A, int lda,
                         const float *B, int ldb,
                         float beta, float *C, int ldc,
                         float *fake1, int szefake1, int flag1,
                         float *fake2, int szefake2, int flag2);
void QUARK_CORE_sgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         float alpha, const float *A, int lda,
                         const float **B, int ldb,
                         float beta, float *C, int ldc);
void QUARK_CORE_sgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           float alpha, const float *A, int lda,
                           const float **B, int ldb,
                           float beta, float *C, int ldc,
                           float *fake1, int szefake1, int flag1);
void QUARK_CORE_sgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         float alpha, const float *A, int lda,
                         const float *B, int ldb,
                         float beta, float **C, int ldc);
void QUARK_CORE_sgemm_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           const float *alpha, const float *A, int lda,
                                                            const float *B, int ldb,
                           const float *beta,        float *C, int ldc,
                           const float *Alock,
                           const float *Block,
                           const float *Clock);
void QUARK_CORE_sgemv(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum trans, int m, int n,
                      float alpha, const float *A, int lda,
                                                const float *x, int incx,
                      float beta,        float *y, int incy);
void QUARK_CORE_sgemv_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum trans,
                           int m, int n,
                           const float *alpha, const float *A, int lda,
                                                            const float *x, int incx,
                           const float *beta,        float *y, int incy,
                           const float *Alock,
                           const float *xlock,
                           const float *ylock);
void QUARK_CORE_sgeqp3_init( Quark *quark, Quark_Task_Flags *task_flags,
                             int n, int *jpvt );
void QUARK_CORE_sgeqp3_larfg(Quark *quark, Quark_Task_Flags *task_flags,
                             PLASMA_desc A, int ii, int jj, int i, int j,
                             float *tau, float *beta );
void QUARK_CORE_sgeqp3_norms( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, int ioff, int joff, float *norms1, float *norms2 );
void QUARK_CORE_sgeqp3_pivot( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A,
                              float *F, int ldf,
                              int jj, int k, int *jpvt,
                              float *norms1, float *norms2, int *info );
void QUARK_CORE_sgeqp3_tntpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              float *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_sgeqp3_update( Quark *quark, Quark_Task_Flags *task_flags,
                               float *Ajj, int lda1,
                               float *Ajk, int lda2,
                               float *Fk,  int ldf,
                               int joff, int k, int koff, int nb,
                               float *norms1, float *norms2, int *info );
void QUARK_CORE_sgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A, int lda,
                       float *T, int ldt);
void QUARK_CORE_sgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       const int *IPIV,
                       const float *L, int ldl,
                       float *A, int lda);
void QUARK_CORE_sgessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, const float *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF );
void QUARK_CORE_sgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       float *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_sgetrf_incpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int ib, int nb,
                              float *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_sgetrf_nopiv(Quark *quark, Quark_Task_Flags *task_flags,
                             int m, int n, int ib, int nb,
                             float *A, int lda,
                             PLASMA_sequence *sequence, PLASMA_request *request,
                             int iinfo);
void QUARK_CORE_sgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_sgetrf_data_t *data, int m, int n, int nb,
                              float *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_sgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_sgetrf_data_t *data,
                              PLASMA_desc A, float *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_sgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, float *A, int szeA);
void QUARK_CORE_sgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, float *A, int szeA,
                           float *fake, int szeF, int paramF);
void QUARK_CORE_sgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, float *A, int szeA,
                           float *fake1, int szeF1, int paramF1,
                           float *fake2, int szeF2, int paramF2);
void QUARK_CORE_ssymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      float alpha, const float *A, int lda,
                      const float *B, int ldb,
                      float beta, float *C, int ldc);
void QUARK_CORE_ssygst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, PLASMA_enum uplo, int N,
                       float *A, int LDA,
                       float *B, int LDB,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_ssyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      float alpha, const float *A, int lda,
                      float beta, float *C, int ldc);
void QUARK_CORE_ssyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       float alpha, const float *A, int lda,
                       const float *B, int LDB,
                       float beta, float *C, int ldc);
void QUARK_CORE_ssyrfb(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo,
                       int n, int k, int ib, int nb,
                       const float *A, int lda,
                       const float *T, int ldt,
                       float *C, int ldc);
void QUARK_CORE_shessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const float *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF );
void QUARK_CORE_slacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int mb,
                       const float *A, int lda,
                       float *B, int ldb);
void QUARK_CORE_slacpy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, int m, int n, int nb,
                          const float *A, int lda,
                          float *B, int ldb,
                          float *fake1, int szefake1, int flag1);
void QUARK_CORE_slacpy_pivot(Quark *quark, Quark_Task_Flags *task_flags,
                             const PLASMA_desc descA,
                             PLASMA_enum direct,
                             int k1, int k2, const int *ipiv,
                             int *rankin, int *rankout,
                             float *A, int lda,
                             int pos, int init);
void QUARK_CORE_slange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       const float *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_slange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int M, int N,
                          const float *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
#ifdef COMPLEX
void QUARK_CORE_slansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const float *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_slansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, PLASMA_enum uplo, int N,
                          const float *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
#endif
void QUARK_CORE_slansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const float *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_slansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, PLASMA_enum uplo, int N,
                          const float *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
void QUARK_CORE_slantr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                       const float *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_slantr_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                          const float *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
void QUARK_CORE_slascal(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int m, int n, int nb,
                        float alpha, float *A, int lda);
void QUARK_CORE_slaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n1, int n2, float alpha,
                       float beta, float *tileA, int ldtilea);
void QUARK_CORE_slaset2(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int n1, int n2, float alpha,
                        float *tileA, int ldtilea);
void QUARK_CORE_slaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, float *A, int lda,
                       int i1,  int i2, const int *ipiv, int inc);
void QUARK_CORE_slaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, float *A, int lda,
                          int i1,  int i2, const int *ipiv, int inc,
                          float *fake1, int szefake1, int flag1,
                          float *fake2, int szefake2, int flag2);
void QUARK_CORE_slaswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, float *A,
                              int i1,  int i2, const int *ipiv, int inc, float *fakepanel);
void QUARK_CORE_slaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, float *A,
                                 int i1,  int i2, const int *ipiv, int inc,
                                 float *fake1, int szefake1, int flag1,
                                 float *fake2, int szefake2, int flag2);
void QUARK_CORE_slaswpc_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_desc descA, float *A,
                               int i1,  int i2, const int *ipiv, int inc, float *fakepanel);
void QUARK_CORE_slatro(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int mb,
                       const float *A, int lda,
                       float *B, int ldb);
void QUARK_CORE_slatro_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int mb,
                          const float *A, int lda,
                                float *B, int ldb,
                          float *fake1, int szefake1, int flag1);
void QUARK_CORE_slauum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       float *A, int lda);
void QUARK_CORE_splgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       float bump, int m, int n, float *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_splgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       float bump, int m, int n, float *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_splrnt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, float *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_spltmg(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum mtxtype, int m, int n, float *A, int lda,
                        int gM, int gN, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_spltmg_chebvand( Quark *quark, Quark_Task_Flags *task_flags,
                                 int M, int N, float *A, int LDA,
                                 int gN, int m0, int n0,
                                 float *W );
void QUARK_CORE_spltmg_circul( Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, float *A, int LDA,
                               int gM, int m0, int n0,
                               const float *W );
void QUARK_CORE_spltmg_fiedler(Quark *quark, Quark_Task_Flags *task_flags,
                               int m, int n,
                               const float *X, int incX,
                               const float *Y, int incY,
                               float *A, int lda);
void QUARK_CORE_spltmg_hankel( Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_enum uplo, int M, int N, float *A, int LDA,
                               int m0, int n0, int nb,
                               const float *V1,
                               const float *V2);
void QUARK_CORE_spltmg_toeppd1(Quark *quark, Quark_Task_Flags *task_flags,
                               int gM, int m0, int M,
                               float *W,
                               unsigned long long int seed);
void QUARK_CORE_spltmg_toeppd2(Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, int K, int m0, int n0,
                               const float *W,
                               float *A, int LDA );
void QUARK_CORE_spotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       float *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_ssetvar(Quark *quark, Quark_Task_Flags *task_flags,
                        const float *alpha, float *x,
                        float *Alock);
void QUARK_CORE_sshift( Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int m, int n, int L,
                        float *A);
void QUARK_CORE_sshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L,
                        float *A, float *W);
void QUARK_CORE_sssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       const float *L1, int ldl1,
                       const float *L2, int ldl2,
                       const int *IPIV);
void QUARK_CORE_sstedc(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       float *D, float *E,
                       float *Z, int ldz);
void QUARK_CORE_sstedc_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum compz, int n,
                          float *D, float *E,
                          float *Z, int ldz,
                          void *fake1, int szefake1, int flag1,
                          void *fake2, int szefake2, int flag2);
void QUARK_CORE_ssteqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       float *D, float *E,
                       float *Z, int ldz);
void QUARK_CORE_ssymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      float alpha, const float *A, int lda,
                      const float *B, int ldb,
                      float beta, float *C, int ldc);
void QUARK_CORE_ssyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      float alpha, const float *A, int lda,
                      float beta, float *C, int ldc);
void QUARK_CORE_ssyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       float alpha, const float *A, int lda,
                       const float *B, int LDB,
                       float beta, float *C, int ldc);
void QUARK_CORE_ssyssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const float *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF );
void QUARK_CORE_sswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       float *A, int szeA);
void QUARK_CORE_sswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, float *Aij,
                              int i1,  int i2, const int *ipiv, int inc,
                              const float *Akk, int ldak);
void QUARK_CORE_stradd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       float  alpha,
                       const float *A, int lda,
                       float  beta,
                       float *B, int ldb);
void QUARK_CORE_strasm(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int m, int n,
                       const float *A, int lda, int szeA,
                       float *work, int szeW);
void QUARK_CORE_strasm_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int m, int n,
                          const float *A, int lda, int szeA,
                          float *work, int szeW,
                          float *fake, int szeF);
void QUARK_CORE_strdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                        int n,
                        int nb,
                        float *A,
                        int lda,
                        float *V,
                        float *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_strmm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      float alpha, const float *A, int lda,
                      float *B, int ldb);
void QUARK_CORE_strmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                         int m, int n, int nb,
                         float alpha, const float *A, int lda,
                         float **B, int ldb);
void QUARK_CORE_strsm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      float alpha, const float *A, int lda,
                      float *B, int ldb);
void QUARK_CORE_strssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, PLASMA_enum diag,
                           int m, int n, const float *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF );
void QUARK_CORE_strtri(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum diag, int n, int nb,
                       float *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_stslqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt);
void QUARK_CORE_stsmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       const float *V, int ldv,
                       const float *T, int ldt);
void QUARK_CORE_stsmlq_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              const float *V, int ldv,
                              const float *T, int ldt);
void QUARK_CORE_stsmlq_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              float *A3, int lda3,
                              const float *V, int ldv,
                              const float *T, int ldt);
void QUARK_CORE_stsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       const float *V, int ldv,
                       const float *T, int ldt);
void QUARK_CORE_stsmqr_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              const float *V, int ldv,
                              const float *T, int ldt);
void QUARK_CORE_stsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              float *A3, int lda3,
                              const float *V, int ldv,
                              const float *T, int ldt);
void QUARK_CORE_stsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt);
void QUARK_CORE_ststrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *U, int ldu,
                       float *A, int lda,
                       float *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_sttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       const float *V, int ldv,
                       const float *T, int ldt);
void QUARK_CORE_sttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt);
void QUARK_CORE_sttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       const float *V, int ldv,
                       const float *T, int ldt);
void QUARK_CORE_sttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt);
void QUARK_CORE_spamm(Quark *quark, Quark_Task_Flags *task_flags,
                       int op, PLASMA_enum side, PLASMA_enum storev,
                       int m, int n, int k, int l,
                       const float *A1, int lda1,
                       float *A2, int lda2,
                       const float *V, int ldv,
                       float *W, int ldw);
void QUARK_CORE_splssq( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, const float *A, float *result );
void QUARK_CORE_sormlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int ib,  int nb, int k,
                       const float *A, int lda,
                       const float *T, int ldt,
                       float *C, int ldc);
void QUARK_CORE_sormqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const float *A, int lda,
                       const float *T, int ldt,
                       float *C, int ldc);


void QUARK_CORE_slascl(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum type, int kl, int ku, float cfrom, float cto,
                       int m, int n, float *A, int lda);
void QUARK_CORE_slascl_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                            PLASMA_enum type, int kl, int ku, float *cfrom, float *cto,
                            int m, int n, float *A, int lda,
                            void *fake, int szefake, int flag);
void QUARK_CORE_slaed0_lascl( Quark *quark, Quark_Task_Flags *task_flags,
                              int n, float *scale, float *D, float *E);
void QUARK_CORE_slaed0_betaapprox(Quark *quark, Quark_Task_Flags *task_flags,
                                  int subpbs, const int *subpbs_info,
                                  float *D, const float *E);

#ifndef COMPLEX
void QUARK_CORE_slaed2_computeK(Quark *quark, Quark_Task_Flags *task_flags,
                                int *K1, int n, int n1,
                                float *beta, float *D, float *Q, int LDQ,
                                float *Z, float *DLAMBDA, float *W,
                                int *INDX, int *INDXC, int *INDXP, int *INDXQ,
                                int *COLTYP,
                                float **Qmerge, int wsmode,
                                int *K2);

void QUARK_CORE_slaed1_pipelined(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, const int *K,
                                 const int *INDX, const int *ctot,
                                 float *D, const float *beta,
                                 float *Q, int LDQ, float *Q2,
                                 const float *DLAMBDA, const float *W, float *Wred,
                                 int start, int end);
void QUARK_CORE_slaed2_compressq(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, int start, int end,
                                 const int *INDX, const int *ctot,
                                 const float *Q, int LDQ,
                                 float *Q2, int *K);
void QUARK_CORE_slaed4_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                            int n, const int *K,
                            float *D, const float *beta,
                            float **Q, const int *LDQ,
                            const float *DLAMBDA, const float *W, const int *INDX,
                            int start, int end,
                            PLASMA_sequence *sequence, PLASMA_request *request,
                            void *fakeQ, int flagfQ);
void QUARK_CORE_slaed3_compW_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                                  int n, const int *K,
                                  float **Q, const int *LDQ,
                                  const float *DLAMBDA, float *W,
                                  const int *INDX,
                                  int start, int end,
                                  void *fakeQ, int flagfQ,
                                  void *fakeW, int flagfW);

void QUARK_CORE_slaed3_reduceW(Quark *quark, Quark_Task_Flags *task_flags,
                               int n, int n1, const int *K, int l,
                               const float *Q, int LDQ,
                               const float *Wred, float *W);
void QUARK_CORE_slaed3_reduceW_p2(Quark *quark, Quark_Task_Flags *task_flags,
                                  int n, int n1, const int *K, int l,
                                  float **Q, const int *LDQ,
                                  const float *Wred, float *W);

void QUARK_CORE_slaed2_copydef(Quark *quark, Quark_Task_Flags *task_flags,
                               int n, int n1, const int *K, const int *ctot,
                               float *Q, int LDQ, const float *Q2,
                               int start, int end);
void QUARK_CORE_slaed3_computevectors(Quark *quark, Quark_Task_Flags *task_flags,
                                      int wsmode, int n, const int *K,
                      const int *il_nondef, const int *iu_nondef,
                                      float *Q, int LDQ, float *W, const int *INDXC,
                                      float **WSglobal, float **WSlocal,
                                      int start, int end );
void QUARK_CORE_slaed3_wscopy( Quark *quark, Quark_Task_Flags *task_flags,
                               const int *K, const int *il_nondef, const int *iu_nondef,
                               const float *Q, int LDQ, float **WORK,
                               int start, int end );
void QUARK_CORE_slaed3_updatevectors(Quark *quark, Quark_Task_Flags *task_flags,
                                     int oper, int wsmode, int n, int n1, int *K,
                     int *il_nondef, int *iu_nondef,
                                     float *D, float *Q, int LDQ, float *Q2,
                                     int *INDXQ, int *COLTYP, float **WORK,
                                     int start, int end, float **WORKDEP);
void QUARK_CORE_slaed3_pipelined(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, int *K, int *il_nondef, int *iu_nondef,
                                 float *D, float *Q, int LDQ, float *Q2,
                                 int *INDXC, int *INDXQ, int *COLTYP, float *W,
                                 int start, int end2);

void QUARK_CORE_sDC_fakedep(Quark *quark, Quark_Task_Flags *task_flags,
                            int nb_tasks, int nb, float *Q, int LDQ, float *W);
#endif

void QUARK_CORE_sswap(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, float *Q,
                      int LDQ, float *work,
                      int *perm, int begin, int end);
#ifdef COMPLEX
void QUARK_CORE_slag2c(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n,
                       const float *Q, int LDQ,
                       float *Z, int LDZ);
#endif
void QUARK_CORE_slaed3_freebigwork(Quark *quark, Quark_Task_Flags *task_flags,
                   int *K_bis, int largework, float **WORK);
void QUARK_CORE_slaset_identity(Quark *quark, Quark_Task_Flags *task_flags,
                int n, int start, int size,
                float *A);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_sasum_quark(Quark *quark);
void CORE_sasum_f1_quark(Quark *quark);
void CORE_sgeadd_quark(Quark *quark);
void CORE_sbrdalg1_quark(Quark *quark);
void CORE_sgelqt_quark(Quark *quark);
void CORE_sgemm_quark(Quark *quark);
void CORE_sgemm_tile_quark(Quark *quark);
void CORE_sgemv_quark(Quark *quark);
void CORE_sgemv_tile_quark(Quark *quark);
void CORE_sgeqp3_init_quark(Quark *quark);
void CORE_sgeqp3_larfg_quark(Quark *quark);
void CORE_sgeqp3_norms_quark(Quark *quark);
void CORE_sgeqp3_pivot_quark(Quark *quark);
void CORE_sgeqp3_tntpiv_quark(Quark *quark);
void CORE_sgeqp3_update_quark(Quark *quark);
void CORE_sgeqrt_quark(Quark *quark);
void CORE_sgessm_quark(Quark *quark);
void CORE_sgessq_quark(Quark *quark);
void CORE_sgessq_f1_quark(Quark *quark);
void CORE_sgetrf_quark(Quark *quark);
void CORE_sgetrf_incpiv_quark(Quark *quark);
void CORE_sgetrf_nopiv_quark(Quark* quark);
void CORE_sgetrf_reclap_quark(Quark *quark);
void CORE_sgetrf_rectil_quark(Quark* quark);
void CORE_sgetrip_quark(Quark *quark);
void CORE_sgetrip_f1_quark(Quark *quark);
void CORE_sgetrip_f2_quark(Quark *quark);
#ifdef COMPLEX
void CORE_ssymm_quark(Quark *quark);
void CORE_ssyrk_quark(Quark *quark);
void CORE_ssyr2k_quark(Quark *quark);
#endif
void CORE_ssygst_quark(Quark *quark);
void CORE_ssyrfb_quark(Quark *quark);
void CORE_shessq_quark(Quark *quark);
void CORE_shessq_f1_quark(Quark *quark);
void CORE_slacpy_quark(Quark *quark);
void CORE_slacpy_f1_quark(Quark *quark);
void CORE_slacpy_pivot_quark(Quark *quark);
void CORE_slatro_quark(Quark *quark);
void CORE_slatro_f1_quark(Quark *quark);
void CORE_slange_quark(Quark *quark);
void CORE_slange_f1_quark(Quark *quark);
#ifdef COMPLEX
void CORE_slansy_quark(Quark *quark);
void CORE_slansy_f1_quark(Quark *quark);
#endif
void CORE_slansy_quark(Quark *quark);
void CORE_slansy_f1_quark(Quark *quark);
void CORE_slaset_quark(Quark *quark);
void CORE_slaset2_quark(Quark *quark);
void CORE_slatro_quark(Quark *quark);
void CORE_slauum_quark(Quark *quark);
void CORE_spamm_quark(Quark *quark);
void CORE_splgsy_quark(Quark *quark);
void CORE_splgsy_quark(Quark *quark);
void CORE_splrnt_quark(Quark *quark);
void CORE_spltmg_quark(Quark *quark);
void CORE_splssq_quark(Quark *quark);
void CORE_spotrf_quark(Quark *quark);
void CORE_ssetvar_quark(Quark *quark);
void CORE_sshift_quark(Quark *quark);
void CORE_sshiftw_quark(Quark *quark);
void CORE_sssssm_quark(Quark *quark);
void CORE_ssymm_quark(Quark *quark);
void CORE_ssyrk_quark(Quark *quark);
void CORE_ssyr2k_quark(Quark *quark);
void CORE_ssyssq_quark(Quark *quark);
void CORE_ssyssq_f1_quark(Quark *quark);
void CORE_sswpab_quark(Quark *quark);
void CORE_sswptr_ontile_quark(Quark *quark);
void CORE_strdalg1_quark(Quark *quark);
void CORE_strmm_quark(Quark *quark);
void CORE_strsm_quark(Quark *quark);
void CORE_strtri_quark(Quark *quark);
void CORE_stslqt_quark(Quark *quark);
void CORE_stsmlq_quark(Quark *quark);
void CORE_stsmlq_sytra1_quark(Quark *quark);
void CORE_stsmlq_corner_quark(Quark *quark);
void CORE_stsmqr_quark(Quark *quark);
void CORE_stsmqr_sytra1_quark(Quark *quark);
void CORE_stsmqr_corner_quark(Quark *quark);
void CORE_stsqrt_quark(Quark *quark);
void CORE_ststrf_quark(Quark *quark);
void CORE_sttmqr_quark(Quark *quark);
void CORE_sttqrt_quark(Quark *quark);
void CORE_sttmlq_quark(Quark *quark);
void CORE_sttlqt_quark(Quark *quark);
void CORE_sormlq_quark(Quark *quark);
void CORE_sormqr_quark(Quark *quark);
void CORE_slaswp_quark(Quark* quark);
void CORE_slaswp_f2_quark(Quark* quark);
void CORE_slaswp_ontile_quark(Quark *quark);
void CORE_slaswp_ontile_f2_quark(Quark *quark);
void CORE_slaswpc_ontile_quark(Quark *quark);
void CORE_strmm_p2_quark(Quark* quark);
void CORE_sgemm_f2_quark(Quark* quark);
void CORE_sgemm_p2_quark(Quark* quark);
void CORE_sgemm_p2f1_quark(Quark* quark);
void CORE_sgemm_p3_quark(Quark* quark);

#endif /* defined(QUARK_H) */

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
