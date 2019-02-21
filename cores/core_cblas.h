/**
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
                 const parsec_complex32_t *A, int lda, float *work);
void CORE_cbrdalg1( PLASMA_enum uplo,
                    int n,
                    int nb,
                    parsec_complex32_t *A,
                    int lda,
                    parsec_complex32_t *VQ,
                    parsec_complex32_t *TAUQ,
                    parsec_complex32_t *VP,
                    parsec_complex32_t *TAUP,
                    int Vblksiz, int wantz,
                    int i, int sweepid, int m, int grsiz,
                    parsec_complex32_t *work);
int CORE_cgbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex32_t *V, parsec_complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_cgbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex32_t *V, parsec_complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_cgblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex32_t *V, parsec_complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_cgeadd(PLASMA_enum trans, int M, int N,
                      parsec_complex32_t alpha,
                const parsec_complex32_t *A, int LDA,
                      parsec_complex32_t beta,
                      parsec_complex32_t *B, int LDB);
int  CORE_cgelqt(int M, int N, int IB,
                 parsec_complex32_t *A, int LDA,
                 parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *TAU,
                 parsec_complex32_t *WORK);
void CORE_cgemm(PLASMA_enum transA, PLASMA_enum transB,
                int M, int N, int K,
                parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                                          const parsec_complex32_t *B, int LDB,
                parsec_complex32_t beta,        parsec_complex32_t *C, int LDC);
void CORE_cgemv(PLASMA_enum trans, int M, int N,
                parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                                          const parsec_complex32_t *x, int incx,
                parsec_complex32_t beta,        parsec_complex32_t *y, int incy);
void CORE_cgeqp3_init( int n, int *jpvt );
void CORE_cgeqp3_larfg( PLASMA_desc A, int ii, int jj, int i, int j,
                        parsec_complex32_t *tau, parsec_complex32_t *beta );
void CORE_cgeqp3_norms( PLASMA_desc A, int ioff, int joff, float *norms1, float *norms2 );
void CORE_cgeqp3_pivot( PLASMA_desc A, parsec_complex32_t *F, int ldf,
                        int jj, int k, int *jpvt,
                        float *norms1, float *norms2, int *info );
int  CORE_cgeqp3_tntpiv(int m, int n,
                        parsec_complex32_t *A, int lda,
                        int *IPIV, parsec_complex32_t *tau,
                        int *iwork);
void CORE_cgeqp3_update( const parsec_complex32_t *Ajj, int lda1,
                         parsec_complex32_t       *Ajk, int lda2,
                         const parsec_complex32_t *Fk,  int ldf,
                         int joff, int k, int koff, int nb,
                         float *norms1, float *norms2,
                         int *info );
int  CORE_cgeqrt(int M, int N, int IB,
                 parsec_complex32_t *A, int LDA,
                 parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *TAU, parsec_complex32_t *WORK);
int  CORE_cgessm(int M, int N, int K, int IB,
                 const int *IPIV,
                 const parsec_complex32_t *L, int LDL,
                 parsec_complex32_t *A, int LDA);
int  CORE_cgessq(int M, int N,
                 const parsec_complex32_t *A, int LDA,
                 float *scale, float *sumsq);
int  CORE_cgetf2_nopiv(int m, int n,
                      parsec_complex32_t *A, int lda);
int  CORE_cgetrf(int M, int N,
                 parsec_complex32_t *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_cgetrf_incpiv(int M, int N, int IB,
                        parsec_complex32_t *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_cgetrf_nopiv(int m, int n, int ib,
                      parsec_complex32_t *A, int lda);
int  CORE_cgetrf_reclap(CORE_cgetrf_data_t *data, int M, int N,
                        parsec_complex32_t *A, int LDA,
                        int *IPIV, int *info);
CORE_cgetrf_data_t *CORE_cgetrf_reclap_init(int nbthrd);
int  CORE_cgetrf_rectil(CORE_cgetrf_data_t *data, const PLASMA_desc A, int *IPIV, int *info);
CORE_cgetrf_data_t *CORE_cgetrf_rectil_init(int nbthrd);
void CORE_cgetrip(int m, int n, parsec_complex32_t *A,
                  parsec_complex32_t *work);
int CORE_chbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex32_t *V, parsec_complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_chblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex32_t *V, parsec_complex32_t *TAU,
                int st, int ed, int eltsize);
int CORE_chbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex32_t *V, parsec_complex32_t *TAU,
                int st, int ed, int eltsize);
void CORE_chbtype1cb(int N, int NB,
                     parsec_complex32_t *A, int LDA,
                     parsec_complex32_t *V, parsec_complex32_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     parsec_complex32_t *WORK);
void CORE_chbtype2cb(int N, int NB,
                     parsec_complex32_t *A, int LDA,
                     parsec_complex32_t *V, parsec_complex32_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     parsec_complex32_t *WORK);
void CORE_chbtype3cb(int N, int NB,
                     parsec_complex32_t *A, int LDA,
                     const parsec_complex32_t *V, const parsec_complex32_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     parsec_complex32_t *WORK);
void CORE_cgbtype1cb(PLASMA_enum uplo, int N, int NB,
                parsec_complex32_t *A, int LDA,
                parsec_complex32_t *VQ, parsec_complex32_t *TAUQ,
                parsec_complex32_t *VP, parsec_complex32_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                parsec_complex32_t *WORK);
void CORE_cgbtype2cb(PLASMA_enum uplo, int N, int NB,
                parsec_complex32_t *A, int LDA,
                parsec_complex32_t *VQ, parsec_complex32_t *TAUQ,
                parsec_complex32_t *VP, parsec_complex32_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                parsec_complex32_t *WORK);
void CORE_cgbtype3cb(PLASMA_enum uplo, int N, int NB,
                parsec_complex32_t *A, int LDA,
                parsec_complex32_t *VQ, parsec_complex32_t *TAUQ,
                parsec_complex32_t *VP, parsec_complex32_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                parsec_complex32_t *WORK);
void CORE_chegst(int itype, PLASMA_enum uplo, int N,
                 parsec_complex32_t *A, int LDA,
                 parsec_complex32_t *B, int LDB, int *INFO);
#ifdef COMPLEX
void CORE_chemm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                                          const parsec_complex32_t *B, int LDB,
                parsec_complex32_t beta,        parsec_complex32_t *C, int LDC);
void CORE_cherk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                float alpha, const parsec_complex32_t *A, int LDA,
                float beta,        parsec_complex32_t *C, int LDC);
void CORE_cher2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                                           const parsec_complex32_t *B, int LDB,
                 float beta,                    parsec_complex32_t *C, int LDC);
int  CORE_chessq(PLASMA_enum uplo, int N,
                 const parsec_complex32_t *A, int LDA,
                 float *scale, float *sumsq);
#endif
int  CORE_cherfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 const parsec_complex32_t *A,    int LDA,
                 const parsec_complex32_t *T,    int LDT,
                       parsec_complex32_t *C,    int LDC,
                       parsec_complex32_t *WORK, int LDWORK);
void CORE_clacpy(PLASMA_enum uplo, int M, int N,
                 const parsec_complex32_t *A, int LDA,
                       parsec_complex32_t *B, int LDB);
int CORE_clacpy_pivot( const PLASMA_desc descA,
                       PLASMA_enum direct,
                       int k1, int k2, const int *ipiv,
                       int *rankin, int *rankout,
                       parsec_complex32_t *A, int lda,
                       int init);
void CORE_clange(int norm, int M, int N,
                 const parsec_complex32_t *A, int LDA,
                 float *work, float *normA);
#ifdef COMPLEX
void CORE_clanhe(int norm, PLASMA_enum uplo, int N,
                 const parsec_complex32_t *A, int LDA,
                 float *work, float *normA);
#endif
void CORE_clansy(int norm, PLASMA_enum uplo, int N,
                 const parsec_complex32_t *A, int LDA,
                 float *work, float *normA);
void CORE_clantr(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N,
                 const parsec_complex32_t *A, int LDA,
                 float *work, float *normA);
int CORE_clarfb_gemm(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                     int M, int N, int K,
                     const parsec_complex32_t *V, int LDV,
                     const parsec_complex32_t *T, int LDT,
                           parsec_complex32_t *C, int LDC,
                           parsec_complex32_t *WORK, int LDWORK);
int CORE_clarfx2(PLASMA_enum side, int N,
                 parsec_complex32_t V,
                 parsec_complex32_t TAU,
                 parsec_complex32_t *C1, int LDC1,
                 parsec_complex32_t *C2, int LDC2);
int CORE_clarfx2c(PLASMA_enum uplo,
                  parsec_complex32_t V,
                  parsec_complex32_t TAU,
                  parsec_complex32_t *C1,
                  parsec_complex32_t *C2,
                  parsec_complex32_t *C3);
int CORE_clarfx2ce(PLASMA_enum uplo,
                   parsec_complex32_t *V,
                   parsec_complex32_t *TAU,
                   parsec_complex32_t *C1,
                   parsec_complex32_t *C2,
                   parsec_complex32_t *C3);
void CORE_clarfy(int N,
                 parsec_complex32_t *A, int LDA,
                 const parsec_complex32_t *V,
                 const parsec_complex32_t *TAU,
                 parsec_complex32_t *WORK);
int  CORE_clascal(PLASMA_enum uplo, int m, int n,
                  parsec_complex32_t alpha, parsec_complex32_t *A, int lda);
void CORE_claset(PLASMA_enum uplo, int n1, int n2,
                 parsec_complex32_t alpha, parsec_complex32_t beta,
                 parsec_complex32_t *tileA, int ldtilea);
void CORE_claset2(PLASMA_enum uplo, int n1, int n2, parsec_complex32_t alpha,
                  parsec_complex32_t *tileA, int ldtilea);
void CORE_claswp(int N, parsec_complex32_t *A, int LDA,
                 int I1,  int I2, const int *IPIV, int INC);
int  CORE_claswp_ontile( PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_claswpc_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_clatro(PLASMA_enum uplo, PLASMA_enum trans,
                 int M, int N,
                 const parsec_complex32_t *A, int LDA,
                       parsec_complex32_t *B, int LDB);
void CORE_clauum(PLASMA_enum uplo, int N, parsec_complex32_t *A, int LDA);
int CORE_cpamm(int op, PLASMA_enum side, PLASMA_enum storev,
               int M, int N, int K, int L,
               const parsec_complex32_t *A1, int LDA1,
                     parsec_complex32_t *A2, int LDA2,
               const parsec_complex32_t *V, int LDV,
                     parsec_complex32_t *W, int LDW);
int  CORE_cparfb(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                 int M1, int N1, int M2, int N2, int K, int L,
                       parsec_complex32_t *A1, int LDA1,
                       parsec_complex32_t *A2, int LDA2,
                 const parsec_complex32_t *V, int LDV,
                 const parsec_complex32_t *T, int LDT,
                       parsec_complex32_t *WORK, int LDWORK);
int CORE_cpemv(PLASMA_enum trans, PLASMA_enum storev,
               int M, int N, int L,
               parsec_complex32_t ALPHA,
               const parsec_complex32_t *A, int LDA,
               const parsec_complex32_t *X, int INCX,
               parsec_complex32_t BETA,
               parsec_complex32_t *Y, int INCY,
               parsec_complex32_t *WORK);
void CORE_cplghe(float bump, int m, int n, parsec_complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_cplgsy(parsec_complex32_t bump, int m, int n, parsec_complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_cplrnt(int m, int n, parsec_complex32_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
int  CORE_cpltmg(PLASMA_enum mtxtype, int m, int n, parsec_complex32_t *A, int lda,
                  int gM, int gN, int m0, int n0, unsigned long long int seed );
int  CORE_cpltmg_chebvand( int M, int N, parsec_complex32_t *A, int LDA,
                           int gN, int m0, int n0,
                           parsec_complex32_t *W );
int  CORE_cpltmg_circul( int M, int N, parsec_complex32_t *A, int LDA,
                         int gM, int m0, int n0,
                         const parsec_complex32_t *V );
void CORE_cpltmg_condexq( int M, int N, parsec_complex32_t *Q, int LDQ );
void CORE_cpltmg_fiedler(int m, int n,
                         const parsec_complex32_t *X, int incX,
                         const parsec_complex32_t *Y, int incY,
                         parsec_complex32_t *A, int lda);
int  CORE_cpltmg_hankel( PLASMA_enum uplo, int M, int N, parsec_complex32_t *A, int LDA,
                         int m0, int n0, int nb,
                         const parsec_complex32_t *V1,
                         const parsec_complex32_t *V2 );
void CORE_cpltmg_toeppd1( int gM, int m0, int M, parsec_complex32_t *W,
                          unsigned long long int seed );
void CORE_cpltmg_toeppd2( int M, int N, int K, int m0, int n0,
                          const parsec_complex32_t *W,
                          parsec_complex32_t *A, int LDA );
void CORE_cpotrf(PLASMA_enum uplo, int N, parsec_complex32_t *A, int LDA, int *INFO);
void CORE_csetvar(const parsec_complex32_t *alpha, parsec_complex32_t *x);
void CORE_cshift(int s, int m, int n, int L,
                 parsec_complex32_t *A);
void CORE_cshiftw(int s, int cl, int m, int n, int L,
                  parsec_complex32_t *A, parsec_complex32_t *W);
int  CORE_cssssm(int M1, int N1, int M2, int N2, int K, int IB,
                       parsec_complex32_t *A1, int LDA1,
                       parsec_complex32_t *A2, int LDA2,
                 const parsec_complex32_t *L1, int LDL1,
                 const parsec_complex32_t *L2, int LDL2,
                 const int *IPIV);
int CORE_cstedc(PLASMA_enum compz, int n,
                float *D, float *E,
                parsec_complex32_t *Z, int LDZ,
                parsec_complex32_t *WORK, int LWORK,
#ifdef COMPLEX
                float *RWORK, int LRWORK,
#endif
                int *IWORK, int LIWORK);
int CORE_csteqr(PLASMA_enum compz, int n,
                float *D, float *E,
                parsec_complex32_t *Z, int LDZ,
                float *WORK);
void CORE_csymm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                                          const parsec_complex32_t *B, int LDB,
                parsec_complex32_t beta,        parsec_complex32_t *C, int LDC);
void CORE_csyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                parsec_complex32_t beta,        parsec_complex32_t *C, int LDC);
void CORE_csyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                                           const parsec_complex32_t *B, int LDB,
                 parsec_complex32_t beta,        parsec_complex32_t *C, int LDC);
int  CORE_csyssq(PLASMA_enum uplo, int N,
                 const parsec_complex32_t *A, int LDA,
                 float *scale, float *sumsq);
void CORE_cswpab(int i, int n1, int n2,
                 parsec_complex32_t *A, parsec_complex32_t *work);
int  CORE_cswptr_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc,
                        const parsec_complex32_t *Akk, int ldak);
int CORE_ctradd(PLASMA_enum uplo, PLASMA_enum trans, int M, int N,
                      parsec_complex32_t alpha,
                const parsec_complex32_t *A, int LDA,
                      parsec_complex32_t beta,
                      parsec_complex32_t *B, int LDB);
void CORE_ctrasm(PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N, const parsec_complex32_t *A, int lda, float *work);
void CORE_ctrdalg1(int n,
                        int nb,
                        parsec_complex32_t *A,
                        int lda,
                        parsec_complex32_t *V,
                        parsec_complex32_t *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        parsec_complex32_t *work);
void CORE_ctrmm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                                                parsec_complex32_t *B, int LDB);
void CORE_ctrsm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                parsec_complex32_t alpha, const parsec_complex32_t *A, int LDA,
                                                parsec_complex32_t *B, int LDB);
int  CORE_ctrssq(PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                 const parsec_complex32_t *A, int LDA,
                 float *scale, float *sumsq);
void CORE_ctrtri(PLASMA_enum uplo, PLASMA_enum diag, int N,
                 parsec_complex32_t *A, int LDA, int *info);
int  CORE_ctslqt(int M, int N, int IB,
                 parsec_complex32_t *A1, int LDA1,
                 parsec_complex32_t *A2, int LDA2,
                 parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *TAU, parsec_complex32_t *WORK);
int  CORE_ctsmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 parsec_complex32_t *A1, int LDA1,
                 parsec_complex32_t *A2, int LDA2,
                 const parsec_complex32_t *V, int LDV,
                 const parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *WORK, int LDWORK);
int CORE_ctsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        parsec_complex32_t *A1, int lda1,
                        parsec_complex32_t *A2, int lda2,
                        parsec_complex32_t *A3, int lda3,
                        const parsec_complex32_t *V, int ldv,
                        const parsec_complex32_t *T, int ldt,
                        parsec_complex32_t *WORK, int ldwork);
int CORE_ctsmlq_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        parsec_complex32_t *A1, int lda1,
                        parsec_complex32_t *A2, int lda2,
                        const parsec_complex32_t *V, int ldv,
                        const parsec_complex32_t *T, int ldt,
                        parsec_complex32_t *WORK, int ldwork);
int  CORE_ctsmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 parsec_complex32_t *A1, int LDA1,
                 parsec_complex32_t *A2, int LDA2,
                 const parsec_complex32_t *V, int LDV,
                 const parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *WORK, int LDWORK);
int CORE_ctsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        parsec_complex32_t *A1, int lda1,
                        parsec_complex32_t *A2, int lda2,
                        parsec_complex32_t *A3, int lda3,
                        const parsec_complex32_t *V, int ldv,
                        const parsec_complex32_t *T, int ldt,
                        parsec_complex32_t *WORK, int ldwork);
int CORE_ctsmqr_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        parsec_complex32_t *A1, int lda1,
                        parsec_complex32_t *A2, int lda2,
                        const parsec_complex32_t *V, int ldv,
                        const parsec_complex32_t *T, int ldt,
                        parsec_complex32_t *WORK, int ldwork);
int  CORE_ctsqrt(int M, int N, int IB,
                 parsec_complex32_t *A1, int LDA1,
                 parsec_complex32_t *A2, int LDA2,
                 parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *TAU, parsec_complex32_t *WORK);
int  CORE_ctstrf(int M, int N, int IB, int NB,
                 parsec_complex32_t *U, int LDU,
                 parsec_complex32_t *A, int LDA,
                 parsec_complex32_t *L, int LDL,
                 int *IPIV, parsec_complex32_t *WORK,
                 int LDWORK, int *INFO);
int  CORE_cttmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 parsec_complex32_t *A1, int LDA1,
                 parsec_complex32_t *A2, int LDA2,
                 const parsec_complex32_t *V, int LDV,
                 const parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *WORK, int LDWORK);
int  CORE_cttqrt(int M, int N, int IB,
                 parsec_complex32_t *A1, int LDA1,
                 parsec_complex32_t *A2, int LDA2,
                 parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *TAU,
                 parsec_complex32_t *WORK);
int  CORE_cttmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 parsec_complex32_t *A1, int LDA1,
                 parsec_complex32_t *A2, int LDA2,
                 const parsec_complex32_t *V, int LDV,
                 const parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *WORK, int LDWORK);
int  CORE_cttlqt(int M, int N, int IB,
                 parsec_complex32_t *A1, int LDA1,
                 parsec_complex32_t *A2, int LDA2,
                 parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *TAU,
                 parsec_complex32_t *WORK);
int  CORE_cunmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const parsec_complex32_t *V, int LDV,
                 const parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *C, int LDC,
                 parsec_complex32_t *WORK, int LDWORK);
int  CORE_cunmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const parsec_complex32_t *V, int LDV,
                 const parsec_complex32_t *T, int LDT,
                 parsec_complex32_t *C, int LDC,
                 parsec_complex32_t *WORK, int LDWORK);

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
void CORE_cswap(int m, int n, parsec_complex32_t *Q, int ldq,
                const parsec_complex32_t *work, const int *perm,
                int start, int end);
int CORE_clascl(PLASMA_enum type, int kl, int ku, float cfrom, float cto,
                int m, int n, parsec_complex32_t *A, int lda);
#ifdef COMPLEX
int CORE_slag2c(int m, int n, const float *Q, int LDQ,
                 parsec_complex32_t *Z, int LDZ);
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
void QUARK_CORE_scasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                       const parsec_complex32_t *A, int lda, int szeA,
                       float *work, int szeW);
void QUARK_CORE_scasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                          const parsec_complex32_t *A, int lda, int szeA,
                          float *work, int szeW,
                          float *fake, int szeF);
void QUARK_CORE_cgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum trans, int m, int n, int nb,
                       parsec_complex32_t  alpha,
                       const parsec_complex32_t *A, int lda,
                       parsec_complex32_t  beta,
                       parsec_complex32_t *B, int ldb);
void QUARK_CORE_cbrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo,
                        int n, int nb,
                        parsec_complex32_t *A,
                        int lda,
                        parsec_complex32_t *VQ,
                        parsec_complex32_t *TAUQ,
                        parsec_complex32_t *VP,
                        parsec_complex32_t *TAUP,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_cgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex32_t *A, int lda,
                       parsec_complex32_t *T, int ldt);
void QUARK_CORE_cgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int k, int nb,
                      parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                      const parsec_complex32_t *B, int ldb,
                      parsec_complex32_t beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_cgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum transA, PLASMA_enum transB,
                        int m, int n, int k, int nb,
                        parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                        const parsec_complex32_t *B, int ldb,
                        parsec_complex32_t beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_cgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                         const parsec_complex32_t *B, int ldb,
                         parsec_complex32_t beta, parsec_complex32_t *C, int ldc,
                         parsec_complex32_t *fake1, int szefake1, int flag1,
                         parsec_complex32_t *fake2, int szefake2, int flag2);
void QUARK_CORE_cgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                         const parsec_complex32_t **B, int ldb,
                         parsec_complex32_t beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_cgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                           const parsec_complex32_t **B, int ldb,
                           parsec_complex32_t beta, parsec_complex32_t *C, int ldc,
                           parsec_complex32_t *fake1, int szefake1, int flag1);
void QUARK_CORE_cgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                         const parsec_complex32_t *B, int ldb,
                         parsec_complex32_t beta, parsec_complex32_t **C, int ldc);
void QUARK_CORE_cgemm_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           const parsec_complex32_t *alpha, const parsec_complex32_t *A, int lda,
                                                            const parsec_complex32_t *B, int ldb,
                           const parsec_complex32_t *beta,        parsec_complex32_t *C, int ldc,
                           const parsec_complex32_t *Alock,
                           const parsec_complex32_t *Block,
                           const parsec_complex32_t *Clock);
void QUARK_CORE_cgemv(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum trans, int m, int n,
                      parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                                                const parsec_complex32_t *x, int incx,
                      parsec_complex32_t beta,        parsec_complex32_t *y, int incy);
void QUARK_CORE_cgemv_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum trans,
                           int m, int n,
                           const parsec_complex32_t *alpha, const parsec_complex32_t *A, int lda,
                                                            const parsec_complex32_t *x, int incx,
                           const parsec_complex32_t *beta,        parsec_complex32_t *y, int incy,
                           const parsec_complex32_t *Alock,
                           const parsec_complex32_t *xlock,
                           const parsec_complex32_t *ylock);
void QUARK_CORE_cgeqp3_init( Quark *quark, Quark_Task_Flags *task_flags,
                             int n, int *jpvt );
void QUARK_CORE_cgeqp3_larfg(Quark *quark, Quark_Task_Flags *task_flags,
                             PLASMA_desc A, int ii, int jj, int i, int j,
                             parsec_complex32_t *tau, parsec_complex32_t *beta );
void QUARK_CORE_cgeqp3_norms( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, int ioff, int joff, float *norms1, float *norms2 );
void QUARK_CORE_cgeqp3_pivot( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A,
                              parsec_complex32_t *F, int ldf,
                              int jj, int k, int *jpvt,
                              float *norms1, float *norms2, int *info );
void QUARK_CORE_cgeqp3_tntpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              parsec_complex32_t *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_cgeqp3_update( Quark *quark, Quark_Task_Flags *task_flags,
                               parsec_complex32_t *Ajj, int lda1,
                               parsec_complex32_t *Ajk, int lda2,
                               parsec_complex32_t *Fk,  int ldf,
                               int joff, int k, int koff, int nb,
                               float *norms1, float *norms2, int *info );
void QUARK_CORE_cgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex32_t *A, int lda,
                       parsec_complex32_t *T, int ldt);
void QUARK_CORE_cgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       const int *IPIV,
                       const parsec_complex32_t *L, int ldl,
                       parsec_complex32_t *A, int lda);
void QUARK_CORE_cgessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, const parsec_complex32_t *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF );
void QUARK_CORE_cgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       parsec_complex32_t *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_cgetrf_incpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int ib, int nb,
                              parsec_complex32_t *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_cgetrf_nopiv(Quark *quark, Quark_Task_Flags *task_flags,
                             int m, int n, int ib, int nb,
                             parsec_complex32_t *A, int lda,
                             PLASMA_sequence *sequence, PLASMA_request *request,
                             int iinfo);
void QUARK_CORE_cgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_cgetrf_data_t *data, int m, int n, int nb,
                              parsec_complex32_t *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_cgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_cgetrf_data_t *data,
                              PLASMA_desc A, parsec_complex32_t *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_cgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, parsec_complex32_t *A, int szeA);
void QUARK_CORE_cgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, parsec_complex32_t *A, int szeA,
                           parsec_complex32_t *fake, int szeF, int paramF);
void QUARK_CORE_cgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, parsec_complex32_t *A, int szeA,
                           parsec_complex32_t *fake1, int szeF1, int paramF1,
                           parsec_complex32_t *fake2, int szeF2, int paramF2);
void QUARK_CORE_chemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                      const parsec_complex32_t *B, int ldb,
                      parsec_complex32_t beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_chegst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, PLASMA_enum uplo, int N,
                       parsec_complex32_t *A, int LDA,
                       parsec_complex32_t *B, int LDB,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_cherk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      float alpha, const parsec_complex32_t *A, int lda,
                      float beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_cher2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                       const parsec_complex32_t *B, int LDB,
                       float beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_cherfb(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo,
                       int n, int k, int ib, int nb,
                       const parsec_complex32_t *A, int lda,
                       const parsec_complex32_t *T, int ldt,
                       parsec_complex32_t *C, int ldc);
void QUARK_CORE_chessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const parsec_complex32_t *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF );
void QUARK_CORE_clacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int mb,
                       const parsec_complex32_t *A, int lda,
                       parsec_complex32_t *B, int ldb);
void QUARK_CORE_clacpy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, int m, int n, int nb,
                          const parsec_complex32_t *A, int lda,
                          parsec_complex32_t *B, int ldb,
                          parsec_complex32_t *fake1, int szefake1, int flag1);
void QUARK_CORE_clacpy_pivot(Quark *quark, Quark_Task_Flags *task_flags,
                             const PLASMA_desc descA,
                             PLASMA_enum direct,
                             int k1, int k2, const int *ipiv,
                             int *rankin, int *rankout,
                             parsec_complex32_t *A, int lda,
                             int pos, int init);
void QUARK_CORE_clange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       const parsec_complex32_t *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_clange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int M, int N,
                          const parsec_complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
#ifdef COMPLEX
void QUARK_CORE_clanhe(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const parsec_complex32_t *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_clanhe_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, PLASMA_enum uplo, int N,
                          const parsec_complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
#endif
void QUARK_CORE_clansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const parsec_complex32_t *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_clansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, PLASMA_enum uplo, int N,
                          const parsec_complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
void QUARK_CORE_clantr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                       const parsec_complex32_t *A, int LDA, int szeA,
                       int szeW, float *result);
void QUARK_CORE_clantr_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                          const parsec_complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF);
void QUARK_CORE_clascal(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int m, int n, int nb,
                        parsec_complex32_t alpha, parsec_complex32_t *A, int lda);
void QUARK_CORE_claset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n1, int n2, parsec_complex32_t alpha,
                       parsec_complex32_t beta, parsec_complex32_t *tileA, int ldtilea);
void QUARK_CORE_claset2(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int n1, int n2, parsec_complex32_t alpha,
                        parsec_complex32_t *tileA, int ldtilea);
void QUARK_CORE_claswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, parsec_complex32_t *A, int lda,
                       int i1,  int i2, const int *ipiv, int inc);
void QUARK_CORE_claswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, parsec_complex32_t *A, int lda,
                          int i1,  int i2, const int *ipiv, int inc,
                          parsec_complex32_t *fake1, int szefake1, int flag1,
                          parsec_complex32_t *fake2, int szefake2, int flag2);
void QUARK_CORE_claswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, parsec_complex32_t *A,
                              int i1,  int i2, const int *ipiv, int inc, parsec_complex32_t *fakepanel);
void QUARK_CORE_claswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, parsec_complex32_t *A,
                                 int i1,  int i2, const int *ipiv, int inc,
                                 parsec_complex32_t *fake1, int szefake1, int flag1,
                                 parsec_complex32_t *fake2, int szefake2, int flag2);
void QUARK_CORE_claswpc_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_desc descA, parsec_complex32_t *A,
                               int i1,  int i2, const int *ipiv, int inc, parsec_complex32_t *fakepanel);
void QUARK_CORE_clatro(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int mb,
                       const parsec_complex32_t *A, int lda,
                       parsec_complex32_t *B, int ldb);
void QUARK_CORE_clatro_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int mb,
                          const parsec_complex32_t *A, int lda,
                                parsec_complex32_t *B, int ldb,
                          parsec_complex32_t *fake1, int szefake1, int flag1);
void QUARK_CORE_clauum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       parsec_complex32_t *A, int lda);
void QUARK_CORE_cplghe(Quark *quark, Quark_Task_Flags *task_flags,
                       float bump, int m, int n, parsec_complex32_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_cplgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       parsec_complex32_t bump, int m, int n, parsec_complex32_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_cplrnt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, parsec_complex32_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_cpltmg(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum mtxtype, int m, int n, parsec_complex32_t *A, int lda,
                        int gM, int gN, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_cpltmg_chebvand( Quark *quark, Quark_Task_Flags *task_flags,
                                 int M, int N, parsec_complex32_t *A, int LDA,
                                 int gN, int m0, int n0,
                                 parsec_complex32_t *W );
void QUARK_CORE_cpltmg_circul( Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, parsec_complex32_t *A, int LDA,
                               int gM, int m0, int n0,
                               const parsec_complex32_t *W );
void QUARK_CORE_cpltmg_fiedler(Quark *quark, Quark_Task_Flags *task_flags,
                               int m, int n,
                               const parsec_complex32_t *X, int incX,
                               const parsec_complex32_t *Y, int incY,
                               parsec_complex32_t *A, int lda);
void QUARK_CORE_cpltmg_hankel( Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_enum uplo, int M, int N, parsec_complex32_t *A, int LDA,
                               int m0, int n0, int nb,
                               const parsec_complex32_t *V1,
                               const parsec_complex32_t *V2);
void QUARK_CORE_cpltmg_toeppd1(Quark *quark, Quark_Task_Flags *task_flags,
                               int gM, int m0, int M,
                               parsec_complex32_t *W,
                               unsigned long long int seed);
void QUARK_CORE_cpltmg_toeppd2(Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, int K, int m0, int n0,
                               const parsec_complex32_t *W,
                               parsec_complex32_t *A, int LDA );
void QUARK_CORE_cpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       parsec_complex32_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_csetvar(Quark *quark, Quark_Task_Flags *task_flags,
                        const parsec_complex32_t *alpha, parsec_complex32_t *x,
                        parsec_complex32_t *Alock);
void QUARK_CORE_cshift( Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int m, int n, int L,
                        parsec_complex32_t *A);
void QUARK_CORE_cshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L,
                        parsec_complex32_t *A, parsec_complex32_t *W);
void QUARK_CORE_cssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       const parsec_complex32_t *L1, int ldl1,
                       const parsec_complex32_t *L2, int ldl2,
                       const int *IPIV);
void QUARK_CORE_cstedc(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       float *D, float *E,
                       parsec_complex32_t *Z, int ldz);
void QUARK_CORE_cstedc_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum compz, int n,
                          float *D, float *E,
                          parsec_complex32_t *Z, int ldz,
                          void *fake1, int szefake1, int flag1,
                          void *fake2, int szefake2, int flag2);
void QUARK_CORE_csteqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       float *D, float *E,
                       parsec_complex32_t *Z, int ldz);
void QUARK_CORE_csymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                      const parsec_complex32_t *B, int ldb,
                      parsec_complex32_t beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_csyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                      parsec_complex32_t beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_csyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                       const parsec_complex32_t *B, int LDB,
                       parsec_complex32_t beta, parsec_complex32_t *C, int ldc);
void QUARK_CORE_csyssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const parsec_complex32_t *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF );
void QUARK_CORE_cswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       parsec_complex32_t *A, int szeA);
void QUARK_CORE_cswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, parsec_complex32_t *Aij,
                              int i1,  int i2, const int *ipiv, int inc,
                              const parsec_complex32_t *Akk, int ldak);
void QUARK_CORE_ctradd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       parsec_complex32_t  alpha,
                       const parsec_complex32_t *A, int lda,
                       parsec_complex32_t  beta,
                       parsec_complex32_t *B, int ldb);
void QUARK_CORE_ctrasm(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int m, int n,
                       const parsec_complex32_t *A, int lda, int szeA,
                       float *work, int szeW);
void QUARK_CORE_ctrasm_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int m, int n,
                          const parsec_complex32_t *A, int lda, int szeA,
                          float *work, int szeW,
                          float *fake, int szeF);
void QUARK_CORE_ctrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                        int n,
                        int nb,
                        parsec_complex32_t *A,
                        int lda,
                        parsec_complex32_t *V,
                        parsec_complex32_t *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_ctrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                      parsec_complex32_t *B, int ldb);
void QUARK_CORE_ctrmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                         int m, int n, int nb,
                         parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                         parsec_complex32_t **B, int ldb);
void QUARK_CORE_ctrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      parsec_complex32_t alpha, const parsec_complex32_t *A, int lda,
                      parsec_complex32_t *B, int ldb);
void QUARK_CORE_ctrssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, PLASMA_enum diag,
                           int m, int n, const parsec_complex32_t *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF );
void QUARK_CORE_ctrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum diag, int n, int nb,
                       parsec_complex32_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_ctslqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       parsec_complex32_t *T, int ldt);
void QUARK_CORE_ctsmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       const parsec_complex32_t *V, int ldv,
                       const parsec_complex32_t *T, int ldt);
void QUARK_CORE_ctsmlq_hetra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              parsec_complex32_t *A1, int lda1,
                              parsec_complex32_t *A2, int lda2,
                              const parsec_complex32_t *V, int ldv,
                              const parsec_complex32_t *T, int ldt);
void QUARK_CORE_ctsmlq_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              parsec_complex32_t *A1, int lda1,
                              parsec_complex32_t *A2, int lda2,
                              parsec_complex32_t *A3, int lda3,
                              const parsec_complex32_t *V, int ldv,
                              const parsec_complex32_t *T, int ldt);
void QUARK_CORE_ctsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       const parsec_complex32_t *V, int ldv,
                       const parsec_complex32_t *T, int ldt);
void QUARK_CORE_ctsmqr_hetra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              parsec_complex32_t *A1, int lda1,
                              parsec_complex32_t *A2, int lda2,
                              const parsec_complex32_t *V, int ldv,
                              const parsec_complex32_t *T, int ldt);
void QUARK_CORE_ctsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              parsec_complex32_t *A1, int lda1,
                              parsec_complex32_t *A2, int lda2,
                              parsec_complex32_t *A3, int lda3,
                              const parsec_complex32_t *V, int ldv,
                              const parsec_complex32_t *T, int ldt);
void QUARK_CORE_ctsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       parsec_complex32_t *T, int ldt);
void QUARK_CORE_ctstrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex32_t *U, int ldu,
                       parsec_complex32_t *A, int lda,
                       parsec_complex32_t *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_cttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       const parsec_complex32_t *V, int ldv,
                       const parsec_complex32_t *T, int ldt);
void QUARK_CORE_cttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       parsec_complex32_t *T, int ldt);
void QUARK_CORE_cttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       const parsec_complex32_t *V, int ldv,
                       const parsec_complex32_t *T, int ldt);
void QUARK_CORE_cttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       parsec_complex32_t *T, int ldt);
void QUARK_CORE_cpamm(Quark *quark, Quark_Task_Flags *task_flags,
                       int op, PLASMA_enum side, PLASMA_enum storev,
                       int m, int n, int k, int l,
                       const parsec_complex32_t *A1, int lda1,
                       parsec_complex32_t *A2, int lda2,
                       const parsec_complex32_t *V, int ldv,
                       parsec_complex32_t *W, int ldw);
void QUARK_CORE_cplssq( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, const float *A, float *result );
void QUARK_CORE_cunmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int ib,  int nb, int k,
                       const parsec_complex32_t *A, int lda,
                       const parsec_complex32_t *T, int ldt,
                       parsec_complex32_t *C, int ldc);
void QUARK_CORE_cunmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const parsec_complex32_t *A, int lda,
                       const parsec_complex32_t *T, int ldt,
                       parsec_complex32_t *C, int ldc);


void QUARK_CORE_clascl(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum type, int kl, int ku, float cfrom, float cto,
                       int m, int n, parsec_complex32_t *A, int lda);
void QUARK_CORE_clascl_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                            PLASMA_enum type, int kl, int ku, float *cfrom, float *cto,
                            int m, int n, parsec_complex32_t *A, int lda,
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

void QUARK_CORE_cswap(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, parsec_complex32_t *Q,
                      int LDQ, parsec_complex32_t *work,
                      int *perm, int begin, int end);
#ifdef COMPLEX
void QUARK_CORE_slag2c(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n,
                       const float *Q, int LDQ,
                       parsec_complex32_t *Z, int LDZ);
#endif
void QUARK_CORE_slaed3_freebigwork(Quark *quark, Quark_Task_Flags *task_flags,
                   int *K_bis, int largework, float **WORK);
void QUARK_CORE_claset_identity(Quark *quark, Quark_Task_Flags *task_flags,
                int n, int start, int size,
                parsec_complex32_t *A);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_scasum_quark(Quark *quark);
void CORE_scasum_f1_quark(Quark *quark);
void CORE_cgeadd_quark(Quark *quark);
void CORE_cbrdalg1_quark(Quark *quark);
void CORE_cgelqt_quark(Quark *quark);
void CORE_cgemm_quark(Quark *quark);
void CORE_cgemm_tile_quark(Quark *quark);
void CORE_cgemv_quark(Quark *quark);
void CORE_cgemv_tile_quark(Quark *quark);
void CORE_cgeqp3_init_quark(Quark *quark);
void CORE_cgeqp3_larfg_quark(Quark *quark);
void CORE_cgeqp3_norms_quark(Quark *quark);
void CORE_cgeqp3_pivot_quark(Quark *quark);
void CORE_cgeqp3_tntpiv_quark(Quark *quark);
void CORE_cgeqp3_update_quark(Quark *quark);
void CORE_cgeqrt_quark(Quark *quark);
void CORE_cgessm_quark(Quark *quark);
void CORE_cgessq_quark(Quark *quark);
void CORE_cgessq_f1_quark(Quark *quark);
void CORE_cgetrf_quark(Quark *quark);
void CORE_cgetrf_incpiv_quark(Quark *quark);
void CORE_cgetrf_nopiv_quark(Quark* quark);
void CORE_cgetrf_reclap_quark(Quark *quark);
void CORE_cgetrf_rectil_quark(Quark* quark);
void CORE_cgetrip_quark(Quark *quark);
void CORE_cgetrip_f1_quark(Quark *quark);
void CORE_cgetrip_f2_quark(Quark *quark);
#ifdef COMPLEX
void CORE_chemm_quark(Quark *quark);
void CORE_cherk_quark(Quark *quark);
void CORE_cher2k_quark(Quark *quark);
#endif
void CORE_chegst_quark(Quark *quark);
void CORE_cherfb_quark(Quark *quark);
void CORE_chessq_quark(Quark *quark);
void CORE_chessq_f1_quark(Quark *quark);
void CORE_clacpy_quark(Quark *quark);
void CORE_clacpy_f1_quark(Quark *quark);
void CORE_clacpy_pivot_quark(Quark *quark);
void CORE_clatro_quark(Quark *quark);
void CORE_clatro_f1_quark(Quark *quark);
void CORE_clange_quark(Quark *quark);
void CORE_clange_f1_quark(Quark *quark);
#ifdef COMPLEX
void CORE_clanhe_quark(Quark *quark);
void CORE_clanhe_f1_quark(Quark *quark);
#endif
void CORE_clansy_quark(Quark *quark);
void CORE_clansy_f1_quark(Quark *quark);
void CORE_claset_quark(Quark *quark);
void CORE_claset2_quark(Quark *quark);
void CORE_clatro_quark(Quark *quark);
void CORE_clauum_quark(Quark *quark);
void CORE_cpamm_quark(Quark *quark);
void CORE_cplghe_quark(Quark *quark);
void CORE_cplgsy_quark(Quark *quark);
void CORE_cplrnt_quark(Quark *quark);
void CORE_cpltmg_quark(Quark *quark);
void CORE_cplssq_quark(Quark *quark);
void CORE_cpotrf_quark(Quark *quark);
void CORE_csetvar_quark(Quark *quark);
void CORE_cshift_quark(Quark *quark);
void CORE_cshiftw_quark(Quark *quark);
void CORE_cssssm_quark(Quark *quark);
void CORE_csymm_quark(Quark *quark);
void CORE_csyrk_quark(Quark *quark);
void CORE_csyr2k_quark(Quark *quark);
void CORE_csyssq_quark(Quark *quark);
void CORE_csyssq_f1_quark(Quark *quark);
void CORE_cswpab_quark(Quark *quark);
void CORE_cswptr_ontile_quark(Quark *quark);
void CORE_ctrdalg1_quark(Quark *quark);
void CORE_ctrmm_quark(Quark *quark);
void CORE_ctrsm_quark(Quark *quark);
void CORE_ctrtri_quark(Quark *quark);
void CORE_ctslqt_quark(Quark *quark);
void CORE_ctsmlq_quark(Quark *quark);
void CORE_ctsmlq_hetra1_quark(Quark *quark);
void CORE_ctsmlq_corner_quark(Quark *quark);
void CORE_ctsmqr_quark(Quark *quark);
void CORE_ctsmqr_hetra1_quark(Quark *quark);
void CORE_ctsmqr_corner_quark(Quark *quark);
void CORE_ctsqrt_quark(Quark *quark);
void CORE_ctstrf_quark(Quark *quark);
void CORE_cttmqr_quark(Quark *quark);
void CORE_cttqrt_quark(Quark *quark);
void CORE_cttmlq_quark(Quark *quark);
void CORE_cttlqt_quark(Quark *quark);
void CORE_cunmlq_quark(Quark *quark);
void CORE_cunmqr_quark(Quark *quark);
void CORE_claswp_quark(Quark* quark);
void CORE_claswp_f2_quark(Quark* quark);
void CORE_claswp_ontile_quark(Quark *quark);
void CORE_claswp_ontile_f2_quark(Quark *quark);
void CORE_claswpc_ontile_quark(Quark *quark);
void CORE_ctrmm_p2_quark(Quark* quark);
void CORE_cgemm_f2_quark(Quark* quark);
void CORE_cgemm_p2_quark(Quark* quark);
void CORE_cgemm_p2f1_quark(Quark* quark);
void CORE_cgemm_p3_quark(Quark* quark);

#endif /* defined(QUARK_H) */

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
