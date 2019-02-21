/**
 *
 * @file core_dblas.h
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
 * @generated d Fri Apr  1 11:02:20 2016
 *
 **/
#ifndef _PLASMA_CORE_DBLAS_H_
#define _PLASMA_CORE_DBLAS_H_

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

struct CORE_dgetrf_data_s;
typedef struct CORE_dgetrf_data_s CORE_dgetrf_data_t;

/** ****************************************************************************
 *  Declarations of serial kernels - alphabetical order
 **/
void CORE_dasum(int storev, PLASMA_enum uplo, int M, int N,
                 const double *A, int lda, double *work);
void CORE_dbrdalg1( PLASMA_enum uplo,
                    int n,
                    int nb,
                    double *A,
                    int lda,
                    double *VQ,
                    double *TAUQ,
                    double *VP,
                    double *TAUP,
                    int Vblksiz, int wantz,
                    int i, int sweepid, int m, int grsiz,
                    double *work);
int CORE_dgbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, double *V, double *TAU,
                int st, int ed, int eltsize);
int CORE_dgbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, double *V, double *TAU,
                int st, int ed, int eltsize);
int CORE_dgblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, double *V, double *TAU,
                int st, int ed, int eltsize);
int CORE_dgeadd(PLASMA_enum trans, int M, int N,
                      double alpha,
                const double *A, int LDA,
                      double beta,
                      double *B, int LDB);
int  CORE_dgelqt(int M, int N, int IB,
                 double *A, int LDA,
                 double *T, int LDT,
                 double *TAU,
                 double *WORK);
void CORE_dgemm(PLASMA_enum transA, PLASMA_enum transB,
                int M, int N, int K,
                double alpha, const double *A, int LDA,
                                          const double *B, int LDB,
                double beta,        double *C, int LDC);
void CORE_dgemv(PLASMA_enum trans, int M, int N,
                double alpha, const double *A, int LDA,
                                          const double *x, int incx,
                double beta,        double *y, int incy);
void CORE_dgeqp3_init( int n, int *jpvt );
void CORE_dgeqp3_larfg( PLASMA_desc A, int ii, int jj, int i, int j,
                        double *tau, double *beta );
void CORE_dgeqp3_norms( PLASMA_desc A, int ioff, int joff, double *norms1, double *norms2 );
void CORE_dgeqp3_pivot( PLASMA_desc A, double *F, int ldf,
                        int jj, int k, int *jpvt,
                        double *norms1, double *norms2, int *info );
int  CORE_dgeqp3_tntpiv(int m, int n,
                        double *A, int lda,
                        int *IPIV, double *tau,
                        int *iwork);
void CORE_dgeqp3_update( const double *Ajj, int lda1,
                         double       *Ajk, int lda2,
                         const double *Fk,  int ldf,
                         int joff, int k, int koff, int nb,
                         double *norms1, double *norms2,
                         int *info );
int  CORE_dgeqrt(int M, int N, int IB,
                 double *A, int LDA,
                 double *T, int LDT,
                 double *TAU, double *WORK);
int  CORE_dgessm(int M, int N, int K, int IB,
                 const int *IPIV,
                 const double *L, int LDL,
                 double *A, int LDA);
int  CORE_dgessq(int M, int N,
                 const double *A, int LDA,
                 double *scale, double *sumsq);
int  CORE_dgetf2_nopiv(int m, int n,
                      double *A, int lda);
int  CORE_dgetrf(int M, int N,
                 double *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_dgetrf_incpiv(int M, int N, int IB,
                        double *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_dgetrf_nopiv(int m, int n, int ib,
                      double *A, int lda);
int  CORE_dgetrf_reclap(CORE_dgetrf_data_t *data, int M, int N,
                        double *A, int LDA,
                        int *IPIV, int *info);
CORE_dgetrf_data_t *CORE_dgetrf_reclap_init(int nbthrd);
int  CORE_dgetrf_rectil(CORE_dgetrf_data_t *data, const PLASMA_desc A, int *IPIV, int *info);
CORE_dgetrf_data_t *CORE_dgetrf_rectil_init(int nbthrd);
void CORE_dgetrip(int m, int n, double *A,
                  double *work);
int CORE_dhbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, double *V, double *TAU,
                int st, int ed, int eltsize);
int CORE_dhblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, double *V, double *TAU,
                int st, int ed, int eltsize);
int CORE_dhbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, double *V, double *TAU,
                int st, int ed, int eltsize);
void CORE_dsbtype1cb(int N, int NB,
                     double *A, int LDA,
                     double *V, double *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     double *WORK);
void CORE_dsbtype2cb(int N, int NB,
                     double *A, int LDA,
                     double *V, double *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     double *WORK);
void CORE_dsbtype3cb(int N, int NB,
                     double *A, int LDA,
                     const double *V, const double *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     double *WORK);
void CORE_dgbtype1cb(PLASMA_enum uplo, int N, int NB,
                double *A, int LDA,
                double *VQ, double *TAUQ,
                double *VP, double *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                double *WORK);
void CORE_dgbtype2cb(PLASMA_enum uplo, int N, int NB,
                double *A, int LDA,
                double *VQ, double *TAUQ,
                double *VP, double *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                double *WORK);
void CORE_dgbtype3cb(PLASMA_enum uplo, int N, int NB,
                double *A, int LDA,
                double *VQ, double *TAUQ,
                double *VP, double *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                double *WORK);
void CORE_dsygst(int itype, PLASMA_enum uplo, int N,
                 double *A, int LDA,
                 double *B, int LDB, int *INFO);
#ifdef COMPLEX
void CORE_dsymm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                double alpha, const double *A, int LDA,
                                          const double *B, int LDB,
                double beta,        double *C, int LDC);
void CORE_dsyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                double alpha, const double *A, int LDA,
                double beta,        double *C, int LDC);
void CORE_dsyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 double alpha, const double *A, int LDA,
                                           const double *B, int LDB,
                 double beta,                    double *C, int LDC);
int  CORE_dhessq(PLASMA_enum uplo, int N,
                 const double *A, int LDA,
                 double *scale, double *sumsq);
#endif
int  CORE_dsyrfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 const double *A,    int LDA,
                 const double *T,    int LDT,
                       double *C,    int LDC,
                       double *WORK, int LDWORK);
void CORE_dlacpy(PLASMA_enum uplo, int M, int N,
                 const double *A, int LDA,
                       double *B, int LDB);
int CORE_dlacpy_pivot( const PLASMA_desc descA,
                       PLASMA_enum direct,
                       int k1, int k2, const int *ipiv,
                       int *rankin, int *rankout,
                       double *A, int lda,
                       int init);
void CORE_dlange(int norm, int M, int N,
                 const double *A, int LDA,
                 double *work, double *normA);
#ifdef COMPLEX
void CORE_dlansy(int norm, PLASMA_enum uplo, int N,
                 const double *A, int LDA,
                 double *work, double *normA);
#endif
void CORE_dlansy(int norm, PLASMA_enum uplo, int N,
                 const double *A, int LDA,
                 double *work, double *normA);
void CORE_dlantr(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N,
                 const double *A, int LDA,
                 double *work, double *normA);
int CORE_dlarfb_gemm(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                     int M, int N, int K,
                     const double *V, int LDV,
                     const double *T, int LDT,
                           double *C, int LDC,
                           double *WORK, int LDWORK);
int CORE_dlarfx2(PLASMA_enum side, int N,
                 double V,
                 double TAU,
                 double *C1, int LDC1,
                 double *C2, int LDC2);
int CORE_dlarfx2c(PLASMA_enum uplo,
                  double V,
                  double TAU,
                  double *C1,
                  double *C2,
                  double *C3);
int CORE_dlarfx2ce(PLASMA_enum uplo,
                   double *V,
                   double *TAU,
                   double *C1,
                   double *C2,
                   double *C3);
void CORE_dlarfy(int N,
                 double *A, int LDA,
                 const double *V,
                 const double *TAU,
                 double *WORK);
int  CORE_dlascal(PLASMA_enum uplo, int m, int n,
                  double alpha, double *A, int lda);
void CORE_dlaset(PLASMA_enum uplo, int n1, int n2,
                 double alpha, double beta,
                 double *tileA, int ldtilea);
void CORE_dlaset2(PLASMA_enum uplo, int n1, int n2, double alpha,
                  double *tileA, int ldtilea);
void CORE_dlaswp(int N, double *A, int LDA,
                 int I1,  int I2, const int *IPIV, int INC);
int  CORE_dlaswp_ontile( PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_dlaswpc_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_dlatro(PLASMA_enum uplo, PLASMA_enum trans,
                 int M, int N,
                 const double *A, int LDA,
                       double *B, int LDB);
void CORE_dlauum(PLASMA_enum uplo, int N, double *A, int LDA);
int CORE_dpamm(int op, PLASMA_enum side, PLASMA_enum storev,
               int M, int N, int K, int L,
               const double *A1, int LDA1,
                     double *A2, int LDA2,
               const double *V, int LDV,
                     double *W, int LDW);
int  CORE_dparfb(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                 int M1, int N1, int M2, int N2, int K, int L,
                       double *A1, int LDA1,
                       double *A2, int LDA2,
                 const double *V, int LDV,
                 const double *T, int LDT,
                       double *WORK, int LDWORK);
int CORE_dpemv(PLASMA_enum trans, PLASMA_enum storev,
               int M, int N, int L,
               double ALPHA,
               const double *A, int LDA,
               const double *X, int INCX,
               double BETA,
               double *Y, int INCY,
               double *WORK);
void CORE_dplgsy(double bump, int m, int n, double *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_dplgsy(double bump, int m, int n, double *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_dplrnt(int m, int n, double *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
int  CORE_dpltmg(PLASMA_enum mtxtype, int m, int n, double *A, int lda,
                  int gM, int gN, int m0, int n0, unsigned long long int seed );
int  CORE_dpltmg_chebvand( int M, int N, double *A, int LDA,
                           int gN, int m0, int n0,
                           double *W );
int  CORE_dpltmg_circul( int M, int N, double *A, int LDA,
                         int gM, int m0, int n0,
                         const double *V );
void CORE_dpltmg_condexq( int M, int N, double *Q, int LDQ );
void CORE_dpltmg_fiedler(int m, int n,
                         const double *X, int incX,
                         const double *Y, int incY,
                         double *A, int lda);
int  CORE_dpltmg_hankel( PLASMA_enum uplo, int M, int N, double *A, int LDA,
                         int m0, int n0, int nb,
                         const double *V1,
                         const double *V2 );
void CORE_dpltmg_toeppd1( int gM, int m0, int M, double *W,
                          unsigned long long int seed );
void CORE_dpltmg_toeppd2( int M, int N, int K, int m0, int n0,
                          const double *W,
                          double *A, int LDA );
void CORE_dpotrf(PLASMA_enum uplo, int N, double *A, int LDA, int *INFO);
void CORE_dsetvar(const double *alpha, double *x);
void CORE_dshift(int s, int m, int n, int L,
                 double *A);
void CORE_dshiftw(int s, int cl, int m, int n, int L,
                  double *A, double *W);
int  CORE_dssssm(int M1, int N1, int M2, int N2, int K, int IB,
                       double *A1, int LDA1,
                       double *A2, int LDA2,
                 const double *L1, int LDL1,
                 const double *L2, int LDL2,
                 const int *IPIV);
int CORE_dstedc(PLASMA_enum compz, int n,
                double *D, double *E,
                double *Z, int LDZ,
                double *WORK, int LWORK,
#ifdef COMPLEX
                double *RWORK, int LRWORK,
#endif
                int *IWORK, int LIWORK);
int CORE_dsteqr(PLASMA_enum compz, int n,
                double *D, double *E,
                double *Z, int LDZ,
                double *WORK);
void CORE_dsymm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                double alpha, const double *A, int LDA,
                                          const double *B, int LDB,
                double beta,        double *C, int LDC);
void CORE_dsyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                double alpha, const double *A, int LDA,
                double beta,        double *C, int LDC);
void CORE_dsyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 double alpha, const double *A, int LDA,
                                           const double *B, int LDB,
                 double beta,        double *C, int LDC);
int  CORE_dsyssq(PLASMA_enum uplo, int N,
                 const double *A, int LDA,
                 double *scale, double *sumsq);
void CORE_dswpab(int i, int n1, int n2,
                 double *A, double *work);
int  CORE_dswptr_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc,
                        const double *Akk, int ldak);
int CORE_dtradd(PLASMA_enum uplo, PLASMA_enum trans, int M, int N,
                      double alpha,
                const double *A, int LDA,
                      double beta,
                      double *B, int LDB);
void CORE_dtrasm(PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N, const double *A, int lda, double *work);
void CORE_dtrdalg1(int n,
                        int nb,
                        double *A,
                        int lda,
                        double *V,
                        double *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        double *work);
void CORE_dtrmm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                double alpha, const double *A, int LDA,
                                                double *B, int LDB);
void CORE_dtrsm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                double alpha, const double *A, int LDA,
                                                double *B, int LDB);
int  CORE_dtrssq(PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                 const double *A, int LDA,
                 double *scale, double *sumsq);
void CORE_dtrtri(PLASMA_enum uplo, PLASMA_enum diag, int N,
                 double *A, int LDA, int *info);
int  CORE_dtslqt(int M, int N, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *T, int LDT,
                 double *TAU, double *WORK);
int  CORE_dtsmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 const double *V, int LDV,
                 const double *T, int LDT,
                 double *WORK, int LDWORK);
int CORE_dtsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *A3, int lda3,
                        const double *V, int ldv,
                        const double *T, int ldt,
                        double *WORK, int ldwork);
int CORE_dtsmlq_sytra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        const double *V, int ldv,
                        const double *T, int ldt,
                        double *WORK, int ldwork);
int  CORE_dtsmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 const double *V, int LDV,
                 const double *T, int LDT,
                 double *WORK, int LDWORK);
int CORE_dtsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *A3, int lda3,
                        const double *V, int ldv,
                        const double *T, int ldt,
                        double *WORK, int ldwork);
int CORE_dtsmqr_sytra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        const double *V, int ldv,
                        const double *T, int ldt,
                        double *WORK, int ldwork);
int  CORE_dtsqrt(int M, int N, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *T, int LDT,
                 double *TAU, double *WORK);
int  CORE_dtstrf(int M, int N, int IB, int NB,
                 double *U, int LDU,
                 double *A, int LDA,
                 double *L, int LDL,
                 int *IPIV, double *WORK,
                 int LDWORK, int *INFO);
int  CORE_dttmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 const double *V, int LDV,
                 const double *T, int LDT,
                 double *WORK, int LDWORK);
int  CORE_dttqrt(int M, int N, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *T, int LDT,
                 double *TAU,
                 double *WORK);
int  CORE_dttmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 const double *V, int LDV,
                 const double *T, int LDT,
                 double *WORK, int LDWORK);
int  CORE_dttlqt(int M, int N, int IB,
                 double *A1, int LDA1,
                 double *A2, int LDA2,
                 double *T, int LDT,
                 double *TAU,
                 double *WORK);
int  CORE_dormlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const double *V, int LDV,
                 const double *T, int LDT,
                 double *C, int LDC,
                 double *WORK, int LDWORK);
int  CORE_dormqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const double *V, int LDV,
                 const double *T, int LDT,
                 double *C, int LDC,
                 double *WORK, int LDWORK);

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
void CORE_dswap(int m, int n, double *Q, int ldq,
                const double *work, const int *perm,
                int start, int end);
int CORE_dlascl(PLASMA_enum type, int kl, int ku, double cfrom, double cto,
                int m, int n, double *A, int lda);
#ifdef COMPLEX
int CORE_dlag2z(int m, int n, const double *Q, int LDQ,
                 double *Z, int LDZ);
#endif

#ifndef COMPLEX
void CORE_dlaed3_freebigwork(int oper, double **WORK);
void CORE_dlaed0_betaapprox(int subpbs, const int *subpbs_info,
                            double *D, const double *E);
int CORE_dlapst(PLASMA_enum type, int n,
                const double *D, int *INDX);
#endif

#if defined(QUARK_H)
/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by PLASMA) - alphabetical order
 **/
void QUARK_CORE_dasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                       const double *A, int lda, int szeA,
                       double *work, int szeW);
void QUARK_CORE_dasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                          const double *A, int lda, int szeA,
                          double *work, int szeW,
                          double *fake, int szeF);
void QUARK_CORE_dgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum trans, int m, int n, int nb,
                       double  alpha,
                       const double *A, int lda,
                       double  beta,
                       double *B, int ldb);
void QUARK_CORE_dbrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo,
                        int n, int nb,
                        double *A,
                        int lda,
                        double *VQ,
                        double *TAUQ,
                        double *VP,
                        double *TAUP,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_dgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt);
void QUARK_CORE_dgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int k, int nb,
                      double alpha, const double *A, int lda,
                      const double *B, int ldb,
                      double beta, double *C, int ldc);
void QUARK_CORE_dgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum transA, PLASMA_enum transB,
                        int m, int n, int k, int nb,
                        double alpha, const double *A, int lda,
                        const double *B, int ldb,
                        double beta, double *C, int ldc);
void QUARK_CORE_dgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         double alpha, const double *A, int lda,
                         const double *B, int ldb,
                         double beta, double *C, int ldc,
                         double *fake1, int szefake1, int flag1,
                         double *fake2, int szefake2, int flag2);
void QUARK_CORE_dgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         double alpha, const double *A, int lda,
                         const double **B, int ldb,
                         double beta, double *C, int ldc);
void QUARK_CORE_dgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           double alpha, const double *A, int lda,
                           const double **B, int ldb,
                           double beta, double *C, int ldc,
                           double *fake1, int szefake1, int flag1);
void QUARK_CORE_dgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         double alpha, const double *A, int lda,
                         const double *B, int ldb,
                         double beta, double **C, int ldc);
void QUARK_CORE_dgemm_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           const double *alpha, const double *A, int lda,
                                                            const double *B, int ldb,
                           const double *beta,        double *C, int ldc,
                           const double *Alock,
                           const double *Block,
                           const double *Clock);
void QUARK_CORE_dgemv(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum trans, int m, int n,
                      double alpha, const double *A, int lda,
                                                const double *x, int incx,
                      double beta,        double *y, int incy);
void QUARK_CORE_dgemv_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum trans,
                           int m, int n,
                           const double *alpha, const double *A, int lda,
                                                            const double *x, int incx,
                           const double *beta,        double *y, int incy,
                           const double *Alock,
                           const double *xlock,
                           const double *ylock);
void QUARK_CORE_dgeqp3_init( Quark *quark, Quark_Task_Flags *task_flags,
                             int n, int *jpvt );
void QUARK_CORE_dgeqp3_larfg(Quark *quark, Quark_Task_Flags *task_flags,
                             PLASMA_desc A, int ii, int jj, int i, int j,
                             double *tau, double *beta );
void QUARK_CORE_dgeqp3_norms( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, int ioff, int joff, double *norms1, double *norms2 );
void QUARK_CORE_dgeqp3_pivot( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A,
                              double *F, int ldf,
                              int jj, int k, int *jpvt,
                              double *norms1, double *norms2, int *info );
void QUARK_CORE_dgeqp3_tntpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              double *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_dgeqp3_update( Quark *quark, Quark_Task_Flags *task_flags,
                               double *Ajj, int lda1,
                               double *Ajk, int lda2,
                               double *Fk,  int ldf,
                               int joff, int k, int koff, int nb,
                               double *norms1, double *norms2, int *info );
void QUARK_CORE_dgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt);
void QUARK_CORE_dgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       const int *IPIV,
                       const double *L, int ldl,
                       double *A, int lda);
void QUARK_CORE_dgessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );
void QUARK_CORE_dgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       double *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_dgetrf_incpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int ib, int nb,
                              double *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_dgetrf_nopiv(Quark *quark, Quark_Task_Flags *task_flags,
                             int m, int n, int ib, int nb,
                             double *A, int lda,
                             PLASMA_sequence *sequence, PLASMA_request *request,
                             int iinfo);
void QUARK_CORE_dgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_dgetrf_data_t *data, int m, int n, int nb,
                              double *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_dgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_dgetrf_data_t *data,
                              PLASMA_desc A, double *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_dgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, double *A, int szeA);
void QUARK_CORE_dgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, double *A, int szeA,
                           double *fake, int szeF, int paramF);
void QUARK_CORE_dgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, double *A, int szeA,
                           double *fake1, int szeF1, int paramF1,
                           double *fake2, int szeF2, int paramF2);
void QUARK_CORE_dsymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      double alpha, const double *A, int lda,
                      const double *B, int ldb,
                      double beta, double *C, int ldc);
void QUARK_CORE_dsygst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, PLASMA_enum uplo, int N,
                       double *A, int LDA,
                       double *B, int LDB,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_dsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      double alpha, const double *A, int lda,
                      double beta, double *C, int ldc);
void QUARK_CORE_dsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       double alpha, const double *A, int lda,
                       const double *B, int LDB,
                       double beta, double *C, int ldc);
void QUARK_CORE_dsyrfb(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo,
                       int n, int k, int ib, int nb,
                       const double *A, int lda,
                       const double *T, int ldt,
                       double *C, int ldc);
void QUARK_CORE_dhessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );
void QUARK_CORE_dlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int mb,
                       const double *A, int lda,
                       double *B, int ldb);
void QUARK_CORE_dlacpy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, int m, int n, int nb,
                          const double *A, int lda,
                          double *B, int ldb,
                          double *fake1, int szefake1, int flag1);
void QUARK_CORE_dlacpy_pivot(Quark *quark, Quark_Task_Flags *task_flags,
                             const PLASMA_desc descA,
                             PLASMA_enum direct,
                             int k1, int k2, const int *ipiv,
                             int *rankin, int *rankout,
                             double *A, int lda,
                             int pos, int init);
void QUARK_CORE_dlange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       const double *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_dlange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int M, int N,
                          const double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
#ifdef COMPLEX
void QUARK_CORE_dlansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const double *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_dlansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, PLASMA_enum uplo, int N,
                          const double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
#endif
void QUARK_CORE_dlansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const double *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_dlansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, PLASMA_enum uplo, int N,
                          const double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
void QUARK_CORE_dlantr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                       const double *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_dlantr_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                          const double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
void QUARK_CORE_dlascal(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int m, int n, int nb,
                        double alpha, double *A, int lda);
void QUARK_CORE_dlaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n1, int n2, double alpha,
                       double beta, double *tileA, int ldtilea);
void QUARK_CORE_dlaset2(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int n1, int n2, double alpha,
                        double *tileA, int ldtilea);
void QUARK_CORE_dlaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, double *A, int lda,
                       int i1,  int i2, const int *ipiv, int inc);
void QUARK_CORE_dlaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, double *A, int lda,
                          int i1,  int i2, const int *ipiv, int inc,
                          double *fake1, int szefake1, int flag1,
                          double *fake2, int szefake2, int flag2);
void QUARK_CORE_dlaswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, double *A,
                              int i1,  int i2, const int *ipiv, int inc, double *fakepanel);
void QUARK_CORE_dlaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, double *A,
                                 int i1,  int i2, const int *ipiv, int inc,
                                 double *fake1, int szefake1, int flag1,
                                 double *fake2, int szefake2, int flag2);
void QUARK_CORE_dlaswpc_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_desc descA, double *A,
                               int i1,  int i2, const int *ipiv, int inc, double *fakepanel);
void QUARK_CORE_dlatro(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int mb,
                       const double *A, int lda,
                       double *B, int ldb);
void QUARK_CORE_dlatro_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int mb,
                          const double *A, int lda,
                                double *B, int ldb,
                          double *fake1, int szefake1, int flag1);
void QUARK_CORE_dlauum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       double *A, int lda);
void QUARK_CORE_dplgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       double bump, int m, int n, double *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_dplgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       double bump, int m, int n, double *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_dplrnt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, double *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_dpltmg(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum mtxtype, int m, int n, double *A, int lda,
                        int gM, int gN, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_dpltmg_chebvand( Quark *quark, Quark_Task_Flags *task_flags,
                                 int M, int N, double *A, int LDA,
                                 int gN, int m0, int n0,
                                 double *W );
void QUARK_CORE_dpltmg_circul( Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, double *A, int LDA,
                               int gM, int m0, int n0,
                               const double *W );
void QUARK_CORE_dpltmg_fiedler(Quark *quark, Quark_Task_Flags *task_flags,
                               int m, int n,
                               const double *X, int incX,
                               const double *Y, int incY,
                               double *A, int lda);
void QUARK_CORE_dpltmg_hankel( Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_enum uplo, int M, int N, double *A, int LDA,
                               int m0, int n0, int nb,
                               const double *V1,
                               const double *V2);
void QUARK_CORE_dpltmg_toeppd1(Quark *quark, Quark_Task_Flags *task_flags,
                               int gM, int m0, int M,
                               double *W,
                               unsigned long long int seed);
void QUARK_CORE_dpltmg_toeppd2(Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, int K, int m0, int n0,
                               const double *W,
                               double *A, int LDA );
void QUARK_CORE_dpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_dsetvar(Quark *quark, Quark_Task_Flags *task_flags,
                        const double *alpha, double *x,
                        double *Alock);
void QUARK_CORE_dshift( Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int m, int n, int L,
                        double *A);
void QUARK_CORE_dshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L,
                        double *A, double *W);
void QUARK_CORE_dssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       const double *L1, int ldl1,
                       const double *L2, int ldl2,
                       const int *IPIV);
void QUARK_CORE_dstedc(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       double *D, double *E,
                       double *Z, int ldz);
void QUARK_CORE_dstedc_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum compz, int n,
                          double *D, double *E,
                          double *Z, int ldz,
                          void *fake1, int szefake1, int flag1,
                          void *fake2, int szefake2, int flag2);
void QUARK_CORE_dsteqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       double *D, double *E,
                       double *Z, int ldz);
void QUARK_CORE_dsymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      double alpha, const double *A, int lda,
                      const double *B, int ldb,
                      double beta, double *C, int ldc);
void QUARK_CORE_dsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      double alpha, const double *A, int lda,
                      double beta, double *C, int ldc);
void QUARK_CORE_dsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       double alpha, const double *A, int lda,
                       const double *B, int LDB,
                       double beta, double *C, int ldc);
void QUARK_CORE_dsyssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );
void QUARK_CORE_dswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       double *A, int szeA);
void QUARK_CORE_dswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, double *Aij,
                              int i1,  int i2, const int *ipiv, int inc,
                              const double *Akk, int ldak);
void QUARK_CORE_dtradd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       double  alpha,
                       const double *A, int lda,
                       double  beta,
                       double *B, int ldb);
void QUARK_CORE_dtrasm(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int m, int n,
                       const double *A, int lda, int szeA,
                       double *work, int szeW);
void QUARK_CORE_dtrasm_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int m, int n,
                          const double *A, int lda, int szeA,
                          double *work, int szeW,
                          double *fake, int szeF);
void QUARK_CORE_dtrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                        int n,
                        int nb,
                        double *A,
                        int lda,
                        double *V,
                        double *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_dtrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      double alpha, const double *A, int lda,
                      double *B, int ldb);
void QUARK_CORE_dtrmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                         int m, int n, int nb,
                         double alpha, const double *A, int lda,
                         double **B, int ldb);
void QUARK_CORE_dtrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      double alpha, const double *A, int lda,
                      double *B, int ldb);
void QUARK_CORE_dtrssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, PLASMA_enum diag,
                           int m, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );
void QUARK_CORE_dtrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum diag, int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_dtslqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);
void QUARK_CORE_dtsmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       const double *V, int ldv,
                       const double *T, int ldt);
void QUARK_CORE_dtsmlq_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              const double *V, int ldv,
                              const double *T, int ldt);
void QUARK_CORE_dtsmlq_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              double *A3, int lda3,
                              const double *V, int ldv,
                              const double *T, int ldt);
void QUARK_CORE_dtsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       const double *V, int ldv,
                       const double *T, int ldt);
void QUARK_CORE_dtsmqr_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              const double *V, int ldv,
                              const double *T, int ldt);
void QUARK_CORE_dtsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              double *A3, int lda3,
                              const double *V, int ldv,
                              const double *T, int ldt);
void QUARK_CORE_dtsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);
void QUARK_CORE_dtstrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *U, int ldu,
                       double *A, int lda,
                       double *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_dttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       const double *V, int ldv,
                       const double *T, int ldt);
void QUARK_CORE_dttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);
void QUARK_CORE_dttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       const double *V, int ldv,
                       const double *T, int ldt);
void QUARK_CORE_dttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt);
void QUARK_CORE_dpamm(Quark *quark, Quark_Task_Flags *task_flags,
                       int op, PLASMA_enum side, PLASMA_enum storev,
                       int m, int n, int k, int l,
                       const double *A1, int lda1,
                       double *A2, int lda2,
                       const double *V, int ldv,
                       double *W, int ldw);
void QUARK_CORE_dplssq( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, const double *A, double *result );
void QUARK_CORE_dormlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int ib,  int nb, int k,
                       const double *A, int lda,
                       const double *T, int ldt,
                       double *C, int ldc);
void QUARK_CORE_dormqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const double *A, int lda,
                       const double *T, int ldt,
                       double *C, int ldc);


void QUARK_CORE_dlascl(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum type, int kl, int ku, double cfrom, double cto,
                       int m, int n, double *A, int lda);
void QUARK_CORE_dlascl_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                            PLASMA_enum type, int kl, int ku, double *cfrom, double *cto,
                            int m, int n, double *A, int lda,
                            void *fake, int szefake, int flag);
void QUARK_CORE_dlaed0_lascl( Quark *quark, Quark_Task_Flags *task_flags,
                              int n, double *scale, double *D, double *E);
void QUARK_CORE_dlaed0_betaapprox(Quark *quark, Quark_Task_Flags *task_flags,
                                  int subpbs, const int *subpbs_info,
                                  double *D, const double *E);

#ifndef COMPLEX
void QUARK_CORE_dlaed2_computeK(Quark *quark, Quark_Task_Flags *task_flags,
                                int *K1, int n, int n1,
                                double *beta, double *D, double *Q, int LDQ,
                                double *Z, double *DLAMBDA, double *W,
                                int *INDX, int *INDXC, int *INDXP, int *INDXQ,
                                int *COLTYP,
                                double **Qmerge, int wsmode,
                                int *K2);

void QUARK_CORE_dlaed1_pipelined(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, const int *K,
                                 const int *INDX, const int *ctot,
                                 double *D, const double *beta,
                                 double *Q, int LDQ, double *Q2,
                                 const double *DLAMBDA, const double *W, double *Wred,
                                 int start, int end);
void QUARK_CORE_dlaed2_compressq(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, int start, int end,
                                 const int *INDX, const int *ctot,
                                 const double *Q, int LDQ,
                                 double *Q2, int *K);
void QUARK_CORE_dlaed4_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                            int n, const int *K,
                            double *D, const double *beta,
                            double **Q, const int *LDQ,
                            const double *DLAMBDA, const double *W, const int *INDX,
                            int start, int end,
                            PLASMA_sequence *sequence, PLASMA_request *request,
                            void *fakeQ, int flagfQ);
void QUARK_CORE_dlaed3_compW_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                                  int n, const int *K,
                                  double **Q, const int *LDQ,
                                  const double *DLAMBDA, double *W,
                                  const int *INDX,
                                  int start, int end,
                                  void *fakeQ, int flagfQ,
                                  void *fakeW, int flagfW);

void QUARK_CORE_dlaed3_reduceW(Quark *quark, Quark_Task_Flags *task_flags,
                               int n, int n1, const int *K, int l,
                               const double *Q, int LDQ,
                               const double *Wred, double *W);
void QUARK_CORE_dlaed3_reduceW_p2(Quark *quark, Quark_Task_Flags *task_flags,
                                  int n, int n1, const int *K, int l,
                                  double **Q, const int *LDQ,
                                  const double *Wred, double *W);

void QUARK_CORE_dlaed2_copydef(Quark *quark, Quark_Task_Flags *task_flags,
                               int n, int n1, const int *K, const int *ctot,
                               double *Q, int LDQ, const double *Q2,
                               int start, int end);
void QUARK_CORE_dlaed3_computevectors(Quark *quark, Quark_Task_Flags *task_flags,
                                      int wsmode, int n, const int *K,
                      const int *il_nondef, const int *iu_nondef,
                                      double *Q, int LDQ, double *W, const int *INDXC,
                                      double **WSglobal, double **WSlocal,
                                      int start, int end );
void QUARK_CORE_dlaed3_wscopy( Quark *quark, Quark_Task_Flags *task_flags,
                               const int *K, const int *il_nondef, const int *iu_nondef,
                               const double *Q, int LDQ, double **WORK,
                               int start, int end );
void QUARK_CORE_dlaed3_updatevectors(Quark *quark, Quark_Task_Flags *task_flags,
                                     int oper, int wsmode, int n, int n1, int *K,
                     int *il_nondef, int *iu_nondef,
                                     double *D, double *Q, int LDQ, double *Q2,
                                     int *INDXQ, int *COLTYP, double **WORK,
                                     int start, int end, double **WORKDEP);
void QUARK_CORE_dlaed3_pipelined(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, int *K, int *il_nondef, int *iu_nondef,
                                 double *D, double *Q, int LDQ, double *Q2,
                                 int *INDXC, int *INDXQ, int *COLTYP, double *W,
                                 int start, int end2);

void QUARK_CORE_dDC_fakedep(Quark *quark, Quark_Task_Flags *task_flags,
                            int nb_tasks, int nb, double *Q, int LDQ, double *W);
#endif

void QUARK_CORE_dswap(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, double *Q,
                      int LDQ, double *work,
                      int *perm, int begin, int end);
#ifdef COMPLEX
void QUARK_CORE_dlag2z(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n,
                       const double *Q, int LDQ,
                       double *Z, int LDZ);
#endif
void QUARK_CORE_dlaed3_freebigwork(Quark *quark, Quark_Task_Flags *task_flags,
                   int *K_bis, int largework, double **WORK);
void QUARK_CORE_dlaset_identity(Quark *quark, Quark_Task_Flags *task_flags,
                int n, int start, int size,
                double *A);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_dasum_quark(Quark *quark);
void CORE_dasum_f1_quark(Quark *quark);
void CORE_dgeadd_quark(Quark *quark);
void CORE_dbrdalg1_quark(Quark *quark);
void CORE_dgelqt_quark(Quark *quark);
void CORE_dgemm_quark(Quark *quark);
void CORE_dgemm_tile_quark(Quark *quark);
void CORE_dgemv_quark(Quark *quark);
void CORE_dgemv_tile_quark(Quark *quark);
void CORE_dgeqp3_init_quark(Quark *quark);
void CORE_dgeqp3_larfg_quark(Quark *quark);
void CORE_dgeqp3_norms_quark(Quark *quark);
void CORE_dgeqp3_pivot_quark(Quark *quark);
void CORE_dgeqp3_tntpiv_quark(Quark *quark);
void CORE_dgeqp3_update_quark(Quark *quark);
void CORE_dgeqrt_quark(Quark *quark);
void CORE_dgessm_quark(Quark *quark);
void CORE_dgessq_quark(Quark *quark);
void CORE_dgessq_f1_quark(Quark *quark);
void CORE_dgetrf_quark(Quark *quark);
void CORE_dgetrf_incpiv_quark(Quark *quark);
void CORE_dgetrf_nopiv_quark(Quark* quark);
void CORE_dgetrf_reclap_quark(Quark *quark);
void CORE_dgetrf_rectil_quark(Quark* quark);
void CORE_dgetrip_quark(Quark *quark);
void CORE_dgetrip_f1_quark(Quark *quark);
void CORE_dgetrip_f2_quark(Quark *quark);
#ifdef COMPLEX
void CORE_dsymm_quark(Quark *quark);
void CORE_dsyrk_quark(Quark *quark);
void CORE_dsyr2k_quark(Quark *quark);
#endif
void CORE_dsygst_quark(Quark *quark);
void CORE_dsyrfb_quark(Quark *quark);
void CORE_dhessq_quark(Quark *quark);
void CORE_dhessq_f1_quark(Quark *quark);
void CORE_dlacpy_quark(Quark *quark);
void CORE_dlacpy_f1_quark(Quark *quark);
void CORE_dlacpy_pivot_quark(Quark *quark);
void CORE_dlatro_quark(Quark *quark);
void CORE_dlatro_f1_quark(Quark *quark);
void CORE_dlange_quark(Quark *quark);
void CORE_dlange_f1_quark(Quark *quark);
#ifdef COMPLEX
void CORE_dlansy_quark(Quark *quark);
void CORE_dlansy_f1_quark(Quark *quark);
#endif
void CORE_dlansy_quark(Quark *quark);
void CORE_dlansy_f1_quark(Quark *quark);
void CORE_dlaset_quark(Quark *quark);
void CORE_dlaset2_quark(Quark *quark);
void CORE_dlatro_quark(Quark *quark);
void CORE_dlauum_quark(Quark *quark);
void CORE_dpamm_quark(Quark *quark);
void CORE_dplgsy_quark(Quark *quark);
void CORE_dplgsy_quark(Quark *quark);
void CORE_dplrnt_quark(Quark *quark);
void CORE_dpltmg_quark(Quark *quark);
void CORE_dplssq_quark(Quark *quark);
void CORE_dpotrf_quark(Quark *quark);
void CORE_dsetvar_quark(Quark *quark);
void CORE_dshift_quark(Quark *quark);
void CORE_dshiftw_quark(Quark *quark);
void CORE_dssssm_quark(Quark *quark);
void CORE_dsymm_quark(Quark *quark);
void CORE_dsyrk_quark(Quark *quark);
void CORE_dsyr2k_quark(Quark *quark);
void CORE_dsyssq_quark(Quark *quark);
void CORE_dsyssq_f1_quark(Quark *quark);
void CORE_dswpab_quark(Quark *quark);
void CORE_dswptr_ontile_quark(Quark *quark);
void CORE_dtrdalg1_quark(Quark *quark);
void CORE_dtrmm_quark(Quark *quark);
void CORE_dtrsm_quark(Quark *quark);
void CORE_dtrtri_quark(Quark *quark);
void CORE_dtslqt_quark(Quark *quark);
void CORE_dtsmlq_quark(Quark *quark);
void CORE_dtsmlq_sytra1_quark(Quark *quark);
void CORE_dtsmlq_corner_quark(Quark *quark);
void CORE_dtsmqr_quark(Quark *quark);
void CORE_dtsmqr_sytra1_quark(Quark *quark);
void CORE_dtsmqr_corner_quark(Quark *quark);
void CORE_dtsqrt_quark(Quark *quark);
void CORE_dtstrf_quark(Quark *quark);
void CORE_dttmqr_quark(Quark *quark);
void CORE_dttqrt_quark(Quark *quark);
void CORE_dttmlq_quark(Quark *quark);
void CORE_dttlqt_quark(Quark *quark);
void CORE_dormlq_quark(Quark *quark);
void CORE_dormqr_quark(Quark *quark);
void CORE_dlaswp_quark(Quark* quark);
void CORE_dlaswp_f2_quark(Quark* quark);
void CORE_dlaswp_ontile_quark(Quark *quark);
void CORE_dlaswp_ontile_f2_quark(Quark *quark);
void CORE_dlaswpc_ontile_quark(Quark *quark);
void CORE_dtrmm_p2_quark(Quark* quark);
void CORE_dgemm_f2_quark(Quark* quark);
void CORE_dgemm_p2_quark(Quark* quark);
void CORE_dgemm_p2f1_quark(Quark* quark);
void CORE_dgemm_p3_quark(Quark* quark);

#endif /* defined(QUARK_H) */

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
