/**
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
void CORE_dzasum(int storev, PLASMA_enum uplo, int M, int N,
                 const parsec_complex64_t *A, int lda, double *work);
void CORE_zbrdalg1( PLASMA_enum uplo,
                    int n,
                    int nb,
                    parsec_complex64_t *A,
                    int lda,
                    parsec_complex64_t *VQ,
                    parsec_complex64_t *TAUQ,
                    parsec_complex64_t *VP,
                    parsec_complex64_t *TAUP,
                    int Vblksiz, int wantz,
                    int i, int sweepid, int m, int grsiz,
                    parsec_complex64_t *work);
int CORE_zgbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex64_t *V, parsec_complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zgbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex64_t *V, parsec_complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zgblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex64_t *V, parsec_complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zgeadd(PLASMA_enum trans, int M, int N,
                      parsec_complex64_t alpha,
                const parsec_complex64_t *A, int LDA,
                      parsec_complex64_t beta,
                      parsec_complex64_t *B, int LDB);
int  CORE_zgelqt(int M, int N, int IB,
                 parsec_complex64_t *A, int LDA,
                 parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *TAU,
                 parsec_complex64_t *WORK);
void CORE_zgemm(PLASMA_enum transA, PLASMA_enum transB,
                int M, int N, int K,
                parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                                          const parsec_complex64_t *B, int LDB,
                parsec_complex64_t beta,        parsec_complex64_t *C, int LDC);
void CORE_zgemv(PLASMA_enum trans, int M, int N,
                parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                                          const parsec_complex64_t *x, int incx,
                parsec_complex64_t beta,        parsec_complex64_t *y, int incy);
void CORE_zgeqp3_init( int n, int *jpvt );
void CORE_zgeqp3_larfg( PLASMA_desc A, int ii, int jj, int i, int j,
                        parsec_complex64_t *tau, parsec_complex64_t *beta );
void CORE_zgeqp3_norms( PLASMA_desc A, int ioff, int joff, double *norms1, double *norms2 );
void CORE_zgeqp3_pivot( PLASMA_desc A, parsec_complex64_t *F, int ldf,
                        int jj, int k, int *jpvt,
                        double *norms1, double *norms2, int *info );
int  CORE_zgeqp3_tntpiv(int m, int n,
                        parsec_complex64_t *A, int lda,
                        int *IPIV, parsec_complex64_t *tau,
                        int *iwork);
void CORE_zgeqp3_update( const parsec_complex64_t *Ajj, int lda1,
                         parsec_complex64_t       *Ajk, int lda2,
                         const parsec_complex64_t *Fk,  int ldf,
                         int joff, int k, int koff, int nb,
                         double *norms1, double *norms2,
                         int *info );
int  CORE_zgeqrt(int M, int N, int IB,
                 parsec_complex64_t *A, int LDA,
                 parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *TAU, parsec_complex64_t *WORK);
int  CORE_zgessm(int M, int N, int K, int IB,
                 const int *IPIV,
                 const parsec_complex64_t *L, int LDL,
                 parsec_complex64_t *A, int LDA);
int  CORE_zgessq(int M, int N,
                 const parsec_complex64_t *A, int LDA,
                 double *scale, double *sumsq);
int  CORE_zgetf2_nopiv(int m, int n,
                      parsec_complex64_t *A, int lda);
int  CORE_zgetrf(int M, int N,
                 parsec_complex64_t *A, int LDA,
                 int *IPIV, int *INFO);
int  CORE_zgetrf_incpiv(int M, int N, int IB,
                        parsec_complex64_t *A, int LDA,
                        int *IPIV, int *INFO);
int  CORE_zgetrf_nopiv(int m, int n, int ib,
                      parsec_complex64_t *A, int lda);
int  CORE_zgetrf_reclap(CORE_zgetrf_data_t *data, int M, int N,
                        parsec_complex64_t *A, int LDA,
                        int *IPIV, int *info);
CORE_zgetrf_data_t *CORE_zgetrf_reclap_init(int nbthrd);
int  CORE_zgetrf_rectil(CORE_zgetrf_data_t *data, const PLASMA_desc A, int *IPIV, int *info);
CORE_zgetrf_data_t *CORE_zgetrf_rectil_init(int nbthrd);
void CORE_zgetrip(int m, int n, parsec_complex64_t *A,
                  parsec_complex64_t *work);
int CORE_zhbelr(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex64_t *V, parsec_complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zhblrx(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex64_t *V, parsec_complex64_t *TAU,
                int st, int ed, int eltsize);
int CORE_zhbrce(PLASMA_enum uplo, int N,
                PLASMA_desc *A, parsec_complex64_t *V, parsec_complex64_t *TAU,
                int st, int ed, int eltsize);
void CORE_zhbtype1cb(int N, int NB,
                     parsec_complex64_t *A, int LDA,
                     parsec_complex64_t *V, parsec_complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     parsec_complex64_t *WORK);
void CORE_zhbtype2cb(int N, int NB,
                     parsec_complex64_t *A, int LDA,
                     parsec_complex64_t *V, parsec_complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     parsec_complex64_t *WORK);
void CORE_zhbtype3cb(int N, int NB,
                     parsec_complex64_t *A, int LDA,
                     const parsec_complex64_t *V, const parsec_complex64_t *TAU,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     parsec_complex64_t *WORK);
void CORE_zgbtype1cb(PLASMA_enum uplo, int N, int NB,
                parsec_complex64_t *A, int LDA,
                parsec_complex64_t *VQ, parsec_complex64_t *TAUQ,
                parsec_complex64_t *VP, parsec_complex64_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                parsec_complex64_t *WORK);
void CORE_zgbtype2cb(PLASMA_enum uplo, int N, int NB,
                parsec_complex64_t *A, int LDA,
                parsec_complex64_t *VQ, parsec_complex64_t *TAUQ,
                parsec_complex64_t *VP, parsec_complex64_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                parsec_complex64_t *WORK);
void CORE_zgbtype3cb(PLASMA_enum uplo, int N, int NB,
                parsec_complex64_t *A, int LDA,
                parsec_complex64_t *VQ, parsec_complex64_t *TAUQ,
                parsec_complex64_t *VP, parsec_complex64_t *TAUP,
                int st, int ed, int sweep, int Vblksiz, int WANTZ,
                parsec_complex64_t *WORK);
void CORE_zhegst(int itype, PLASMA_enum uplo, int N,
                 parsec_complex64_t *A, int LDA,
                 parsec_complex64_t *B, int LDB, int *INFO);
#ifdef COMPLEX
void CORE_zhemm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                                          const parsec_complex64_t *B, int LDB,
                parsec_complex64_t beta,        parsec_complex64_t *C, int LDC);
void CORE_zherk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                double alpha, const parsec_complex64_t *A, int LDA,
                double beta,        parsec_complex64_t *C, int LDC);
void CORE_zher2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                                           const parsec_complex64_t *B, int LDB,
                 double beta,                    parsec_complex64_t *C, int LDC);
int  CORE_zhessq(PLASMA_enum uplo, int N,
                 const parsec_complex64_t *A, int LDA,
                 double *scale, double *sumsq);
#endif
int  CORE_zherfb(PLASMA_enum uplo, int N, int K, int IB, int NB,
                 const parsec_complex64_t *A,    int LDA,
                 const parsec_complex64_t *T,    int LDT,
                       parsec_complex64_t *C,    int LDC,
                       parsec_complex64_t *WORK, int LDWORK);
void CORE_zlacpy(PLASMA_enum uplo, int M, int N,
                 const parsec_complex64_t *A, int LDA,
                       parsec_complex64_t *B, int LDB);
int CORE_zlacpy_pivot( const PLASMA_desc descA,
                       PLASMA_enum direct,
                       int k1, int k2, const int *ipiv,
                       int *rankin, int *rankout,
                       parsec_complex64_t *A, int lda,
                       int init);
void CORE_zlange(int norm, int M, int N,
                 const parsec_complex64_t *A, int LDA,
                 double *work, double *normA);
#ifdef COMPLEX
void CORE_zlanhe(int norm, PLASMA_enum uplo, int N,
                 const parsec_complex64_t *A, int LDA,
                 double *work, double *normA);
#endif
void CORE_zlansy(int norm, PLASMA_enum uplo, int N,
                 const parsec_complex64_t *A, int LDA,
                 double *work, double *normA);
void CORE_zlantr(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N,
                 const parsec_complex64_t *A, int LDA,
                 double *work, double *normA);
int CORE_zlarfb_gemm(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                     int M, int N, int K,
                     const parsec_complex64_t *V, int LDV,
                     const parsec_complex64_t *T, int LDT,
                           parsec_complex64_t *C, int LDC,
                           parsec_complex64_t *WORK, int LDWORK);
int CORE_zlarfx2(PLASMA_enum side, int N,
                 parsec_complex64_t V,
                 parsec_complex64_t TAU,
                 parsec_complex64_t *C1, int LDC1,
                 parsec_complex64_t *C2, int LDC2);
int CORE_zlarfx2c(PLASMA_enum uplo,
                  parsec_complex64_t V,
                  parsec_complex64_t TAU,
                  parsec_complex64_t *C1,
                  parsec_complex64_t *C2,
                  parsec_complex64_t *C3);
int CORE_zlarfx2ce(PLASMA_enum uplo,
                   parsec_complex64_t *V,
                   parsec_complex64_t *TAU,
                   parsec_complex64_t *C1,
                   parsec_complex64_t *C2,
                   parsec_complex64_t *C3);
void CORE_zlarfy(int N,
                 parsec_complex64_t *A, int LDA,
                 const parsec_complex64_t *V,
                 const parsec_complex64_t *TAU,
                 parsec_complex64_t *WORK);
int  CORE_zlascal(PLASMA_enum uplo, int m, int n,
                  parsec_complex64_t alpha, parsec_complex64_t *A, int lda);
void CORE_zlaset(PLASMA_enum uplo, int n1, int n2,
                 parsec_complex64_t alpha, parsec_complex64_t beta,
                 parsec_complex64_t *tileA, int ldtilea);
void CORE_zlaset2(PLASMA_enum uplo, int n1, int n2, parsec_complex64_t alpha,
                  parsec_complex64_t *tileA, int ldtilea);
void CORE_zlaswp(int N, parsec_complex64_t *A, int LDA,
                 int I1,  int I2, const int *IPIV, int INC);
int  CORE_zlaswp_ontile( PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_zlaswpc_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc);
int  CORE_zlatro(PLASMA_enum uplo, PLASMA_enum trans,
                 int M, int N,
                 const parsec_complex64_t *A, int LDA,
                       parsec_complex64_t *B, int LDB);
void CORE_zlauum(PLASMA_enum uplo, int N, parsec_complex64_t *A, int LDA);
int CORE_zpamm(int op, PLASMA_enum side, PLASMA_enum storev,
               int M, int N, int K, int L,
               const parsec_complex64_t *A1, int LDA1,
                     parsec_complex64_t *A2, int LDA2,
               const parsec_complex64_t *V, int LDV,
                     parsec_complex64_t *W, int LDW);
int  CORE_zparfb(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                 int M1, int N1, int M2, int N2, int K, int L,
                       parsec_complex64_t *A1, int LDA1,
                       parsec_complex64_t *A2, int LDA2,
                 const parsec_complex64_t *V, int LDV,
                 const parsec_complex64_t *T, int LDT,
                       parsec_complex64_t *WORK, int LDWORK);
int CORE_zpemv(PLASMA_enum trans, PLASMA_enum storev,
               int M, int N, int L,
               parsec_complex64_t ALPHA,
               const parsec_complex64_t *A, int LDA,
               const parsec_complex64_t *X, int INCX,
               parsec_complex64_t BETA,
               parsec_complex64_t *Y, int INCY,
               parsec_complex64_t *WORK);
void CORE_zplghe(double bump, int m, int n, parsec_complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_zplgsy(parsec_complex64_t bump, int m, int n, parsec_complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
void CORE_zplrnt(int m, int n, parsec_complex64_t *A, int lda,
                 int bigM, int m0, int n0, unsigned long long int seed );
int  CORE_zpltmg(PLASMA_enum mtxtype, int m, int n, parsec_complex64_t *A, int lda,
                  int gM, int gN, int m0, int n0, unsigned long long int seed );
int  CORE_zpltmg_chebvand( int M, int N, parsec_complex64_t *A, int LDA,
                           int gN, int m0, int n0,
                           parsec_complex64_t *W );
int  CORE_zpltmg_circul( int M, int N, parsec_complex64_t *A, int LDA,
                         int gM, int m0, int n0,
                         const parsec_complex64_t *V );
void CORE_zpltmg_condexq( int M, int N, parsec_complex64_t *Q, int LDQ );
void CORE_zpltmg_fiedler(int m, int n,
                         const parsec_complex64_t *X, int incX,
                         const parsec_complex64_t *Y, int incY,
                         parsec_complex64_t *A, int lda);
int  CORE_zpltmg_hankel( PLASMA_enum uplo, int M, int N, parsec_complex64_t *A, int LDA,
                         int m0, int n0, int nb,
                         const parsec_complex64_t *V1,
                         const parsec_complex64_t *V2 );
void CORE_zpltmg_toeppd1( int gM, int m0, int M, parsec_complex64_t *W,
                          unsigned long long int seed );
void CORE_zpltmg_toeppd2( int M, int N, int K, int m0, int n0,
                          const parsec_complex64_t *W,
                          parsec_complex64_t *A, int LDA );
void CORE_zpotrf(PLASMA_enum uplo, int N, parsec_complex64_t *A, int LDA, int *INFO);
void CORE_zsetvar(const parsec_complex64_t *alpha, parsec_complex64_t *x);
void CORE_zshift(int s, int m, int n, int L,
                 parsec_complex64_t *A);
void CORE_zshiftw(int s, int cl, int m, int n, int L,
                  parsec_complex64_t *A, parsec_complex64_t *W);
int  CORE_zssssm(int M1, int N1, int M2, int N2, int K, int IB,
                       parsec_complex64_t *A1, int LDA1,
                       parsec_complex64_t *A2, int LDA2,
                 const parsec_complex64_t *L1, int LDL1,
                 const parsec_complex64_t *L2, int LDL2,
                 const int *IPIV);
int CORE_zstedc(PLASMA_enum compz, int n,
                double *D, double *E,
                parsec_complex64_t *Z, int LDZ,
                parsec_complex64_t *WORK, int LWORK,
#ifdef COMPLEX
                double *RWORK, int LRWORK,
#endif
                int *IWORK, int LIWORK);
int CORE_zsteqr(PLASMA_enum compz, int n,
                double *D, double *E,
                parsec_complex64_t *Z, int LDZ,
                double *WORK);
void CORE_zsymm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                                          const parsec_complex64_t *B, int LDB,
                parsec_complex64_t beta,        parsec_complex64_t *C, int LDC);
void CORE_zsyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                parsec_complex64_t beta,        parsec_complex64_t *C, int LDC);
void CORE_zsyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                                           const parsec_complex64_t *B, int LDB,
                 parsec_complex64_t beta,        parsec_complex64_t *C, int LDC);
int  CORE_zsyssq(PLASMA_enum uplo, int N,
                 const parsec_complex64_t *A, int LDA,
                 double *scale, double *sumsq);
void CORE_zswpab(int i, int n1, int n2,
                 parsec_complex64_t *A, parsec_complex64_t *work);
int  CORE_zswptr_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc,
                        const parsec_complex64_t *Akk, int ldak);
int CORE_ztradd(PLASMA_enum uplo, PLASMA_enum trans, int M, int N,
                      parsec_complex64_t alpha,
                const parsec_complex64_t *A, int LDA,
                      parsec_complex64_t beta,
                      parsec_complex64_t *B, int LDB);
void CORE_ztrasm(PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag,
                 int M, int N, const parsec_complex64_t *A, int lda, double *work);
void CORE_ztrdalg1(int n,
                        int nb,
                        parsec_complex64_t *A,
                        int lda,
                        parsec_complex64_t *V,
                        parsec_complex64_t *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        parsec_complex64_t *work);
void CORE_ztrmm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                                                parsec_complex64_t *B, int LDB);
void CORE_ztrsm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                parsec_complex64_t alpha, const parsec_complex64_t *A, int LDA,
                                                parsec_complex64_t *B, int LDB);
int  CORE_ztrssq(PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                 const parsec_complex64_t *A, int LDA,
                 double *scale, double *sumsq);
void CORE_ztrtri(PLASMA_enum uplo, PLASMA_enum diag, int N,
                 parsec_complex64_t *A, int LDA, int *info);
int  CORE_ztslqt(int M, int N, int IB,
                 parsec_complex64_t *A1, int LDA1,
                 parsec_complex64_t *A2, int LDA2,
                 parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *TAU, parsec_complex64_t *WORK);
int  CORE_ztsmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 parsec_complex64_t *A1, int LDA1,
                 parsec_complex64_t *A2, int LDA2,
                 const parsec_complex64_t *V, int LDV,
                 const parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *WORK, int LDWORK);
int CORE_ztsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        parsec_complex64_t *A1, int lda1,
                        parsec_complex64_t *A2, int lda2,
                        parsec_complex64_t *A3, int lda3,
                        const parsec_complex64_t *V, int ldv,
                        const parsec_complex64_t *T, int ldt,
                        parsec_complex64_t *WORK, int ldwork);
int CORE_ztsmlq_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        parsec_complex64_t *A1, int lda1,
                        parsec_complex64_t *A2, int lda2,
                        const parsec_complex64_t *V, int ldv,
                        const parsec_complex64_t *T, int ldt,
                        parsec_complex64_t *WORK, int ldwork);
int  CORE_ztsmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 parsec_complex64_t *A1, int LDA1,
                 parsec_complex64_t *A2, int LDA2,
                 const parsec_complex64_t *V, int LDV,
                 const parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *WORK, int LDWORK);
int CORE_ztsmqr_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        parsec_complex64_t *A1, int lda1,
                        parsec_complex64_t *A2, int lda2,
                        parsec_complex64_t *A3, int lda3,
                        const parsec_complex64_t *V, int ldv,
                        const parsec_complex64_t *T, int ldt,
                        parsec_complex64_t *WORK, int ldwork);
int CORE_ztsmqr_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        parsec_complex64_t *A1, int lda1,
                        parsec_complex64_t *A2, int lda2,
                        const parsec_complex64_t *V, int ldv,
                        const parsec_complex64_t *T, int ldt,
                        parsec_complex64_t *WORK, int ldwork);
int  CORE_ztsqrt(int M, int N, int IB,
                 parsec_complex64_t *A1, int LDA1,
                 parsec_complex64_t *A2, int LDA2,
                 parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *TAU, parsec_complex64_t *WORK);
int  CORE_ztstrf(int M, int N, int IB, int NB,
                 parsec_complex64_t *U, int LDU,
                 parsec_complex64_t *A, int LDA,
                 parsec_complex64_t *L, int LDL,
                 int *IPIV, parsec_complex64_t *WORK,
                 int LDWORK, int *INFO);
int  CORE_zttmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 parsec_complex64_t *A1, int LDA1,
                 parsec_complex64_t *A2, int LDA2,
                 const parsec_complex64_t *V, int LDV,
                 const parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *WORK, int LDWORK);
int  CORE_zttqrt(int M, int N, int IB,
                 parsec_complex64_t *A1, int LDA1,
                 parsec_complex64_t *A2, int LDA2,
                 parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *TAU,
                 parsec_complex64_t *WORK);
int  CORE_zttmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 parsec_complex64_t *A1, int LDA1,
                 parsec_complex64_t *A2, int LDA2,
                 const parsec_complex64_t *V, int LDV,
                 const parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *WORK, int LDWORK);
int  CORE_zttlqt(int M, int N, int IB,
                 parsec_complex64_t *A1, int LDA1,
                 parsec_complex64_t *A2, int LDA2,
                 parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *TAU,
                 parsec_complex64_t *WORK);
int  CORE_zunmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const parsec_complex64_t *V, int LDV,
                 const parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *C, int LDC,
                 parsec_complex64_t *WORK, int LDWORK);
int  CORE_zunmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const parsec_complex64_t *V, int LDV,
                 const parsec_complex64_t *T, int LDT,
                 parsec_complex64_t *C, int LDC,
                 parsec_complex64_t *WORK, int LDWORK);

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
void CORE_zswap(int m, int n, parsec_complex64_t *Q, int ldq,
                const parsec_complex64_t *work, const int *perm,
                int start, int end);
int CORE_zlascl(PLASMA_enum type, int kl, int ku, double cfrom, double cto,
                int m, int n, parsec_complex64_t *A, int lda);
#ifdef COMPLEX
int CORE_dlag2z(int m, int n, const double *Q, int LDQ,
                 parsec_complex64_t *Z, int LDZ);
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
void QUARK_CORE_dzasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                       const parsec_complex64_t *A, int lda, int szeA,
                       double *work, int szeW);
void QUARK_CORE_dzasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int m, int n,
                          const parsec_complex64_t *A, int lda, int szeA,
                          double *work, int szeW,
                          double *fake, int szeF);
void QUARK_CORE_zgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum trans, int m, int n, int nb,
                       parsec_complex64_t  alpha,
                       const parsec_complex64_t *A, int lda,
                       parsec_complex64_t  beta,
                       parsec_complex64_t *B, int ldb);
void QUARK_CORE_zbrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo,
                        int n, int nb,
                        parsec_complex64_t *A,
                        int lda,
                        parsec_complex64_t *VQ,
                        parsec_complex64_t *TAUQ,
                        parsec_complex64_t *VP,
                        parsec_complex64_t *TAUP,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_zgelqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex64_t *A, int lda,
                       parsec_complex64_t *T, int ldt);
void QUARK_CORE_zgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int k, int nb,
                      parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                      const parsec_complex64_t *B, int ldb,
                      parsec_complex64_t beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum transA, PLASMA_enum transB,
                        int m, int n, int k, int nb,
                        parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                        const parsec_complex64_t *B, int ldb,
                        parsec_complex64_t beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                         const parsec_complex64_t *B, int ldb,
                         parsec_complex64_t beta, parsec_complex64_t *C, int ldc,
                         parsec_complex64_t *fake1, int szefake1, int flag1,
                         parsec_complex64_t *fake2, int szefake2, int flag2);
void QUARK_CORE_zgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                         const parsec_complex64_t **B, int ldb,
                         parsec_complex64_t beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                           const parsec_complex64_t **B, int ldb,
                           parsec_complex64_t beta, parsec_complex64_t *C, int ldc,
                           parsec_complex64_t *fake1, int szefake1, int flag1);
void QUARK_CORE_zgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, PLASMA_enum transB,
                         int m, int n, int k, int nb,
                         parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                         const parsec_complex64_t *B, int ldb,
                         parsec_complex64_t beta, parsec_complex64_t **C, int ldc);
void QUARK_CORE_zgemm_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           const parsec_complex64_t *alpha, const parsec_complex64_t *A, int lda,
                                                            const parsec_complex64_t *B, int ldb,
                           const parsec_complex64_t *beta,        parsec_complex64_t *C, int ldc,
                           const parsec_complex64_t *Alock,
                           const parsec_complex64_t *Block,
                           const parsec_complex64_t *Clock);
void QUARK_CORE_zgemv(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum trans, int m, int n,
                      parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                                                const parsec_complex64_t *x, int incx,
                      parsec_complex64_t beta,        parsec_complex64_t *y, int incy);
void QUARK_CORE_zgemv_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum trans,
                           int m, int n,
                           const parsec_complex64_t *alpha, const parsec_complex64_t *A, int lda,
                                                            const parsec_complex64_t *x, int incx,
                           const parsec_complex64_t *beta,        parsec_complex64_t *y, int incy,
                           const parsec_complex64_t *Alock,
                           const parsec_complex64_t *xlock,
                           const parsec_complex64_t *ylock);
void QUARK_CORE_zgeqp3_init( Quark *quark, Quark_Task_Flags *task_flags,
                             int n, int *jpvt );
void QUARK_CORE_zgeqp3_larfg(Quark *quark, Quark_Task_Flags *task_flags,
                             PLASMA_desc A, int ii, int jj, int i, int j,
                             parsec_complex64_t *tau, parsec_complex64_t *beta );
void QUARK_CORE_zgeqp3_norms( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, int ioff, int joff, double *norms1, double *norms2 );
void QUARK_CORE_zgeqp3_pivot( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A,
                              parsec_complex64_t *F, int ldf,
                              int jj, int k, int *jpvt,
                              double *norms1, double *norms2, int *info );
void QUARK_CORE_zgeqp3_tntpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              parsec_complex64_t *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_zgeqp3_update( Quark *quark, Quark_Task_Flags *task_flags,
                               parsec_complex64_t *Ajj, int lda1,
                               parsec_complex64_t *Ajk, int lda2,
                               parsec_complex64_t *Fk,  int ldf,
                               int joff, int k, int koff, int nb,
                               double *norms1, double *norms2, int *info );
void QUARK_CORE_zgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex64_t *A, int lda,
                       parsec_complex64_t *T, int ldt);
void QUARK_CORE_zgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       const int *IPIV,
                       const parsec_complex64_t *L, int ldl,
                       parsec_complex64_t *A, int lda);
void QUARK_CORE_zgessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, const parsec_complex64_t *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );
void QUARK_CORE_zgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       parsec_complex64_t *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_zgetrf_incpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int ib, int nb,
                              parsec_complex64_t *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo);
void QUARK_CORE_zgetrf_nopiv(Quark *quark, Quark_Task_Flags *task_flags,
                             int m, int n, int ib, int nb,
                             parsec_complex64_t *A, int lda,
                             PLASMA_sequence *sequence, PLASMA_request *request,
                             int iinfo);
void QUARK_CORE_zgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_zgetrf_data_t *data, int m, int n, int nb,
                              parsec_complex64_t *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_zgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_zgetrf_data_t *data,
                              PLASMA_desc A, parsec_complex64_t *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread);
void QUARK_CORE_zgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, parsec_complex64_t *A, int szeA);
void QUARK_CORE_zgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, parsec_complex64_t *A, int szeA,
                           parsec_complex64_t *fake, int szeF, int paramF);
void QUARK_CORE_zgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, parsec_complex64_t *A, int szeA,
                           parsec_complex64_t *fake1, int szeF1, int paramF1,
                           parsec_complex64_t *fake2, int szeF2, int paramF2);
void QUARK_CORE_zhemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                      const parsec_complex64_t *B, int ldb,
                      parsec_complex64_t beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zhegst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, PLASMA_enum uplo, int N,
                       parsec_complex64_t *A, int LDA,
                       parsec_complex64_t *B, int LDB,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_zherk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      double alpha, const parsec_complex64_t *A, int lda,
                      double beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zher2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                       const parsec_complex64_t *B, int LDB,
                       double beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zherfb(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo,
                       int n, int k, int ib, int nb,
                       const parsec_complex64_t *A, int lda,
                       const parsec_complex64_t *T, int ldt,
                       parsec_complex64_t *C, int ldc);
void QUARK_CORE_zhessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const parsec_complex64_t *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );
void QUARK_CORE_zlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int mb,
                       const parsec_complex64_t *A, int lda,
                       parsec_complex64_t *B, int ldb);
void QUARK_CORE_zlacpy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, int m, int n, int nb,
                          const parsec_complex64_t *A, int lda,
                          parsec_complex64_t *B, int ldb,
                          parsec_complex64_t *fake1, int szefake1, int flag1);
void QUARK_CORE_zlacpy_pivot(Quark *quark, Quark_Task_Flags *task_flags,
                             const PLASMA_desc descA,
                             PLASMA_enum direct,
                             int k1, int k2, const int *ipiv,
                             int *rankin, int *rankout,
                             parsec_complex64_t *A, int lda,
                             int pos, int init);
void QUARK_CORE_zlange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       const parsec_complex64_t *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_zlange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, int M, int N,
                          const parsec_complex64_t *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
#ifdef COMPLEX
void QUARK_CORE_zlanhe(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const parsec_complex64_t *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_zlanhe_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, PLASMA_enum uplo, int N,
                          const parsec_complex64_t *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
#endif
void QUARK_CORE_zlansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const parsec_complex64_t *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_zlansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          int norm, PLASMA_enum uplo, int N,
                          const parsec_complex64_t *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
void QUARK_CORE_zlantr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                       const parsec_complex64_t *A, int LDA, int szeA,
                       int szeW, double *result);
void QUARK_CORE_zlantr_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                          const parsec_complex64_t *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF);
void QUARK_CORE_zlascal(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int m, int n, int nb,
                        parsec_complex64_t alpha, parsec_complex64_t *A, int lda);
void QUARK_CORE_zlaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n1, int n2, parsec_complex64_t alpha,
                       parsec_complex64_t beta, parsec_complex64_t *tileA, int ldtilea);
void QUARK_CORE_zlaset2(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int n1, int n2, parsec_complex64_t alpha,
                        parsec_complex64_t *tileA, int ldtilea);
void QUARK_CORE_zlaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, parsec_complex64_t *A, int lda,
                       int i1,  int i2, const int *ipiv, int inc);
void QUARK_CORE_zlaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, parsec_complex64_t *A, int lda,
                          int i1,  int i2, const int *ipiv, int inc,
                          parsec_complex64_t *fake1, int szefake1, int flag1,
                          parsec_complex64_t *fake2, int szefake2, int flag2);
void QUARK_CORE_zlaswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, parsec_complex64_t *A,
                              int i1,  int i2, const int *ipiv, int inc, parsec_complex64_t *fakepanel);
void QUARK_CORE_zlaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, parsec_complex64_t *A,
                                 int i1,  int i2, const int *ipiv, int inc,
                                 parsec_complex64_t *fake1, int szefake1, int flag1,
                                 parsec_complex64_t *fake2, int szefake2, int flag2);
void QUARK_CORE_zlaswpc_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_desc descA, parsec_complex64_t *A,
                               int i1,  int i2, const int *ipiv, int inc, parsec_complex64_t *fakepanel);
void QUARK_CORE_zlatro(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int mb,
                       const parsec_complex64_t *A, int lda,
                       parsec_complex64_t *B, int ldb);
void QUARK_CORE_zlatro_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int mb,
                          const parsec_complex64_t *A, int lda,
                                parsec_complex64_t *B, int ldb,
                          parsec_complex64_t *fake1, int szefake1, int flag1);
void QUARK_CORE_zlauum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       parsec_complex64_t *A, int lda);
void QUARK_CORE_zplghe(Quark *quark, Quark_Task_Flags *task_flags,
                       double bump, int m, int n, parsec_complex64_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_zplgsy(Quark *quark, Quark_Task_Flags *task_flags,
                       parsec_complex64_t bump, int m, int n, parsec_complex64_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_zplrnt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, parsec_complex64_t *A, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_zpltmg(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum mtxtype, int m, int n, parsec_complex64_t *A, int lda,
                        int gM, int gN, int m0, int n0, unsigned long long int seed );
void QUARK_CORE_zpltmg_chebvand( Quark *quark, Quark_Task_Flags *task_flags,
                                 int M, int N, parsec_complex64_t *A, int LDA,
                                 int gN, int m0, int n0,
                                 parsec_complex64_t *W );
void QUARK_CORE_zpltmg_circul( Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, parsec_complex64_t *A, int LDA,
                               int gM, int m0, int n0,
                               const parsec_complex64_t *W );
void QUARK_CORE_zpltmg_fiedler(Quark *quark, Quark_Task_Flags *task_flags,
                               int m, int n,
                               const parsec_complex64_t *X, int incX,
                               const parsec_complex64_t *Y, int incY,
                               parsec_complex64_t *A, int lda);
void QUARK_CORE_zpltmg_hankel( Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_enum uplo, int M, int N, parsec_complex64_t *A, int LDA,
                               int m0, int n0, int nb,
                               const parsec_complex64_t *V1,
                               const parsec_complex64_t *V2);
void QUARK_CORE_zpltmg_toeppd1(Quark *quark, Quark_Task_Flags *task_flags,
                               int gM, int m0, int M,
                               parsec_complex64_t *W,
                               unsigned long long int seed);
void QUARK_CORE_zpltmg_toeppd2(Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, int K, int m0, int n0,
                               const parsec_complex64_t *W,
                               parsec_complex64_t *A, int LDA );
void QUARK_CORE_zpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       parsec_complex64_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_zsetvar(Quark *quark, Quark_Task_Flags *task_flags,
                        const parsec_complex64_t *alpha, parsec_complex64_t *x,
                        parsec_complex64_t *Alock);
void QUARK_CORE_zshift( Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int m, int n, int L,
                        parsec_complex64_t *A);
void QUARK_CORE_zshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L,
                        parsec_complex64_t *A, parsec_complex64_t *W);
void QUARK_CORE_zssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       const parsec_complex64_t *L1, int ldl1,
                       const parsec_complex64_t *L2, int ldl2,
                       const int *IPIV);
void QUARK_CORE_zstedc(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       double *D, double *E,
                       parsec_complex64_t *Z, int ldz);
void QUARK_CORE_zstedc_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum compz, int n,
                          double *D, double *E,
                          parsec_complex64_t *Z, int ldz,
                          void *fake1, int szefake1, int flag1,
                          void *fake2, int szefake2, int flag2);
void QUARK_CORE_zsteqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       double *D, double *E,
                       parsec_complex64_t *Z, int ldz);
void QUARK_CORE_zsymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                      const parsec_complex64_t *B, int ldb,
                      parsec_complex64_t beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                      parsec_complex64_t beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                       const parsec_complex64_t *B, int LDB,
                       parsec_complex64_t beta, parsec_complex64_t *C, int ldc);
void QUARK_CORE_zsyssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const parsec_complex64_t *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );
void QUARK_CORE_zswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       parsec_complex64_t *A, int szeA);
void QUARK_CORE_zswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, parsec_complex64_t *Aij,
                              int i1,  int i2, const int *ipiv, int inc,
                              const parsec_complex64_t *Akk, int ldak);
void QUARK_CORE_ztradd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       parsec_complex64_t  alpha,
                       const parsec_complex64_t *A, int lda,
                       parsec_complex64_t  beta,
                       parsec_complex64_t *B, int ldb);
void QUARK_CORE_ztrasm(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int m, int n,
                       const parsec_complex64_t *A, int lda, int szeA,
                       double *work, int szeW);
void QUARK_CORE_ztrasm_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int m, int n,
                          const parsec_complex64_t *A, int lda, int szeA,
                          double *work, int szeW,
                          double *fake, int szeF);
void QUARK_CORE_ztrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                        int n,
                        int nb,
                        parsec_complex64_t *A,
                        int lda,
                        parsec_complex64_t *V,
                        parsec_complex64_t *TAU,
                        int Vblksiz, int wantz,
                        int i, int sweepid, int m, int grsiz,
                        int *PCOL, int *ACOL, int *MCOL);
void QUARK_CORE_ztrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                      parsec_complex64_t *B, int ldb);
void QUARK_CORE_ztrmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                         int m, int n, int nb,
                         parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                         parsec_complex64_t **B, int ldb);
void QUARK_CORE_ztrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      parsec_complex64_t alpha, const parsec_complex64_t *A, int lda,
                      parsec_complex64_t *B, int ldb);
void QUARK_CORE_ztrssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, PLASMA_enum diag,
                           int m, int n, const parsec_complex64_t *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF );
void QUARK_CORE_ztrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum diag, int n, int nb,
                       parsec_complex64_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo);
void QUARK_CORE_ztslqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       parsec_complex64_t *T, int ldt);
void QUARK_CORE_ztsmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       const parsec_complex64_t *V, int ldv,
                       const parsec_complex64_t *T, int ldt);
void QUARK_CORE_ztsmlq_hetra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              parsec_complex64_t *A1, int lda1,
                              parsec_complex64_t *A2, int lda2,
                              const parsec_complex64_t *V, int ldv,
                              const parsec_complex64_t *T, int ldt);
void QUARK_CORE_ztsmlq_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              parsec_complex64_t *A1, int lda1,
                              parsec_complex64_t *A2, int lda2,
                              parsec_complex64_t *A3, int lda3,
                              const parsec_complex64_t *V, int ldv,
                              const parsec_complex64_t *T, int ldt);
void QUARK_CORE_ztsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       const parsec_complex64_t *V, int ldv,
                       const parsec_complex64_t *T, int ldt);
void QUARK_CORE_ztsmqr_hetra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              parsec_complex64_t *A1, int lda1,
                              parsec_complex64_t *A2, int lda2,
                              const parsec_complex64_t *V, int ldv,
                              const parsec_complex64_t *T, int ldt);
void QUARK_CORE_ztsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              parsec_complex64_t *A1, int lda1,
                              parsec_complex64_t *A2, int lda2,
                              parsec_complex64_t *A3, int lda3,
                              const parsec_complex64_t *V, int ldv,
                              const parsec_complex64_t *T, int ldt);
void QUARK_CORE_ztsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       parsec_complex64_t *T, int ldt);
void QUARK_CORE_ztstrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex64_t *U, int ldu,
                       parsec_complex64_t *A, int lda,
                       parsec_complex64_t *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo);
void QUARK_CORE_zttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       const parsec_complex64_t *V, int ldv,
                       const parsec_complex64_t *T, int ldt);
void QUARK_CORE_zttqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       parsec_complex64_t *T, int ldt);
void QUARK_CORE_zttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       const parsec_complex64_t *V, int ldv,
                       const parsec_complex64_t *T, int ldt);
void QUARK_CORE_zttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       parsec_complex64_t *T, int ldt);
void QUARK_CORE_zpamm(Quark *quark, Quark_Task_Flags *task_flags,
                       int op, PLASMA_enum side, PLASMA_enum storev,
                       int m, int n, int k, int l,
                       const parsec_complex64_t *A1, int lda1,
                       parsec_complex64_t *A2, int lda2,
                       const parsec_complex64_t *V, int ldv,
                       parsec_complex64_t *W, int ldw);
void QUARK_CORE_zplssq( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, const double *A, double *result );
void QUARK_CORE_zunmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int ib,  int nb, int k,
                       const parsec_complex64_t *A, int lda,
                       const parsec_complex64_t *T, int ldt,
                       parsec_complex64_t *C, int ldc);
void QUARK_CORE_zunmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const parsec_complex64_t *A, int lda,
                       const parsec_complex64_t *T, int ldt,
                       parsec_complex64_t *C, int ldc);


void QUARK_CORE_zlascl(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum type, int kl, int ku, double cfrom, double cto,
                       int m, int n, parsec_complex64_t *A, int lda);
void QUARK_CORE_zlascl_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                            PLASMA_enum type, int kl, int ku, double *cfrom, double *cto,
                            int m, int n, parsec_complex64_t *A, int lda,
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

void QUARK_CORE_zswap(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n, parsec_complex64_t *Q,
                      int LDQ, parsec_complex64_t *work,
                      int *perm, int begin, int end);
#ifdef COMPLEX
void QUARK_CORE_dlag2z(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n,
                       const double *Q, int LDQ,
                       parsec_complex64_t *Z, int LDZ);
#endif
void QUARK_CORE_dlaed3_freebigwork(Quark *quark, Quark_Task_Flags *task_flags,
                   int *K_bis, int largework, double **WORK);
void QUARK_CORE_zlaset_identity(Quark *quark, Quark_Task_Flags *task_flags,
                int n, int start, int size,
                parsec_complex64_t *A);

/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_dzasum_quark(Quark *quark);
void CORE_dzasum_f1_quark(Quark *quark);
void CORE_zgeadd_quark(Quark *quark);
void CORE_zbrdalg1_quark(Quark *quark);
void CORE_zgelqt_quark(Quark *quark);
void CORE_zgemm_quark(Quark *quark);
void CORE_zgemm_tile_quark(Quark *quark);
void CORE_zgemv_quark(Quark *quark);
void CORE_zgemv_tile_quark(Quark *quark);
void CORE_zgeqp3_init_quark(Quark *quark);
void CORE_zgeqp3_larfg_quark(Quark *quark);
void CORE_zgeqp3_norms_quark(Quark *quark);
void CORE_zgeqp3_pivot_quark(Quark *quark);
void CORE_zgeqp3_tntpiv_quark(Quark *quark);
void CORE_zgeqp3_update_quark(Quark *quark);
void CORE_zgeqrt_quark(Quark *quark);
void CORE_zgessm_quark(Quark *quark);
void CORE_zgessq_quark(Quark *quark);
void CORE_zgessq_f1_quark(Quark *quark);
void CORE_zgetrf_quark(Quark *quark);
void CORE_zgetrf_incpiv_quark(Quark *quark);
void CORE_zgetrf_nopiv_quark(Quark* quark);
void CORE_zgetrf_reclap_quark(Quark *quark);
void CORE_zgetrf_rectil_quark(Quark* quark);
void CORE_zgetrip_quark(Quark *quark);
void CORE_zgetrip_f1_quark(Quark *quark);
void CORE_zgetrip_f2_quark(Quark *quark);
#ifdef COMPLEX
void CORE_zhemm_quark(Quark *quark);
void CORE_zherk_quark(Quark *quark);
void CORE_zher2k_quark(Quark *quark);
#endif
void CORE_zhegst_quark(Quark *quark);
void CORE_zherfb_quark(Quark *quark);
void CORE_zhessq_quark(Quark *quark);
void CORE_zhessq_f1_quark(Quark *quark);
void CORE_zlacpy_quark(Quark *quark);
void CORE_zlacpy_f1_quark(Quark *quark);
void CORE_zlacpy_pivot_quark(Quark *quark);
void CORE_zlatro_quark(Quark *quark);
void CORE_zlatro_f1_quark(Quark *quark);
void CORE_zlange_quark(Quark *quark);
void CORE_zlange_f1_quark(Quark *quark);
#ifdef COMPLEX
void CORE_zlanhe_quark(Quark *quark);
void CORE_zlanhe_f1_quark(Quark *quark);
#endif
void CORE_zlansy_quark(Quark *quark);
void CORE_zlansy_f1_quark(Quark *quark);
void CORE_zlaset_quark(Quark *quark);
void CORE_zlaset2_quark(Quark *quark);
void CORE_zlatro_quark(Quark *quark);
void CORE_zlauum_quark(Quark *quark);
void CORE_zpamm_quark(Quark *quark);
void CORE_zplghe_quark(Quark *quark);
void CORE_zplgsy_quark(Quark *quark);
void CORE_zplrnt_quark(Quark *quark);
void CORE_zpltmg_quark(Quark *quark);
void CORE_zplssq_quark(Quark *quark);
void CORE_zpotrf_quark(Quark *quark);
void CORE_zsetvar_quark(Quark *quark);
void CORE_zshift_quark(Quark *quark);
void CORE_zshiftw_quark(Quark *quark);
void CORE_zssssm_quark(Quark *quark);
void CORE_zsymm_quark(Quark *quark);
void CORE_zsyrk_quark(Quark *quark);
void CORE_zsyr2k_quark(Quark *quark);
void CORE_zsyssq_quark(Quark *quark);
void CORE_zsyssq_f1_quark(Quark *quark);
void CORE_zswpab_quark(Quark *quark);
void CORE_zswptr_ontile_quark(Quark *quark);
void CORE_ztrdalg1_quark(Quark *quark);
void CORE_ztrmm_quark(Quark *quark);
void CORE_ztrsm_quark(Quark *quark);
void CORE_ztrtri_quark(Quark *quark);
void CORE_ztslqt_quark(Quark *quark);
void CORE_ztsmlq_quark(Quark *quark);
void CORE_ztsmlq_hetra1_quark(Quark *quark);
void CORE_ztsmlq_corner_quark(Quark *quark);
void CORE_ztsmqr_quark(Quark *quark);
void CORE_ztsmqr_hetra1_quark(Quark *quark);
void CORE_ztsmqr_corner_quark(Quark *quark);
void CORE_ztsqrt_quark(Quark *quark);
void CORE_ztstrf_quark(Quark *quark);
void CORE_zttmqr_quark(Quark *quark);
void CORE_zttqrt_quark(Quark *quark);
void CORE_zttmlq_quark(Quark *quark);
void CORE_zttlqt_quark(Quark *quark);
void CORE_zunmlq_quark(Quark *quark);
void CORE_zunmqr_quark(Quark *quark);
void CORE_zlaswp_quark(Quark* quark);
void CORE_zlaswp_f2_quark(Quark* quark);
void CORE_zlaswp_ontile_quark(Quark *quark);
void CORE_zlaswp_ontile_f2_quark(Quark *quark);
void CORE_zlaswpc_ontile_quark(Quark *quark);
void CORE_ztrmm_p2_quark(Quark* quark);
void CORE_zgemm_f2_quark(Quark* quark);
void CORE_zgemm_p2_quark(Quark* quark);
void CORE_zgemm_p2f1_quark(Quark* quark);
void CORE_zgemm_p3_quark(Quark* quark);

#endif /* defined(QUARK_H) */

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
