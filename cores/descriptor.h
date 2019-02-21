/**
 *
 * @file cores/descriptor.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_DESCRIPTOR_H_
#define _PLASMA_DESCRIPTOR_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Tile matrix cores/descriptor
 *
 *  Matrices are stored in a contiguous data chunk containing in order
 *  A11, A21, A12, A22 with :
 *
 *           n1      n2
 *      +----------+---+
 *      |          |   |    With m1 = lm - (lm%mb)
 *      |          |   |         m2 = lm%mb
 *  m1  |    A11   |A12|         n1 = ln - (ln%nb)
 *      |          |   |         n2 = ln%nb
 *      |          |   |
 *      +----------+---+
 *  m2  |    A21   |A22|
 *      +----------+---+
 *
 */
typedef struct plasma_desc_t {
    void *mat;          /**< pointer to the beginning of the matrix                           */
    size_t A21;         /**< pointer to the beginning of the matrix A21                       */
    size_t A12;         /**< pointer to the beginning of the matrix A12                       */
    size_t A22;         /**< pointer to the beginning of the matrix A22                       */
    PLASMA_enum dtyp;   /**< precision of the matrix                                          */
    int mb;             /**< number of rows in a tile                                         */
    int nb;             /**< number of columns in a tile                                      */
    int bsiz;           /**< size in elements including padding                               */
    int lm;             /**< number of rows of the entire matrix                              */
    int ln;             /**< number of columns of the entire matrix                           */
    int lm1;            /**< number of tile rows of the A11 matrix - derived parameter        */
    int ln1;            /**< number of tile columns of the A11 matrix - derived parameter     */
    int lmt;            /**< number of tile rows of the entire matrix - derived parameter     */
    int lnt;            /**< number of tile columns of the entire matrix - derived parameter  */
    int i;              /**< row index to the beginning of the submatrix                      */
    int j;              /**< column index to the beginning of the submatrix                   */
    int m;              /**< number of rows of the submatrix                                  */
    int n;              /**< number of columns of the submatrix                               */
    int mt;             /**< number of tile rows of the submatrix - derived parameter         */
    int nt;             /**< number of tile columns of the submatrix - derived parameter      */
} PLASMA_desc;


/***************************************************************************//**
 *  Internal routines used by coreblas library
 **/
PLASMA_desc plasma_desc_init(PLASMA_enum dtyp, int mb, int nb, int bsiz,
                             int lm, int ln, int i, int j, int m, int n);
PLASMA_desc plasma_desc_submatrix(PLASMA_desc descA, int i, int j, int m, int n);


static inline int plasma_element_size(int type)
{
    switch(type) {
    case PlasmaByte:          return          1;
    case PlasmaInteger:       return   sizeof(int);
    case PlasmaRealFloat:     return   sizeof(float);
    case PlasmaRealDouble:    return   sizeof(double);
    case PlasmaComplexFloat:  return 2*sizeof(float);
    case PlasmaComplexDouble: return 2*sizeof(double);
    default:
        fprintf(stderr, "plasma_element_size: invalide type parameter\n");
        return -1;
    }
}

/***************************************************************************//**
 *  Internal function to return adress of block (m,n)
 **/
inline static void *plasma_getaddr(PLASMA_desc A, int m, int n)
{
    size_t mm = m+A.i/A.mb;
    size_t nn = n+A.j/A.nb;
    size_t eltsize = plasma_element_size(A.dtyp);
    size_t offset = 0;

    if (mm < (size_t)(A.lm1)) {
        if (nn < (size_t)(A.ln1))
            offset = A.bsiz*(mm + (size_t)A.lm1 * nn);
        else
            offset = A.A12 + ((size_t)A.mb * (A.ln%A.nb) * mm);
    }
    else {
        if (nn < (size_t)(A.ln1))
            offset = A.A21 + ((size_t)A.nb * (A.lm%A.mb) * nn);
        else
            offset = A.A22;
    }

    return (void*)((char*)A.mat + (offset*eltsize) );
}

/***************************************************************************//**
 *  Internal function to return adress of element A(m,n)
 **/
inline static void *plasma_geteltaddr( const PLASMA_desc *A, int m, int n, int eltsize)
{
    size_t mm = m/A->mb;
    size_t nn = n/A->nb;
    size_t offset = 0;

    if (mm < (size_t)(A->lm1)) {
        if (nn < (size_t)(A->ln1))
            offset = A->bsiz*(mm+A->lm1*nn) + m%A->mb + A->mb * (size_t)(n%A->nb);
        else
            offset = A->A12 + (A->mb*(A->ln%A->nb)*mm) + m%A->mb + A->mb * (size_t)(n%A->nb);
    }
    else {
        if (nn < (size_t)(A->ln1))
            offset = A->A21 + ((A->lm%A->mb)*A->nb*nn) + m%A->mb + (A->lm%A->mb) * (size_t)(n%A->nb);
        else
            offset = A->A22 + m%A->mb  + (A->lm%A->mb) * (size_t)(n%A->nb);
    }
    return (void*)((char*)A->mat + (offset*eltsize) );
}

#ifdef __cplusplus
}
#endif

#endif
