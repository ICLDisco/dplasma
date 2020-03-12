/*
 * Copyright (c) 2009-2017 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */

#include "dplasma/complex.h"

#ifndef DPLASMA_HAVE_COMPLEX_H

#include <math.h>

float cabsf(float _Complex z)
{
    float *zp = (float *)&z;
    float x,y,ans,temp;

    /* Numerical Recipes in C, Second Edition, p 949 */
    x=fabsf(zp[0]);
    y=fabsf(zp[1]);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrtf(1.0f+temp*temp);
    } else {
        temp=x/y;
        ans=y*sqrtf(1.0f+temp*temp);
    }
    return ans;
}

double cabs(double _Complex z)
{
    double *zp = (double *)&z;
    double x,y,ans,temp;

    /* Numerical Recipes in C, Second Edition, p 949 */
    x=fabs(zp[0]);
    y=fabs(zp[1]);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
    } else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
    }
    return ans;
}

double cimag(dplasma_complex64_t z)
{
    return ((double *)&z)[1];
}

double creal(dplasma_complex64_t z)
{
    return ((double *)&z)[0];
}

dplasma_complex64_t conj(dplasma_complex64_t z)
{
    double *zp, *vp;
    dplasma_complex64_t v;

    zp = (double *)&z;
    vp = (double *)&v;
    vp[0] = zp[0];
    vp[1] = -zp[1];
    return v;
}

#endif /* DPLASMA_HAS_COMPLEX */
