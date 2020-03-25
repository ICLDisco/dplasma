/**
 *
 * @file primes.h
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

#ifndef PRIMES_H
#define PRIMES_H

#define IMBALANCE_THRESHOLD 10
#define PWR_MAXSIZE   32
#define PRIME_MAXSIZE 10
#define SIZE_MG       1024
#define SIZE_LEADERS  1023

#ifndef min
#define min(a,b) ((a<b)?a:b)
#endif

#ifndef max
#define max(a,b) ((a>b)?a:b)
#endif

/**
 * @ingroup InPlaceTransformation
 *
 * struct primedec - Structure to describes an integer as a product of prime
 * powers.
 *
 */
struct primedec
{
  int p;
  int e;
  int pe;
};

typedef struct primedec primedec_t;

int lcm(int a, int b);
int gcd(int a, int b);
int modpow(int x, int n, int m);
void factor(int n, primedec_t *pr, int *nf);

int64_t maxval(int n, int *T);
int64_t sum   (int n, int *T);

#endif /* PRIMES_H */
