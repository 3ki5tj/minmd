#ifndef RNG_H__
#define RNG_H__

#include "def.h"
#include "pcg32.h"
#include <math.h>
#include <time.h>

#define RNG_TYPE_PCG 0

typedef struct {
  int type;
  void *rng;
} rng_t;


rng_t *rng_init(int type, uint64_t seed)
{
  rng_t *r = (rng_t *) calloc(1, sizeof(rng_t));
  if (r == NULL) {
    fprintf(stderr, "(rng_init) Error: no memory\n");
    exit(1);
  }

  r->type = type;

  if (type == RNG_TYPE_PCG) {
    pcg32_random_t *pcg = pcg32_random_init(seed, 0);
    r->rng = pcg;
  }

  return r;
}

uint32_t rng_uint32(rng_t *r)
{
  if (r->type == RNG_TYPE_PCG) {
    pcg32_random_t *pcg = (pcg32_random_t *)(r->rng);
    return pcg32_random_uint32(pcg);
  } else {
    fprintf(stderr, "rng_uint32 Error: unknown RNG type %d\n", r->type);
    return 0;
  }
}

double rng_rand01(rng_t *r)
{
  return rng_uint32(r) / 4294967296.0;
}

/* return a normally distributed random number with zero mean and unit variance
 * using the ratio method
 * */
double rng_gauss(rng_t *r)
{
  double x, y, u, v, q;
  do {
    /* 1. Generate two independent uniform deviates U and V
     * 2. Compute X = sqrt(8/e) (V - 0.5) / U
     * 3. If X^2 <= -4 ln(U), then accept X, otherwise start over the algorithm.
     * Ref.
     * https://en.wikipedia.org/wiki/Normal_distribution#Computational_methods
     * */
    u = 1 - rng_rand01(r);
    v = 1.7156 * (rng_rand01(r) - .5);
    /* the following steps are used to avoid the evaluation of log(u) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));
  return v/u;
}


/* return gamma deviated random number
 * p(x) = x^{k-1} exp(-x) / (k-1)! */
double rng_gamma(rng_t *r, double k)
{
  int lt1 = 0;
  double a, b, x, v, u;

  if (k <= 0) return 0;
  if (k < 1) {
    lt1 = 1;
    k += 1;
  }
  a = k - 1./3;
  b = 1./3/sqrt(a);

  for (;;) {
    do {
      x = rng_gauss(r);
      v = 1 + b * x;
    } while (v <= 0);
    v *= v * v;
    x *= x;
    u = rng_rand01(r);
    if (u <= 1 - 0.331*x*x) {
      break;
    }
    u = log(u);
    if (u <= 0.5*x + a*(1-v+log(v))) {
      break;
    }
  }

  x = a*v;
  if (lt1) {
    x *= pow(1-rng_rand01(r), 1./(k-1));
  }
  return x;
}

double rng_chisqr(rng_t *r, double k)
{
  return 2*rng_gamma(r, k*0.5);
}

void rng_free(rng_t *r)
{
  if (r->type == RNG_TYPE_PCG) {
    pcg32_random_free( (pcg32_random_t *)(r->rng) );
  }
  free(r);
}

#endif /* RNG_H__ */

