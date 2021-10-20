#ifndef VEC3D_H__
#define VEC3D_H__

#include <math.h>
#include "rng.h"

#ifndef DIM
#define DIM 3
#endif

typedef real vec[DIM];

INLINE real *vec_zero(real *x)
{
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  return x;
}


INLINE real *vec_neg(real *x)
{
  x[0] = -x[0];
  x[1] = -x[1];
  x[2] = -x[2];
  return x;
}


INLINE real *vec_copy(real *x, const real *y)
{
  x[0] = y[0];
  x[1] = y[1];
  x[2] = y[2];
  return x;
}


INLINE void vec_swap(real *x, real *y)
{
  real z[DIM];
  z[0] = x[0], x[0] = y[0], y[0] = z[0];
  z[1] = x[1], x[1] = y[1], y[1] = z[1];
  z[2] = x[2], x[2] = y[2], y[2] = z[2];
}


INLINE real *vec_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return x;
}


INLINE real *vec_vsmul(real *x, const real *y, real s)
{
  x[0] = y[0] * s;
  x[1] = y[1] * s;
  x[2] = y[2] * s;
  return x;
}


#define vec_inc(x, y) vec_sinc(x, y,  1)
#define vec_dec(x, y) vec_sinc(x, y, -1)

INLINE real *vec_sinc(real *x, const real *y, real s)
{
  x[0] += y[0] * s;
  x[1] += y[1] * s;
  x[2] += y[2] * s;
  return x;
}


#define vec_add(x, y, z) vec_sadd(x, y, z,  1)
#define vec_diff(x, y, z) vec_sadd(x, y, z, -1)

INLINE real *vec_sadd(real *x, const real *y, const real *z, real s)
{
  x[0] = y[0] + z[0] * s;
  x[1] = y[1] + z[1] * s;
  x[2] = y[2] + z[2] * s;
  return x;
}


INLINE real *vec_lincomb(real *x, const real *y, const real *z, real s, real t)
{
  x[0] = y[0] * s + z[0] * t;
  x[1] = y[1] * s + z[1] * t;
  x[2] = y[2] * s + z[2] * t;
  return x;
}


#define vec_sqr(x) vec_dot(x, x)

INLINE real vec_dot(const real *x, const real *y)
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}


INLINE real vec_norm(const real *x)
{
  return (real) sqrt(vec_sqr(x));
}


INLINE real vec_distx(real *dx, const real *x, const real *y)
{
  return vec_norm(vec_diff(dx, x, y));
}


INLINE real vec_dist(const real *x, const real *y)
{
  real dx[DIM];
  return vec_distx(dx, x, y);
}


INLINE real vec_dist2x(real *dx, const real *x, const real *y)
{
  return vec_sqr(vec_diff(dx, x, y));
}


INLINE real vec_dist2(const real *x, const real *y)
{
  real dx[DIM];
  return vec_dist2x(dx, x, y);
}


/* (in-place) normalize the vector */
INLINE real *vec_normalize(real *v)
{
  double s = vec_sqr(v);
  s = (s >= 0) ? 1./sqrt(s) : 1;
  return vec_smul(v, s);
}


INLINE real *vec_cross(real *z, const real *x, const real *y)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
  return z;
}


/* wrap a periodic x to (0, l) */
INLINE real *vec_wrap(real *x, real l)
{
  x[0] = (real) fmod(x[0] + 1000.*l, l);
  x[1] = (real) fmod(x[1] + 1000.*l, l);
  x[2] = (real) fmod(x[2] + 1000.*l, l);
  return x;
}


INLINE real *vec_rand_gauss(real *v, real mag, rng_t *rng)
{
  v[0] = mag * rng_gauss(rng);
  v[1] = mag * rng_gauss(rng);
  v[2] = mag * rng_gauss(rng);
  return v;
}


/* a randomly oriented unit vector */
INLINE real *vec_rand_dir(real *v, rng_t *r)
{
  double a, b, sqr, s;

  do {
    a = 2 * rng_rand01(r) - 1;
    b = 2 * rng_rand01(r) - 1;
    sqr = a*a + b*b;
  } while (sqr >= 1);

  s = 2.*sqrt(1 - sqr);
  v[0] = (real) (a * s);
  v[1] = (real) (b * s);
  v[2] = (real) (1 - 2*sqr);
  return v;
}


#endif /* VEC3D_H__ */

