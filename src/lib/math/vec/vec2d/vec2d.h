#ifndef VEC2D_H__
#define VEC2D_H__

#include "def.h"
#include <math.h>

#define DIM 2

typedef real vec[DIM];

INLINE real *vec_zero(real *x)
{
  x[0] = 0;
  x[1] = 0;
  return x;
}


INLINE real *vec_neg(real *x)
{
  x[0] = -x[0];
  x[1] = -x[1];
  return x;
}


INLINE real *vec_copy(real *x, const real *y)
{
  x[0] = y[0];
  x[1] = y[1];
  return x;
}


INLINE void vec_swap(real *x, real *y)
{
  real z[DIM];
  z[0] = x[0], x[0] = y[0], y[0] = z[0];
  z[1] = x[1], x[1] = y[1], y[1] = z[1];
}


INLINE real *vec_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  return x;
}


INLINE real *vec_vsmul(real *x, const real *y, real s)
{
  x[0] = y[0] * s;
  x[1] = y[1] * s;
  return x;
}


#define vec_inc(x, y) vec_sinc(x, y,  1)
#define vec_dec(x, y) vec_sinc(x, y, -1)

INLINE real *vec_sinc(real *x, const real *y, real s)
{
  x[0] += y[0] * s;
  x[1] += y[1] * s;
  return x;
}


#define vec_add(x, y, z) vec_sadd(x, y, z,  1)
#define vec_diff(x, y, z) vec_sadd(x, y, z, -1)

INLINE real *vec_sadd(real *x, const real *y, const real *z, real s)
{
  x[0] = y[0] + z[0] * s;
  x[1] = y[1] + z[1] * s;
  return x;
}


INLINE real *vec_lincomb(real *x, const real *y, const real *z, real s, real t)
{
  x[0] = y[0] * s + z[0] * t;
  x[1] = y[1] * s + z[1] * t;
  return x;
}


#define vec_sqr(x) vec_dot(x, x)

INLINE real vec_dot(const real *x, const real *y)
{
  return x[0]*y[0] + x[1]*y[1];
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


INLINE real vec_cross(real *x, real *y)
{
  return x[0]*y[1] - x[1]*y[0];
}


/* wrap a periodic x to (0, l) */
INLINE real *vec_wrap(real *x, real l)
{
  x[0] = (real) fmod(x[0] + 1000.*l, l);
  x[1] = (real) fmod(x[1] + 1000.*l, l);
  return x;
}




#endif /* VEC2D_H__ */

