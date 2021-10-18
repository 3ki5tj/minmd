#ifndef STAT_ACCUM_H__
#define STAT_ACCUM_H__

#include "def.h"
#include "minmd_utils.h"
#include <math.h>


typedef struct {
  long long count;
  double sum;
  double svar;
} stat_accum_t;


INLINE stat_accum_t *stat_accum_init()
{
  stat_accum_t *acc;

  XNEW(acc, 1);
  acc->count = 0;
  acc->sum = 0;
  acc->svar = 0;
}


INLINE void stat_accum_free(stat_accum_t *acc)
{
  free(acc);
}


INLINE void stat_accum_add(stat_accum_t *acc, double x)
{
  long long n = acc->count;
  if (n > 0) {
    double xmean = acc->sum / n;
    double xdiff = x - xmean;
    acc->svar += xdiff*xdiff*n/(n+1);
  }
  acc->count = n + 1;
  acc->sum += x;
}


INLINE double stat_accum_get_mean(const stat_accum_t *acc)
{
  return (acc->count > 0) ? acc->sum / acc->count : 0;
}


INLINE double stat_accum_get_var(const stat_accum_t *acc, double *mean)
{
  long long n = acc->count;
  double xmean = 0.0, xvar = 0.0;

  if (n > 0) {
    xmean = acc->sum / n;
    xvar = acc->svar / n;
  }

  if (mean != NULL) {
    *mean = xmean;
  }
  return xvar;
}



INLINE double stat_accum_get_std(const stat_accum_t *acc, double *mean)
{
  double var = stat_accum_get_var(acc, mean);
  return sqrt(var);
}


#endif /* STAT_ACCUM_H__ */

