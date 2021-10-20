#ifndef VRESCALING_H__
#define VRESCALING_H__

#include "utils.h"
#include "rng.h"
#include "mdutils.h"
#include <math.h>
#include "thermostat_basic.h"


typedef struct {
  rng_t *rng;
} vrescaling_param_t;

typedef struct {
  double expndt; /* exp(-dt) */
} vrescaling_data_t;


thermostat_t *vrescaling_init(thermostat_param_t *param)
{
  thermostat_t *ts;
  vrescaling_data_t *vr_data;

  XNEW(ts, 1);
  ts->type = THERMOSTAT_TYPE_VRESCALING;
  ts->param = param;

  XNEW(vr_data, 1);
  vr_data->expndt = exp(-param->dt);
  ts->data = vr_data;
  return ts;
}


void vrescaling_free(thermostat_t *ts)
{
  thermostat_free(ts);
}


real vrescaling_apply(thermostat_t *ts)
{
  thermostat_param_t *tsp = ts->param;
  vrescaling_param_t *vrp = (vrescaling_param_t *) tsp->algo_param;
  vrescaling_data_t *vrd = (vrescaling_data_t *) ts->data;
  double dt = tsp->dt,
         c = vrd->expndt,
         ek1, ek2, s, r, r2;
  int n = tsp->n, i;
  rng_t *rng = vrp->rng;

  ek1 = mdutils_ekin(tsp->n, tsp->mass, tsp->v);
  r = rng_gauss(rng);
  r2 = rng_chisqr(rng, tsp->n_dof - 1);
  ek2 = ek1
      + (1 - c) * ((r2 + r * r) * tsp->tp * 0.5 - ek1)
      + 2 * r * sqrt(c * (1 - c) * ek1 * tsp->tp * 0.5);
  if (ek2 < 0) {
    ek2 = 0;
  }
  s = sqrt(ek2 / ek1);
  for (i = 0; i < n; i++) {
    vec_smul(tsp->v[i], (real) s);
  }
  return (real) ek2;
}

#endif /* VRESCALING_H__ */


