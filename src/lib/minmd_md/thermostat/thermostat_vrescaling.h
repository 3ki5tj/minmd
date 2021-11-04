#ifndef THERMOSTAT_VRESCALING_H__
#define THERMOSTAT_VRESCALING_H__

#include "utils.h"
#include "rng.h"
#include "mdutils.h"
#include <math.h>
#include "thermostat_basic.h"


typedef struct {
  rng_t *rng;
} thermostat_vrescaling_param_t;

typedef struct {
  double halfkt; /* 0.5*kT */
  double expndt; /* exp(-dt) */
} thermostat_vrescaling_data_t;


INLINE thermostat_vrescaling_data_t *thermostat_vrescaling_data_new(
    const thermostat_param_t *param)
{
  thermostat_vrescaling_data_t *vrd;
  XNEW(vrd, 1);
  vrd->halfkt = 0.5 * param->boltz * param->tp;
  vrd->expndt = exp(-param->dt);
  return vrd;
}


INLINE void thermostat_vrescaling_data_delete(thermostat_vrescaling_data_t *vsd)
{
  free(vsd);
}


thermostat_t *thermostat_vrescaling_new(thermostat_param_t *param)
{
  thermostat_t *ts;

  XNEW(ts, 1);
  ts->type = THERMOSTAT_TYPE_VRESCALING;

  /* clone parameters */
  XCLONE(ts->param, param, sizeof(*param));
  XCLONE(ts->param->iparam, param->iparam, sizeof(thermostat_vrescaling_param_t));

  /* initialize data */
  thermostat_vrescaling_data_t *vrd = thermostat_vrescaling_data_new(ts->param);
  ts->data = thermostat_data_new(ts->param, vrd);

  return ts;
}


void thermostat_vrescaling_delete(thermostat_t *ts)
{
  free(ts->param->iparam);
  free(ts->param);

  thermostat_vrescaling_data_delete((thermostat_vrescaling_data_t *)ts->data->idata);
  thermostat_data_delete(ts->data);

  free(ts);
}


real thermostat_vrescaling_apply(thermostat_t *ts)
{
  thermostat_param_t *tsp = ts->param;
  thermostat_vrescaling_param_t *vrp = (thermostat_vrescaling_param_t *) tsp->iparam;
  thermostat_data_t *tsd = ts->data;
  thermostat_vrescaling_data_t *vrd = (thermostat_vrescaling_data_t *) tsd->idata;
  double dt = tsp->dt,
         c = vrd->expndt,
         halfkt = vrd->halfkt,
         ek1, ek2, s, r, r2;
  int n = tsp->n, i;
  rng_t *rng = vrp->rng;

  ek1 = mdutils_ekin(tsp->n, tsp->mass, tsp->v);
  r = rng_gauss(rng);
  r2 = rng_chisqr(rng, tsp->n_dof - 1);
  ek2 = ek1
      + (1 - c) * ((r2 + r * r) * halfkt - ek1)
      + 2 * r * sqrt(c * (1 - c) * ek1 * halfkt);
  if (ek2 < 0) {
    ek2 = 0;
  }
  s = sqrt(ek2 / ek1);
  for (i = 0; i < n; i++) {
    vec_smul(tsp->v[i], (real) s);
  }
  return (real) ek2;
}


#endif /* THERMOSTAT_VRESCALING_H__ */


