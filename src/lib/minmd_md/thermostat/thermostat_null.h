#ifndef THERMOSTAT_NULL_H__
#define THERMOSTAT_NULL_H__

/* dummy thermostat that does nothing */

#include "thermostat_basic.h"
#include "utils.h"
#include "mdutils.h"


thermostat_t *thermostat_null_init(thermostat_param_t *param)
{
  thermostat_t *ts;

  XNEW(ts, 1);
  ts->type = THERMOSTAT_TYPE_NULL;

  /* clone parameters */
  XCLONE(ts->param, param, sizeof(*param));
  ts->param->algo_param = NULL;
  ts->data = NULL;

  /* initialize data */
  ts->data = thermostat_data_init(ts->param, NULL);

  return ts;
}


void thermostat_null_free(thermostat_t *ts)
{
  free(ts->param);
  free(ts->data);
  free(ts);
}


real thermostat_null_apply(thermostat_t *ts)
{
  thermostat_param_t *tsp = ts->param;

  return mdutils_ekin(tsp->n, tsp->mass, tsp->v);
}

#endif /* THERMOSTAT_NULL_H__ */


