#ifndef THERMOSTAT_NULL_H__
#define THERMOSTAT_NULL_H__

/* dummy thermostat that does nothing */

#include "thermostat_basic.h"
#include "utils.h"
#include "mdutils.h"


thermostat_t *thermostat_null_new(thermostat_param_t *param)
{
  thermostat_t *ts;

  XNEW(ts, 1);
  ts->type = THERMOSTAT_TYPE_NULL;

  /* clone parameters */
  XCLONE(ts->param, param, sizeof(*param));
  ts->param->iparam = NULL;

  /* initialize data */
  ts->data = thermostat_data_new(ts->param, NULL);

  return ts;
}


void thermostat_null_delete(thermostat_t *ts)
{
  free(ts->param);
  thermostat_data_delete(ts->data);
  free(ts);
}


real thermostat_null_apply(thermostat_t *ts)
{
  thermostat_param_t *tsp = ts->param;

  return mdutils_ekin(tsp->n, tsp->mass, tsp->v);
}

#endif /* THERMOSTAT_NULL_H__ */


