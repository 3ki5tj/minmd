#ifndef THERMOSTAT_H__
#define THERMOSTAT_H__

#include "thermostat_null.h"
#include "thermostat_vrescaling.h"

thermostat_t *thermostat_new(int type, thermostat_param_t *param)
{
  if (type == THERMOSTAT_TYPE_NULL) {
    return thermostat_null_new(param);
  } else if (type == THERMOSTAT_TYPE_VRESCALING) {
    return thermostat_vrescaling_new(param);
  } else {
    fprintf(stderr, "Error: thermostat_new() does not support thermostat type %d\n", type);
    return NULL;
  }
}


real thermostat_apply(thermostat_t *ts)
{
  if (ts->type == THERMOSTAT_TYPE_NULL) {
    return thermostat_null_apply(ts);
  } else if (ts->type == THERMOSTAT_TYPE_VRESCALING) {
    return thermostat_vrescaling_apply(ts);
  } else {
    fprintf(stderr, "Error: thermostat_apply() does not support thermostat type %d\n", ts->type);
    return 0;
  }
}


real thermostat_delete(thermostat_t *ts)
{
  if (ts->type == THERMOSTAT_TYPE_NULL) {
    thermostat_null_delete(ts);
  } else if (ts->type == THERMOSTAT_TYPE_VRESCALING) {
    thermostat_vrescaling_delete(ts);
  } else {
    fprintf(stderr, "Error: thermostat_delete() does not support thermostat type %d\n", ts->type);
  }
}

#endif /* THERMOSTAT_H__ */
