#ifndef THERMOSTAT_H__
#define THERMOSTAT_H__

#include "thermostat_null.h"
#include "thermostat_vrescaling.h"

thermostat_t *thermostat_init(int type, thermostat_param_t *param)
{
  if (type == THERMOSTAT_TYPE_NULL) {
    return thermostat_null_init(param);
  } else if (type == THERMOSTAT_TYPE_VRESCALING) {
    return thermostat_vrescaling_init(param);
  } else {
    fprintf(stderr, "Error: thermostat_init() does not support thermostat type %d\n", type);
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


real thermostat_free(thermostat_t *ts)
{
  if (ts->type == THERMOSTAT_TYPE_NULL) {
    thermostat_null_free(ts);
  } else if (ts->type == THERMOSTAT_TYPE_VRESCALING) {
    thermostat_vrescaling_free(ts);
  } else {
    fprintf(stderr, "Error: thermostat_free() does not support thermostat type %d\n", ts->type);
  }
}

#endif /* THERMOSTAT_H__ */
