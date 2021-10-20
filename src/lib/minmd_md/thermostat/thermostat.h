#ifndef THERMOSTAT_H__
#define THERMOSTAT_H__

#include "vrescaling.h"

thermostat_t *thermostat_init(int type, thermostat_param_t *param)
{
  if (type == THERMOSTAT_TYPE_VRESCALING) {
    return vrescaling_init(param);
  } else {
    fprintf(stderr, "Error: thermostat_init() does not support thermostat type %d\n", type);
    return NULL;
  }
}


real thermostat_apply(thermostat_t *ts)
{
  if (ts->type == THERMOSTAT_TYPE_VRESCALING) {
    return vrescaling_apply(ts);
  } else {
    fprintf(stderr, "Error: thermostat_apply() does not support thermostat type %d\n", ts->type);
    return 0;
  }
}

#endif /* THERMOSTAT_H__ */
