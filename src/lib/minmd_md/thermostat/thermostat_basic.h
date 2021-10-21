#ifndef THERMOSTAT_BASIC_H__
#define THERMOSTAT_BASIC_H__

#include "def.h"

/* common thermostat parameters */
typedef struct {
  int n;
  int n_dof; /* number of degrees of freedom */
  double tp; /* target temperature */
  double dt; /* effective time step for each operation */
  double boltz;
  const real *mass;
  real (*v)[DIM];
  void *algo_param; /* algorithm-specific parameters */
} thermostat_param_t;

enum {
  THERMOSTAT_TYPE_NULL = 0,
  THERMOSTAT_TYPE_VRESCALING = 1,
  /* THERMOSTAT_TYPE_NHCHAIN, */
  THERMOSTAT_TYPE_COUNT
};

typedef struct {
  int type;
  thermostat_param_t *param;
  void *data;
} thermostat_t;


#endif /* THERMOSTAT_BASIC_H__ */

