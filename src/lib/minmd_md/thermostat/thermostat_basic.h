#ifndef THERMOSTAT_BASIC_H__
#define THERMOSTAT_BASIC_H__

#include "utils.h"

/* common thermostat parameters */
typedef struct {
  int n;
  int n_dof; /* number of degrees of freedom */
  double tp; /* target temperature */
  double dt; /* effective time step for each operation */
  double boltz;
  const real *mass;
  real (*v)[DIM];
  void *iparam; /* algorithm-specific parameters */
} thermostat_param_t;

enum {
  THERMOSTAT_TYPE_NULL = 0,
  THERMOSTAT_TYPE_VRESCALING = 1,
  /* THERMOSTAT_TYPE_NHCHAIN, */
  THERMOSTAT_TYPE_COUNT
};

typedef struct {
  void *idata;
} thermostat_data_t;

typedef struct {
  int type;
  thermostat_param_t *param;
  thermostat_data_t *data;
} thermostat_t;


INLINE thermostat_data_t *thermostat_data_new(const thermostat_param_t *tsp,
    void *idata)
{
  thermostat_data_t *tsd;

  XNEW(tsd, 1);
  tsd->idata = idata;
}


INLINE void thermostat_data_delete(thermostat_data_t *tsd)
{
  free(tsd);
}


#endif /* THERMOSTAT_BASIC_H__ */

