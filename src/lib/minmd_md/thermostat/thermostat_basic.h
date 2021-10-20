#ifndef THERMOSTAT_BASIC_H__
#define THERMOSTAT_BASIC_H__


/* common thermostat parameters */
typedef struct {
  int n;
  int n_dof;
  double tp; /* target temperature */
  double dt; /* effective time step for each operation */
  const real *mass;
  real (*v)[DIM];
  void *algo_param; /* algorithm-specific parameters */
} thermostat_param_t;

enum {
  THERMOSTAT_TYPE_DEFAULT = 0,
  THERMOSTAT_TYPE_VRESCALING = 0,
  /* THERMOSTAT_TYPE_NHCHAIN, */
  THERMOSTAT_TYPE_COUNT
};

typedef struct {
  int type;
  thermostat_param_t *param;
  void *data;
} thermostat_t;


INLINE void thermostat_free(thermostat_t *ts)
{
  if (ts->data != NULL) {
    free(ts->data);
  }
  free(ts);
}


#endif /* THERMOSTAT_BASIC_H__ */

