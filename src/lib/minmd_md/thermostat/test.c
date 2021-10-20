#include "thermostat.h"



void init_velocities(int n, const real *mass, real (*v)[DIM], rng_t *rng)
{
  int i;
  for (i = 0; i < n; i++) {
    vec_rand_gauss(v[i], 1.0, rng);
  }
  mdutils_remove_com(n, mass, v);
}



void test_vrescaling(int n, const real *mass, real (*v)[DIM])
{
  thermostat_param_t param = {0};
  param.n = n;
  param.n_dof = n*DIM - DIM;
  param.tp = 2.0;
  param.dt = 0.05;
  param.mass = mass;
  param.v = v;

  /* velocity rescaling parameters */
  vrescaling_param_t vr_param = {0};
  rng_t *rng = rng_init(0, time(NULL));
  vr_param.rng = rng;
  param.algo_param = &vr_param;

  thermostat_t *ts = thermostat_init(THERMOSTAT_TYPE_VRESCALING, &param);

  init_velocities(n, mass, v, rng);
 
  int k;
  for (k = 0; k < 50; k++) {
    real ek = thermostat_apply(ts);
    printf("k %d, ek %g vs %g\n", k, ek, 0.5 * param.n_dof * param.tp);
  }

  thermostat_free(ts);
  rng_free(rng);
}


int main(void)
{
  int n = 10, i;
  real *mass;
  real (*v)[DIM];

  /* allocate memory for masses and velocities */
  XNEW(mass, n);
  XNEW(v, n);
  for (i = 0; i < n; i++) {
    mass[i] = 1;
  }

  test_vrescaling(n, mass, v);

  free(v);
  free(mass);
  return 0;
}
