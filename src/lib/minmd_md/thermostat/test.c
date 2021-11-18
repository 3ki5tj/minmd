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
  /* initialize parameters specific for the velocity-rescaling thermostat
   * velocity-rescaling thermostat requires a random number generator
   * so we create one here */
  rng_t *rng = rng_new(0, time(NULL));
  thermostat_vrescaling_param_t vr_param = {
    .rng = rng
  };

  /* initialize parameters for a generic thermostat */
  thermostat_param_t param = {
    .n = n,
    .n_dof = n*DIM - DIM,
    .tp = 2.0,
    .dt = 0.05,
    .boltz = 1.0,
    .mass = mass,
    .v = v,
    .iparam = &vr_param
  };

  /* create a new instance of the thermostat */
  thermostat_t *ts = thermostat_new(THERMOSTAT_TYPE_VRESCALING, &param);

  /* initialize random velocities */
  init_velocities(n, mass, v, rng);
 
  /* apply the thermostat several times */
  int k;
  for (k = 0; k < 50; k++) {
    real ek = thermostat_apply(ts);
    printf("k %d, ek %g vs %g\n", k, ek, 0.5 * param.n_dof * param.tp);
  }

  /* delete the thermostat */
  thermostat_delete(ts);

  rng_delete(rng);
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

  /* test the velocity-rescaling thermostat */
  test_vrescaling(n, mass, v);

  /* free memory */
  free(v);
  free(mass);

  return 0;
}


