#include "lj.h"


md_lj_param_t param = {
  .n = 108,
  .rho = 0.7, /* reduced density */
  .tp = 2.0, /* reduced temperature */
  .rc_def = 2.5,
  .md_dt = 0.002, /* MD time step */
  //.thermostat_type = THERMOSTAT_TYPE_NULL, /* no thermostat, good for testing energy conservation */
  .thermostat_type = THERMOSTAT_TYPE_VRESCALING, /* velocity rescaling thermostat */
  .vr_dt = 0.01, /* effective time step for velocity rescaling thermostat */
};

/* equilibration run parameters */
md_running_param_t equil_param = {
  .nsteps = 10000,
  .verbose = 2,
  .nst_print = 1000,
  .do_stat = 0,
};

/* production run parameters */
md_running_param_t prod_param = {
  .nsteps = 10000,
  .verbose = 0,
  .do_stat = 1,
};



int main(void)
{
  md_lj_t *lj = md_lj_open(&param);

  /* equilibration run */
  md_lj_run(lj, &equil_param);

  /* production run */
  md_lj_run(lj, &prod_param);
  md_lj_print_stat(lj);

  md_lj_free(lj);

  return 0;
}
