#include "ocp.h"


md_ocp_param_t param = {
  .n = 108,
  .rho = 0.7, /* reduced density */
  .tp = 2.0, /* reduced temperature */
  .rc_def = 2.5,
  .md_dt = 0.005, /* MD time step */
  //.thermostat_type = THERMOSTAT_TYPE_NULL, /* no thermostat, good for testing energy conservation */
  .thermostat_type = THERMOSTAT_TYPE_VRESCALING, /* velocity rescaling thermostat */
  .vr_dt = 0.02, /* effective time step for velocity rescaling thermostat */
};

/* equilibration run parameters */
md_running_param_t equil_param = {
  .nsteps = 500,
  .verbose = 2,
  .nst_print = 100,
  .do_stat = 0,
};

/* production run parameters */
md_running_param_t prod_param = {
  .nsteps = 2000,
  .verbose = 2,
  .nst_print = 100,
  .do_stat = 1,
};



int main(void)
{
  md_ocp_t *ocp = md_ocp_init(&param);

  /* equilibration run */
  md_ocp_run(ocp, &equil_param);

  /* production run */
  printf("\n\nStarting production run...\n");
  md_ocp_run(ocp, &prod_param);
  md_ocp_print_stat(ocp);

  md_ocp_free(ocp);

  return 0;
}
