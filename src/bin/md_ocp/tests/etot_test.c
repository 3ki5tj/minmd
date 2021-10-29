/* this program is used to test the energy conservation */

#include "ocp.h"


md_ocp_param_t param = {
  .n = 108,
  .rho = 0.7, /* reduced density */
  .tp = 2.0, /* reduced temperature */
  .rc_def = 2.5,
  .md_dt = 0.005, /* MD time step */
  .thermostat_type = THERMOSTAT_TYPE_VRESCALING,
  .vr_dt = 0.02, /* effective time step for velocity rescaling thermostat */
  .ewald_tol = 1e-4,
};

md_running_param_t equil_param = {
  .nsteps = 200,
  .verbose = 2,
  .nst_print = 10,
  .do_stat = 0,
};

int main(void)
{
  md_ocp_t *ocp = md_ocp_init(&param);
  md_ocp_run(ocp, &equil_param);
 
  /* here we replace the thermostat by switching to type THERMOSTAT_TYPE_NULL */
  thermostat_t *old_thermostat = ocp->thermostat;
  ocp->thermostat = thermostat_init(THERMOSTAT_TYPE_NULL, old_thermostat->param);
  thermostat_free(old_thermostat);

  printf("\n\nStarting testing energy conservation...\n");
  printf("Pay attention to etot_shifted\n");
  md_ocp_run(ocp, &equil_param);
  md_ocp_free(ocp);

  return 0;
}
