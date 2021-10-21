#include "lj.h"
#include "ljeos.h"

md_lj_param_t param = {
  .n = 108,
  .rho = 0.7, /* reduced density */
  .tp = 1.5, /* reduced temperature */
  .rc_def = 2.5,
  .md_dt = 0.002, /* MD time step */
  //.thermostat_type = THERMOSTAT_TYPE_NULL, /* no thermostat, good for testing energy conservation */
  .thermostat_type = THERMOSTAT_TYPE_VRESCALING, /* velocity rescaling thermostat */
  .vr_dt = 0.01, /* effective time step for velocity rescaling thermostat */
};

long long nequil = 100;
long long nsteps = 100;


int main(void)
{
  md_lj_t *lj;

  lj = md_lj_open(&param);

  /* equilibration run */
  md_running_param_t equil_param = {
    .nsteps = nequil,
    .verbose = 2,
    .do_stat = 0,
  };
  md_lj_run(lj, &equil_param);

  /* production run */
  md_running_param_t prod_param = {
    .nsteps = nsteps,
    .verbose = 0,
    .do_stat = 1,
  };
  md_lj_run(lj, &prod_param);
  // md_lj_print_summary(lj, );

  md_lj_free(lj);

  return 0;
}
