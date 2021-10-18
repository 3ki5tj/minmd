#include "lj.h"
#include "ljeos.h"

int n = 108;
real rho = 0.7; /* reduced density */
real tp = 1.5; /* reduced temperature */
real rc_def = 2.5;
real md_dt = 0.002; /* MD time step */

long long nequil = 10000;
long long nsteps = 100000;


int main(void)
{
  lj_md_t *lj;


  lj = lj_md_open(n, rho, rc_def);

  lj_md_run(lj, nequil);
  lj_md_run(lj, nsteps);
  // lj_md_print_summary(lj, );

  return 0;
}
