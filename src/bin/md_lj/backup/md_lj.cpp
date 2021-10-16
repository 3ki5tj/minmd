#include "LJSimul.hpp"

int n = 128;
double rho = 0.7;
double rc_def = 2.5;
long long int nequil = 10000;
long long int nsteps = 100000;
const char fn_checkpoint = "lj.checkpoint";

int main(int argc, char **argv)
{
  LJSimul lj_simul(n, rho, rc_def);

  lj_simul.do_md(nequil);
  lj_simul.do_md(nsteps);

  lj_simul.save_checkpoint(fn_checkpoint);
  lj_simul.print_summary();

  return 0;
}
