/* force consistency test */
#include "ocp.h"


md_ocp_param_t param = {
  .n = 108,
  .rho = 0.7, /* reduced density */
  .tp = 2.0, /* reduced temperature */
  .rc_def = 2.5,
  .ewald_tol = 1e-5,
};


double get_f2(md_ocp_t *ocp)
{
  int i, n = ocp->n;
  double f2 = 0;

  for (i = 0; i < n; i++) {
    f2 += vec_sqr(ocp->f[i]);
  }
  return f2;
}

void force_test(md_ocp_t *ocp)
{
  int j, i, n = ocp->n;
  double f2, ene1, ene2;

  // perturbing the initial lattice structure to add force
  for (j = 0; j < 1000; j++) {
    ene1 = md_ocp_force(ocp);
    f2 = get_f2(ocp);
    printf("energy %g, f2 %g\n", ene1, f2);
    if (f2 > 100.0*n) {
      printf("perturbation done!\n\n\n");
      break;
    }
    // randomize the system a bit;
    for (i = 0; i < n; i++) {
      real dx[DIM];
      vec_rand_gauss(dx, 0.05, ocp->rng);
      vec_sinc(ocp->x[i], dx, 1.0);
    }
  }

  ene1 = md_ocp_force(ocp);

  /* increment the coordinates along the force by an amount
   * such that the energy increases by delta */
  double delta = 0.01, fdel;
  f2 = get_f2(ocp);
  fdel = delta/f2;
  for (i = 0; i < n; i++) {
    vec_sinc(ocp->x[i], ocp->f[i], fdel);
  }
        
  ene2 = md_ocp_force(ocp);

  printf("ene1 %g, ene2 %g (%g), de/delta %g, f2 %g, fdel %g\n",
      ene1, ene2, ocp->epot_data->epot, (ene1 - ene2)/delta, f2, fdel );
}


int main(void)
{
  md_ocp_t *ocp = md_ocp_new(&param);

  force_test(ocp);

  md_ocp_delete(ocp);

  return 0;
}
