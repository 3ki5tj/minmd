/* force consistency test */
#include "lj.h"


md_lj_param_t param = {
  .n = 108,
  .rho = 0.7, /* reduced density */
  .tp = 2.0, /* reduced temperature */
  .rc_def = 2.5,
};


double get_f2(md_lj_t *lj)
{
  int i, n = lj->n;
  double f2 = 0;

  for (i = 0; i < n; i++) {
    f2 += vec_sqr(lj->f[i]);
  }
  return f2;
}

void force_test(md_lj_t *lj)
{
  int j, i, n = lj->n;
  double f2, ene1, ene2;

  // perturbing the initial lattice structure to add force
  for (j = 0; j < 1000; j++) {
    ene1 = md_lj_force(lj);
    f2 = get_f2(lj);
    printf("perturbation %d: energy %g, f2 %g\n", j, ene1, f2);
    if (f2 > 20.0*n) {
      printf("perturbation done!\n\n");
      break;
    }
    // randomize the system a bit;
    for (i = 0; i < n; i++) {
      real dx[DIM];
      vec_rand_gauss(dx, 0.05, lj->rng);
      vec_sinc(lj->x[i], dx, 1.0);
    }
  }

  ene1 = md_lj_force(lj);

  /* increment the coordinates along the force by an amount
   * such that the energy increases by delta */
  double delta = 0.01, fdel;
  f2 = get_f2(lj);
  fdel = delta/f2;
  for (i = 0; i < n; i++) {
    vec_sinc(lj->x[i], lj->f[i], fdel);
  }
        
  ene2 = md_lj_force(lj);

  printf("ene1 %g, ene2 %g, de/delta %g, f2 %g, fdel %g\n",
      ene1, ene2, (ene1 - ene2)/delta, f2, fdel);
}


int main(void)
{
  md_lj_t *lj = md_lj_init(&param);

  force_test(lj);

  md_lj_free(lj);

  return 0;
}
