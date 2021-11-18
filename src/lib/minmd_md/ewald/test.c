#include "ewald.h"



/* test the direct Ewald method on the one-particle one-component plasma system
 * which is a single ion immersed in an oppositely-charged background */
void test_ewald_direct_ocp1(void)
{
  /* system specification */
  real charge[1] = {1.0};
  real x[1][DIM] = {{0.1, 0.2, 0.3}};
  real f[1][DIM] = {{0.0, 0.0, 0.0}};

  /* generic Ewald parameters */
  ewald_param_t ewp = {
    .n = 1,
    .box = {1.0, 1.0, 1.0},
    .sigma = 0.0, /* choose automatically */
    .tol = 1e-8,
    .kee = 1.0,
    .charge = charge,
    .x = x,
    .f = f,
  };
  /* create a new Ewald object */
  ewald_t *ew = ewald_new(EWALD_TYPE_DIRECT, &ewp);

  /* clear the force at the beginning of step */
  ewald_force_options_t ewf_opt = {
    .zero_forces = 1,
  };

  /* get the electrostatic energy from the Ewald method */
  double ene = ewald_force(ew, &ewf_opt);
  printf("energy %.8f/2=%g, real %g, recip %g, self %g, background %g\n\n\n",
      ene*2, ene, ew->data->real_energy, ew->data->recip_energy,
      ew->data->self_energy, ew->data->background_energy);
  
  /* delete the Ewald object */
  ewald_delete(ew);
}



/* test the direct Ewald method on the two-ion system
 * which consists of a pair of oppositely-charged ions
 *
 * Test if the force matches the energy */
void test_ewald_direct_two_ions(void)
{
  /* system specification */
  real charge[2] = {1.0, -1.0};
  real x[2][DIM] = {{0.1, 0.0, 0.0}, {0.5, 0.0, 0.0}};
  real f[2][DIM] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  /* generic Ewald parameters */
  ewald_param_t ewp = {
    .n = 2,
    .box = {1.0, 1.0, 1.0},
    .sigma = 0.0, /* let the algorithm choose automatically */
    .tol = 1e-7,
    .kee = 1.0,
    .charge = charge,
    .x = x,
    .f = f,
  };
  /* create a new Ewald object */
  ewald_t *ew = ewald_new(EWALD_TYPE_DIRECT, &ewp);

  /* since the electrostatic force is the only force for our system
   * we need to zero it at the beginning of every step */
  ewald_force_options_t ewf_opt = {
    .zero_forces = 1,
  };

  /* compute the initial energy */
  double ene1 = ewald_force(ew, &ewf_opt),
         ereal1 = ew->data->real_energy,
         erecip1 = ew->data->recip_energy;
  printf("energy %g, real %g, recip %g, self %g, background %g\n",
      ene1, ereal1, erecip1,
      ew->data->self_energy, ew->data->background_energy);
  printf("f:\n%8.5f %8.5f %8.5f\n%8.5f %8.5f %8.5f\ndist %g\n",
      f[0][0], f[0][1], f[0][2],
      f[1][0], f[1][1], f[1][2],
      x[1][0] - x[0][0]);

  /* move along the force so that the energy drops by delta */
  double f2 = vec_sqr(f[0]) + vec_sqr(f[1]),
         delta = 0.01,
         sdx = delta/f2;
  printf("f2 %g, delta %g, sdx %g, dx %g, %g\n", f2, delta, sdx, sdx*f[0][0], sdx*f[1][0]);
  vec_sinc(x[0], f[0], sdx);
  vec_sinc(x[1], f[1], sdx);
  printf("x[0] %g %g %g; x[1] %g %g %g\n", x[0][0], x[0][1], x[0][2], x[1][0], x[1][1], x[1][2]);

  /* test the new energy */
  double ene2 = ewald_force(ew, &ewf_opt),
         ereal2 = ew->data->real_energy,
         erecip2 = ew->data->recip_energy;
  printf("energy %g, real %g, recip %g, self %g, background %g\n",
      ene2, ereal2, erecip2,
      ew->data->self_energy, ew->data->background_energy);
  printf("de/delta %g, dreal/delta %g, drecip/delta %g\n\n\n",
      (ene2 - ene1)/delta,
      (ereal2 - ereal1)/delta,
      (erecip2 - erecip1)/delta);

  /* delete the Ewald object */
  ewald_delete(ew);
}


int main(void)
{
  /* test on the one-particle one-component plasma system */
  test_ewald_direct_ocp1();

  /* test on a system consisting of two oppositely-charged ions */
  test_ewald_direct_two_ions();

  return 0;
}


