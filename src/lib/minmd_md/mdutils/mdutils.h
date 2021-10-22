#ifndef MDUTILS_H__
#define MDUTILS_H__

#include "def.h"
#include "rng.h"
#include "utils.h"
#include "vec.h"



/* remove the center of mass motion */
void mdutils_remove_com(int n, const real *mass, real (*v)[DIM])
{
  int i;
  real momentum[DIM], total_mass = 0;

  vec_zero(momentum);
  for (i = 0; i < n; i++) {
    vec_sinc(momentum, v[i], mass[i]);
    total_mass += mass[i];
  }

  vec_smul(momentum, 1.0/total_mass);

  for (i = 0; i < n; i++) {
    vec_dec(v[i], momentum);
  }
}


void mdutils_init_velocities(int n, const real *mass, real (*v)[DIM], real tp_ref, rng_t *rng)
{
  int i, d;

  for (i = 0; i < n; i++) {
    real mag = (real) sqrt(tp_ref/mass[i]);
    vec_rand_gauss(v[i], mag, rng);
  }

  mdutils_remove_com(n, mass, v);
}


INLINE real mdutils_ekin(int n, const real *mass, real (*v)[DIM])
{
  int i;
  double ek = 0;

  if (mass == NULL) {
    for (i = 0; i < n; i++) {
      ek += vec_sqr(v[i]);
    }
  } else {
    for (i = 0; i < n; i++) {
      ek += mass[i] * vec_sqr(v[i]);
    }
  }
  return (real) (ek * 0.5);
}




#endif /* MDUTILS_H__ */

