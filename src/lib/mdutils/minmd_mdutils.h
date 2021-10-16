#ifndef MDUTILS_H__
#define MDUTILS_H__

#include "def.h"
#include "rng.h"
#include "minmd_utils.h"

#ifdef MINMD_2D

void mdutils_init_face_centered_lattice(int n, real side_length, real (*x)[DIM], rng_t *r)
{
  int i, j, id, nside;
  real a, noise;

  nside = (int) (pow(2*n, 1.0/DIM) + .999999); /* # of particles per side */
  a = side_length / nside;
  noise = a * 1e-5;
  for (id = 0, i = 0; i < n_side && id < n; i++) {
    for (j = 0; j < n_side && id < n; j++) {
      if ((i+j) % 2 != 0) continue;
      /* add some noise to prevent two atoms happened to
       * be separated by precisely some special cutoff distance,
       * which might be half of the box */
      x[id][0] = (i + .5) * a + noise * (2*rng_rand01(r) - 1);
      x[id][1] = (j + .5) * a + noise * (2*rng_rand01(r) - 1);
      id++;
    }
  }
}


#else

void mdutils_init_face_centered_lattice(int n, real side_length, real (*x)[DIM], rng_t *r)
{
  int i, j, k, id, nside;
  real a, noise;

  nside = (int) (pow(2*n, 1.0/DIM) + .999999); /* # of particles per side */
  a = side_length / nside;
  noise = a * 1e-5;
  for (id = 0, i = 0; i < nside && id < n; i++) {
    for (j = 0; j < nside && id < n; j++) {
      for (k = 0; k < nside && id < n; k++) {
        if ((i+j+k) % 2 != 0) {
          continue;
        }
        /* add some noise to prevent two atoms happened to
         * be separated by precisely some special cutoff distance,
         * which might be half of the box */
        x[id][0] = (real) ((i + .5) * a + noise * (2*rng_rand01(r) - 1));
        x[id][1] = (real) ((j + .5) * a + noise * (2*rng_rand01(r) - 1));
        x[id][2] = (real) ((k + .5) * a + noise * (2*rng_rand01(r) - 1));
        id++;
      }
    }
  }
}


#endif /* MINMD_2D */


/* remove the center of mass motion */
void mdutils_remove_com(int n, const real *mass, real (*v)[DIM])
{
  int i, d;
  real momentum[DIM] = {0}, total_mass = 0;

  for (i = 0; i < n; i++) {
    for (d = 0; d < DIM; d++) {
      momentum[d] += mass[i] * v[i][d];
    }
    total_mass += mass[i];
  }

  for (d = 0; d < DIM; d++) {
    momentum[d] /= total_mass;
    for (i = 0; i < n; i++) {
      v[i][d] -= momentum[d];
    }
  }
}


void mdutils_init_velocities(int n, const real *mass, real (*v)[DIM], real tp_ref, rng_t *r)
{
  int i, d;

  for (i = 0; i < n; i++) {
    real mag = (real) sqrt(tp_ref/mass[i]);
    for (d = 0; d < DIM; d++) {
      v[i][d] = mag * rng_gauss(r);
    }
  }

  mdutils_remove_com(n, mass, v);
}


#endif /* MDUTILS_H__ */
