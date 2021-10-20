#ifndef VRESCALING_H__
#define VRESCALING_H__

#include "utils.h"
#include "rng.h"
#include "mdutils.h"

typedef struct {
  int n;
  int n_dof;
  double dt;
} vrescaling_t;



INLINE vrescaling_t *vrescaling_init(int n, int n_dof, double dt)
{
  vrescaling_t *vr;

  XNEW(vr, 1);
  vr->n = n;
  vr->n_dof = n_dof;
  vr->dt = dt;
  return vr;
}


INLINE void vrescaling_free(vrescaling_t *vr)
{
  free(vr);
}


INLINE real vrescaling_apply(vrescaling_t *vr,
    real (*v)[DIM], const real *mass, real tp, rng_t *rng)
{
  int i, n = vr->n;
  double dt = vr->dt;
  double ek1, ek2, s, c, r, r2;

  c = exp(-dt);
  ek1 = mdutils_ekin(n, mass, v);
  r = rng_gauss(rng);
  r2 = rng_chisqr(rng, vr->n_dof - 1);
  ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp * 0.5 - ek1)
      + 2 * r * sqrt(c * (1 - c) * ek1 * tp * 0.5);
  if (ek2 < 0) {
    ek2 = 0;
  }
  s = sqrt(ek2 / ek1);
  for (i = 0; i < n; i++) {
    vec_smul(v[i], s);
  }
  return ek2;
}

#endif /* VRESCALING_H__ */


