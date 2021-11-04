#include "constraint.h"


constraint_t *init_shake(int n, real *mass,
    real (*x0)[DIM], real (*x1)[DIM], real (*v)[DIM])
{
  int i;

  /* 1. initialize the constraint pair list */
  constraint_pair_param_t *prparam;
  XNEW(prparam, n - 1);
  for (i = 0; i < n - 1; i++) {
    prparam[i].i = i;
    prparam[i].j = i + 1;
    prparam[i].dist = 1.0;
  }
  constraint_pairlist_t pair_list = {
    .n = n - 1,
    .pr = prparam,
  };

  /* 2. SHAKE parameters and generic constraint solver parameters */
  /* 2.1. SHAKE-specific parameters */
  constraint_shake_param_t shake_param = {
    .citmax = 100,
    .ctol = 1e-8,
    .vitmax = 1000,
    .vtol = 1e-8,
    .prlist = &pair_list,
  };
  /* 2.2. generic constraint solver parameters */
  constraint_param_t csp = {
    .n = n,
    .mass = mass,
    .x0 = x0,
    .x1 = x1,
    .v = v,
    .iparam = &shake_param,
  };

  /* 3. create a new constraint solver instance */
  constraint_t *cs = constraint_new(CONSTRAINT_TYPE_SHAKE, &csp);

  /* 4. ok to free the input pair list after the constraint solver is initialized */
  free(prparam);

  return cs;
}


int main(void)
{
  int n = 10, i;
  real *mass;
  real (*x0)[DIM];
  real (*x1)[DIM];
  real (*v)[DIM];
  rng_t *rng = rng_new(0, time(NULL));

  /* 1. initialize masses, coordinates and velocities */
  XNEW(mass, n);
  XNEW(x0, n);
  XNEW(x1, n);
  XNEW(v, n);
  for (i = 0; i < n; i++) {
    mass[i] = 1 + rng_gauss(rng)*0.01;
    x0[i][0] = i*1.0;
    x0[i][1] = 0.0;
    x0[i][2] = 0.0;
    x1[i][0] = i*1.0;
    x1[i][1] = 0.0;
    x1[i][2] = 0.0;
    vec_rand_gauss(v[i], 1.0, rng);
    vec_sinc(x1[i], v[i], 0.02);
    vec_rand_gauss(v[i], 1.0, rng);
  }

  /* 2. initialize constraint solver */
  constraint_t *cs = init_shake(n, mass, x0, x1, v);

  /* 3. display the distance before applying the constraint solver */
  printf("BEFORE:\n");
  real dx[DIM], dv[DIM], dist;
  for (i = 0; i < n - 1; i++) {
    dist = vec_distx(dx, x1[i], x1[i+1]);
    vec_diff(dv, v[i], v[i+1]);
    printf("%d: dist %g, v.d = %g\n", i, dist, vec_dot(dv, dx));
  }
  printf("\n\n");

  /* 4. apply the constraint solver */
  constraint_apply(cs, CONSTRAINT_APPLY_COORDINATES|CONSTRAINT_APPLY_VELOCITIES);

  /* 5. display the distances after applying the constraint solver */
  printf("AFTER:\n");
  for (i = 0; i < n - 1; i++) {
    dist = vec_distx(dx, x1[i], x1[i+1]);
    vec_diff(dv, v[i], v[i+1]);
    printf("%d: dist %g, v.d = %g\n", i, dist, vec_dot(dv, dx));
  }
  printf("\n\n");

  constraint_delete(cs);

  free(x0);
  free(x1);
  free(v);
  free(mass);
  rng_delete(rng);
  return 0;
}
