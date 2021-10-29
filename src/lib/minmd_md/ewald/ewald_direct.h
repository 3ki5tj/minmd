#ifndef EWALD_DIRECT_H__
#define EWALD_DIRECT_H__

#include "ewald_basic.h"


ewald_t *ewald_direct_init(ewald_param_t *param)
{
  ewald_t *ew;
  ewald_data_t *ew_data;

  XNEW(ew, 1);
  ew->type = EWALD_TYPE_DIRECT;

  /* clone parameters */
  XCLONE(ew->param, param, sizeof(*param));
  ew->param->algo_param = NULL;

  /* build data */
  ew->data = ewald_data_init(ew->param);
  ew->data->algo_data = NULL;
  return ew;
}


INLINE void ewald_direct_free(ewald_t *ew)
{
  if (ew->param->algo_param != NULL) {
    free(ew->param->algo_param);
  }
  free(ew->param);

  if (ew->data->algo_data != NULL) {
    free(ew->data->algo_data);
  }
  free(ew->data);
  free(ew);
}


INLINE real ewald_direct_force_real(ewald_param_t *ewp, ewald_data_t *ewd, double *virial)
{
  int n = ewp->n, i, j;
  double rc = ewd->rc,
         alpha = ewd->alpha,
         sqrta = ewd->sqrt_alpha,
         sqrtapi = ewd->sqrt_alpha_over_pi,
         kee = ewp->kee;
  double fscl, ene = 0;
  real (*x)[DIM] = ewp->x, dx[DIM], rij,
       (*f)[DIM] = ewp->f, qi, qj;

  for (i = 0; i < n; i++) {
    qi = ewp->charge[i];
    for (j = i + 1; j < n; j++) {
      qj = ewp->charge[j];
      rij = ewald_dist_func(dx, x[i], x[j], ewp->box);
      if (rij < rc) {
        double kqij = kee*qi*qj,
               ec = erfc(sqrta*rij),
               rij2 = rij*rij,
               eij = kqij*ec/rij;
        ene += eij;
        fscl = (real)((kqij*2*sqrtapi*exp(-alpha*rij2) + eij)/rij2);
        vec_sinc(f[i], dx,  fscl);
        vec_sinc(f[j], dx, -fscl);
      }
    }
  }
  /* in this system, virial and energy are the same */
  *virial = ene; /* fi.ri */

  return ene;
}


INLINE real ewald_direct_force_recip(ewald_param_t *ewp, ewald_data_t *ewd, double *virial)
{
  int ik[DIM], *ikmax = ewd->ikmax;
  real kv[DIM], *kbox = ewd->kbox;
  double k2, ene = 0;
  int i, n = ewp->n;
  const real *q = ewp->charge;
  real (*x)[DIM] = ewp->x,
       (*f)[DIM] = ewp->f;
  double efac = ewp->kee*2*PI/ewd->vol;
  double *ck, *sk, cksum, sksum, c, s, kx;

  XNEW(ck, n);
  XNEW(sk, n);

  for (ik[0] = -ikmax[0]; ik[0] <= ikmax[0]; ik[0]++) {
    kv[0] = ik[0] * ewd->kbox[0];
    for (ik[1] = -ikmax[1]; ik[1] <= ikmax[1]; ik[1]++) {
      kv[1] = ik[1] * ewd->kbox[1];
      for (ik[2] = -ikmax[2]; ik[2] <= ikmax[2]; ik[2]++) {
        kv[2] = ik[2] * ewd->kbox[2];
        if (ik[0] == 0 && ik[1] == 0 && ik[2] == 0) {
          continue;
        }
        cksum = 0;
        sksum = 0;
        for (i = 0; i < n; i++) {
          double qi = q[i];
          kx = vec_dot(kv, x[i]);
          c = qi*cos(kx);
          cksum += c;
          ck[i] = c;
          s = qi*sin(kx);
          sksum += s;
          sk[i] = s;
        }
        k2 = vec_sqr(kv);
        double fac = efac*exp(-k2*ewd->one_over_4alpha)/k2;
        ene += (cksum*cksum + sksum*sksum)*fac;
        /* TODO: compute virial */
        for (i = 0; i < n; i++) {
          real fscl = (real) ((cksum*sk[i] - sksum*ck[i])*2*fac);
          vec_sinc(f[i], kv, fscl);
        }
      }
    }
  }
  free(ck);
  free(sk);
  return ene;
}


real ewald_direct_force(ewald_t *ew, ewald_force_options_t *ewf_opt)
{
  ewald_param_t *ewp = ew->param;
  ewald_data_t *ewd = ew->data;

  if (ewf_opt->zero_forces) {
    int i;
    for (i = 0; i < ewp->n; i++) {
      vec_zero(ewp->f[i]);
    }
  }

  ewd->real_energy = ewald_direct_force_real(ewp, ewd, &ewd->real_virial);
  ewd->recip_energy = ewald_direct_force_recip(ewp, ewd, &ewd->recip_virial);
  //printf("real %g, recip %g\n", ewd->real_energy, ewd->recip_energy);
  ewd->energy = ewd->real_energy
              + ewd->recip_energy
              + ewd->self_energy
              + ewd->background_energy;
  ewd->virial = ewd->real_virial
              + ewd->recip_virial;

  return ewd->energy;
}


#endif /* EWALD_DIRECT_H__ */
