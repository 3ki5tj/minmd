#ifndef EWALD_BASIC_H__
#define EWALD_BASIC_H__


#include "def.h"
#include "utils.h"
#include "vec.h"


/* commone ewald parameters */
typedef struct {
  int n;
  real box[DIM];
  real sigma; /* width of the smeared/cloud charge in the real space */
  real rc; /* user-specified real-space cutoff distance */
  real tol; /* tolerance for truncating the reciprocal space sum */
  real kee; /* e^2/(4 pi epsilon_0) */
  const real *charge;
  real (*x)[DIM];
  real (*f)[DIM];
  void *algo_param;
} ewald_param_t;


typedef struct {
  /* NOTE: we mainly use double-precision floating point numbers for internal data */
  double alpha; /* 1/(2 * sigma^2) */
  double sqrt_alpha;
  double sqrt_alpha_over_pi;
  double one_over_4alpha;
  double qtot;
  double q2tot;
  double vol;
  real kbox[DIM];
  int ikmax[DIM];

  double real_energy;
  double recip_energy;
  double self_energy;
  double background_energy;
  double energy;

  double real_virial;
  double recip_virial;
  double virial;

  /* relative error due to cutoff */
  double rc; /* adjusted real-space cutoff distance */
  double real_space_error;
  double recip_space_error;

  void *algo_data;
} ewald_data_t;


enum {
  EWALD_TYPE_NULL = 0,
  EWALD_TYPE_DIRECT = 1,
  EWALD_TYPE_COUNT
};

typedef struct {
  int type;
  ewald_param_t *param;
  ewald_data_t *data;
} ewald_t;


typedef struct {
  int zero_forces;
} ewald_force_options_t;



INLINE real ewald_dist_func(real *dx, const real *xi, const real *xj, const real *box)
{
  int d;
  for (d = 0; d < DIM; d++) {
    dx[d] = xi[d] - xj[d];
    dx[d] = fmod(dx[d] + 100.5*box[d], box[d]) - 0.5*box[d];
  }
  return vec_norm(dx);
}


INLINE real ewald_choose_rc(const ewald_param_t *ewp)
{
  double rc = ewp->rc, hbox;
  int d;

  if (rc <= 0) {
    rc = 1e30;
    for (d = 0; d < DIM; d++) {
      hbox = ewp->box[d]*0.5;
      if (hbox < rc) {
        rc = hbox;
      }
    }
  }
  return rc;
}


INLINE real ewald_choose_alpha(const ewald_param_t *ewp, double rc, double tol)
{
  double low = 0.5,
         high = 9.2,
         x, y;
  int i;

  /* find x such that erfc(x) >= tol */
  for (i = 0; i < 300; i++) {
    x = sqrt(low*high);
    y = erfc(x);
    //printf("round %d, trying x %g, y %g, low %g, high %g\n", i, x, y, low, high);
    if (y > tol) {
      low = x;
    } else {
      high = x;
    }
    if (high/low - 1 < 0.001) {
      break;
    }
  }

  x = high;
  /* sqrt(alpha)*rc = x; */
  return (x*x)/(rc*rc);
}


INLINE real ewald_total_charge(const ewald_param_t *ewp)
{
  int i;
  real qtot = 0;

  for (i = 0; i < ewp->n; i++) {
    qtot += ewp->charge[i];
  }
  return qtot;
}


INLINE real ewald_total_charge2(const ewald_param_t *ewp)
{
  int i;
  real q2tot = 0;

  //printf("%p, charge %p\n", ewp, ewp->charge);
  for (i = 0; i < ewp->n; i++) {
    q2tot += ewp->charge[i] * ewp->charge[i];
  }
  return q2tot;
}


INLINE real ewald_box_volume(const ewald_param_t *ewp)
{
  const real *box = ewp->box;
  return box[0] * box[1] * box[2];
}


INLINE real ewald_self_energy(const ewald_param_t *ewp, double alpha, double q2tot)
{
  return -ewp->kee * sqrt(alpha/PI) * q2tot;
}

INLINE real ewald_background_energy(const ewald_param_t *ewp, double alpha, double qtot, double vol)
{
  return -ewp->kee * qtot * qtot * PI * 0.5 / (alpha * vol);
}


/* real-space relative error relative to Ke^2/rc */
INLINE double ewald_real_space_error(double rc, double alpha)
{
  return erfc(sqrt(alpha)*rc);
}


/* reciprocal-space relative error relative to Ke^2/rc */
INLINE double ewald_recip_space_error(double k, double rc, double alpha, double vol)
{
  double kk = k * k;
  return 2*PI*exp(-kk/(4*alpha))/kk/vol*rc;
}


INLINE void ewald_recip_lattice(real *kbox, const ewald_param_t *ewp)
{
  int d;

  for (d = 0; d < DIM; d++) {
    kbox[d] = 2*PI/ewp->box[d];
  }
}


INLINE double ewald_ikmax(int *ikmax, double tol, const real *kbox, double rc, double alpha, double vol)
{
  int d, km;
  double errmax = 0; /* maximal cutoff error among three dimensions */

  for (d = 0; d < DIM; d++) {
    for (km = 1; km < 10000; km++) {
      double k = kbox[d] * km,
             err = ewald_recip_space_error(k, rc, alpha, vol);
      if (err < tol) {
        if (err > errmax) {
          errmax = err;
        }
        break;
      }
    }
    ikmax[d] = km;
  }

  return errmax;
}


INLINE ewald_data_t *ewald_data_new(const ewald_param_t *ewp,
    void *algo_data)
{
  ewald_data_t *ewd;

  XNEW(ewd, 1);

  if (ewp->rc > 0) { /* use the user-specified value */
    ewd->rc = ewp->rc;
  } else {
    ewd->rc = ewald_choose_rc(ewp);
  }

  if (ewp->sigma > 0) { /* use the user-specified value */
    ewd->alpha = 0.5/(ewp->sigma*ewp->sigma);
  } else {
    ewd->alpha = ewald_choose_alpha(ewp, ewd->rc, ewp->tol);
  }

  ewd->sqrt_alpha = sqrt(ewd->alpha);
  ewd->sqrt_alpha_over_pi = sqrt(ewd->alpha/PI);
  ewd->one_over_4alpha = 0.25/ewd->alpha;
  ewd->q2tot = ewald_total_charge2(ewp);
  ewd->qtot = ewald_total_charge(ewp);
  ewd->vol = ewald_box_volume(ewp);
  ewd->real_space_error = ewald_real_space_error(ewd->rc, ewd->alpha);
  ewald_recip_lattice(ewd->kbox, ewp);
  ewd->recip_space_error = ewald_ikmax(ewd->ikmax, ewp->tol, ewd->kbox, ewd->rc, ewd->alpha, ewd->vol);
  ewd->real_energy = 0;
  ewd->recip_energy = 0;
  ewd->self_energy = ewald_self_energy(ewp, ewd->alpha, ewd->q2tot);
  ewd->background_energy = ewald_background_energy(ewp, ewd->alpha, ewd->qtot, ewd->vol);
  printf("Ewald:\n"
      "alpha %g (1/sqrt(2*alpha) %g), rc %g\n"
      "real-space error %g, recip-space error %g, tol %g\n"
      "self %g, background %g, ikmax %d %d %d\n\n",
      ewd->alpha, 1/sqrt(2*ewd->alpha), ewd->rc,
      ewd->real_space_error, ewd->recip_space_error, ewp->tol,
      ewd->self_energy, ewd->background_energy,
      ewd->ikmax[0], ewd->ikmax[1], ewd->ikmax[2]);

  ewd->algo_data = algo_data;

  return ewd;
}


INLINE void ewald_data_delete(ewald_data_t *ewd)
{
  free(ewd);
}


#endif /* EWALD_BASIC_H__ */

