#ifndef LJ_H__
#define LJ_H__



#include "minmd_lib.h"
#include "ljeos.h"



typedef struct {
  int n;
  real rho;
  real tp;
  real rc_def;
  real md_dt;
  int thermostat_type;
  real vr_dt;
} md_lj_param_t;


typedef struct {
  long long nsteps;
  int verbose;
  int nst_print;
  int do_stat;
} md_running_param_t;


typedef struct {
  real epot_shifted;
  real epot_shift;
  real epot_tail;
  real epot; /* unshifted with tail correction */
  real virial;
  real pres;
  real pres_tail;
} lj_epot_data_t;

lj_epot_data_t *lj_epot_data_init(real rc, int n, real rho)
{
  lj_epot_data_t *epdata;
  XNEW(epdata, 1);

  real irc = 1 / rc,
       irc3 = irc * irc * irc,
       irc6 = irc3 * irc3;

  epdata->epot_shift = 4 * irc6 * (irc6 - 1);
  epdata->epot_tail = 8*PI*rho*n/9*(irc6 - 3)*irc3;
  epdata->pres_tail = 32*PI*rho*rho/9*(irc6 - 1.5)*irc3;

  return epdata;
}

void lj_epot_data_free(lj_epot_data_t *epdata)
{
  free(epdata);
}


typedef struct {
  stat_accum_t *ekin_accum;
  stat_accum_t *epot_accum;
  stat_accum_t *vir_accum;
} md_lj_stat_t;


md_lj_stat_t *md_lj_stat_init(void)
{
  md_lj_stat_t *stat;
  XNEW(stat, 1);
  stat->ekin_accum = stat_accum_init();
  stat->epot_accum = stat_accum_init();
  stat->vir_accum = stat_accum_init();
  return stat;
}

void md_lj_stat_free(md_lj_stat_t *stat)
{
  stat_accum_free(stat->ekin_accum);
  stat_accum_free(stat->epot_accum);
  stat_accum_free(stat->vir_accum);
  free(stat);
}

void md_lj_stat_add(md_lj_stat_t *stat,
                 const lj_epot_data_t *epot_data,
                 real ekin)
{
  real epot = epot_data->epot,
       vir = epot_data->virial;
  stat_accum_add(stat->ekin_accum, ekin);
  stat_accum_add(stat->epot_accum, epot);
  stat_accum_add(stat->vir_accum, vir);
}



typedef struct {
  int n; /* number of particles */
  int n_dof; /* number of degrees of freedom */
  real rho;
  real vol;
  real l; /* box size */
  real rc;
  real tp;
  real md_dt;

  real *mass; /* masses */
  real (*x)[DIM]; /* positions */
  real (*v)[DIM]; /* velocities */
  real (*f)[DIM]; /* forces */

  lj_epot_data_t *epot_data; /* potential energy data */
  real ekin; /* kinetic energy */
  md_lj_stat_t *stat; /* statistical accumulators */

  rng_t *rng; /* random number generator */
  thermostat_t *thermostat;
} md_lj_t;


void md_lj_init_face_centered_lattice(int, real, real (*)[DIM], rng_t *);
real md_lj_force(md_lj_t *);

md_lj_t *md_lj_init(md_lj_param_t *param)
{
  md_lj_t *lj;
  int i, n = param->n;

  XNEW(lj, 1);
  lj->n = n;
  lj->n_dof = n * DIM - DIM; /* - DIM*(DIM+1)/2; // with angular part */
  lj->rho = param->rho;
  lj->vol = n / lj->rho;
  lj->l = pow(lj->vol, 1.0/DIM);
  lj->tp = param->tp;
  lj->md_dt = param->md_dt;

  /* initialize the cutoff distance */
  lj->rc = (param->rc_def < lj->l*0.5) ? param->rc_def : lj->l*0.5;

  /* initialize the potential energy */
  lj->epot_data = lj_epot_data_init(lj->rc, n, lj->rho);

  lj->rng = rng_init(0, 0);

  XNEW(lj->mass, n);
  for (i = 0; i < n; i++) {
    lj->mass[i] = 1;
  }

  XNEW(lj->x, n);
  XNEW(lj->v, n);
  XNEW(lj->f, n);

  /* initialize positions */
  md_lj_init_face_centered_lattice(n, lj->l, lj->x, lj->rng);
  md_lj_force(lj);

  /* initialize velocities at the reference temperature */
  mdutils_init_velocities(n, lj->mass, lj->v, param->tp, lj->rng);
  lj->ekin = mdutils_ekin(n, lj->mass, lj->v);

  /* initialize the thermostat */
  thermostat_vrescaling_param_t vrp = {
    .rng = lj->rng
  };
  thermostat_param_t ts_param = {
    .n = n,
    .n_dof = n*DIM - DIM,
    .tp = param->tp,
    .dt = param->vr_dt,
    .boltz = 1,
    .mass = lj->mass,
    .v = lj->v,
    .algo_param = &vrp,
  };
  lj->thermostat = thermostat_init(param->thermostat_type, &ts_param);

  /* initialize the statistical accumulators */
  lj->stat = md_lj_stat_init();

  return lj;
}


void md_lj_free(md_lj_t *lj)
{
  thermostat_free(lj->thermostat);
  md_lj_stat_free(lj->stat);
  free(lj->mass);
  free(lj->x);
  free(lj->v);
  free(lj->f);
  lj_epot_data_free(lj->epot_data);
  rng_free(lj->rng);
  free(lj);
}



void md_lj_init_face_centered_lattice(int n, real side_length, real (*x)[DIM], rng_t *r)
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


INLINE real *md_lj_vec_pbc(real *x, real l, real invl)
{
  int d;
  for (d = 0; d < DIM; d++) {
    x[d] -= ((int)(x[d]*invl + 100.5) - 100.)*l;
  }
}


INLINE real md_lj_pbc_dist2(real *dx, real *xi, real *xj, real l, real invl)
{
 vec_diff(dx, xi, xj);
 md_lj_vec_pbc(dx, l, invl);
 return vec_sqr(dx);
}


real md_lj_force(md_lj_t *lj)
{
  real dx[DIM], fi[DIM], dr2, ir2, ir6, fs;
  real (*x)[DIM] = lj->x;
  real (*f)[DIM] = lj->f;
  real rc2 = lj->rc * lj->rc;
  real l = lj->l, invl = 1/l;
  int i, j, npr = 0, n = lj->n;

  for (i = 0; i < n; i++) {
    vec_zero(f[i]);
  }

  double epot6 = 0, epot12 = 0;
  for (i = 0; i < n - 1; i++) {
    vec_zero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = md_lj_pbc_dist2(dx, x[i], x[j], l, invl);
      if (dr2 > rc2) {
        continue;
      }
      ir2 = 1/dr2;
      ir6 = ir2 * ir2 * ir2;
      fs = ir6 * (48*ir6 - 24); /* f*r */
      fs /= dr2; /* f*r / r^2 */
      vec_sinc(fi,   dx,  fs);
      vec_sinc(f[j], dx, -fs);
      epot6 += ir6;
      epot12 += ir6 * ir6;
      npr++;
    }
    vec_inc(f[i], fi);
  }

  epot6 *= 4;
  epot12 *= 4;
  real epot_pr = (real) (epot12 - epot6);
  lj_epot_data_t *epdata = lj->epot_data;
  epdata->epot_shifted = epot_pr - npr * epdata->epot_shift;
  epdata->virial = (real) (12 * epot12 - 6 * epot6);
  epdata->epot = epot_pr + epdata->epot_tail; /* unshifted with tail correction */
  return epdata->epot;
}


void md_lj_vv(md_lj_t *lj)
{
  int i, n = lj->n;
  double dt = lj->md_dt,
         dth = dt*0.5;

  /* velocity verlet part 1 */
  for (i = 0; i < n; i++) {
    vec_sinc(lj->v[i], lj->f[i], dth/lj->mass[i]);
    vec_sinc(lj->x[i], lj->v[i], dt);
  }

  md_lj_force(lj);

  /* velocity verlet part 2 */
  for (i = 0; i < n; i++) {
    vec_sinc(lj->v[i], lj->f[i], dth/lj->mass[i]);
  }
}


void md_lj_step(md_lj_t *lj, const md_running_param_t *mrp)
{
  lj->ekin = thermostat_apply(lj->thermostat);
  md_lj_vv(lj);
  lj->ekin = thermostat_apply(lj->thermostat);
}


void md_lj_run(md_lj_t *lj, const md_running_param_t *mrp)
{
  long long step;

  for (step = 1; step <= mrp->nsteps; step++) {
    md_lj_step(lj, mrp);

    if (mrp->do_stat) {
      md_lj_stat_add(lj->stat, lj->epot_data, lj->ekin);
    }

    if (mrp->verbose > 1 && step % mrp->nst_print == 0) {
      fprintf(stderr, "step %8lld: tp %6.4f, epot/n %6.3f",
          step, lj->ekin*2/lj->n_dof, lj->epot_data->epot/lj->n);
      if (lj->thermostat->type == THERMOSTAT_TYPE_NULL) {
        /* total energy with the potentital part shifted, should be conserved */
        double etot_shifted = lj->epot_data->epot_shifted + lj->ekin;
        fprintf(stderr, ", epot_shifted %g, etot_shifted %g",
            lj->epot_data->epot_shifted, etot_shifted);
      }
      fprintf(stderr, "\n");
    }
  }
}



void md_lj_print_stat(md_lj_t *lj)
{
  md_lj_stat_t *stat = lj->stat;
  int n = lj->n, n_dof = lj->n_dof;
  double ekin_mean, ekin_std;
  double epot_mean, epot_std;
  double vir_mean, vir_std;
  double pres_mean, pres_std;
  ekin_std = stat_accum_get_std(stat->ekin_accum, &ekin_mean);
  epot_std = stat_accum_get_std(stat->epot_accum, &epot_mean);
  vir_std = stat_accum_get_std(stat->vir_accum, &vir_mean);
  pres_mean = lj->rho * lj->tp
            + vir_mean / (DIM * lj->vol)
            + lj->epot_data->pres_tail;
  pres_std = vir_std / (DIM * lj->vol);

  // get the equation of state data
  real epot_ref, pres_ref, fex_ref, muex_ref;
  epot_ref = ljeos3d_get(lj->rho, lj->tp, &pres_ref, &fex_ref, &muex_ref);

  printf("  nsteps %8lld\n"
         "       n %8d\n"
         "      tp %8.4f+/-%6.4f (cf: %8.4f)\n"
         "  epot/n %8.4f+/-%6.4f (cf: %8.4f)\n"
         "    pres %8.4f+/-%6.4f (cf: %8.4f)\n"
         "(NOTE: numbers after \"+/-\" are standard deviations, not estimated errors)\n",
      stat->ekin_accum->count, n,
      2*ekin_mean/n_dof, 2*ekin_std/n_dof, lj->tp,
      epot_mean/n, epot_std/n, epot_ref,
      pres_mean, pres_std, pres_ref);
}



#endif /* LJ_H__ */


