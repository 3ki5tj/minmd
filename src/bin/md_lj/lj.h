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
  real epot_pr, epot6, epot12;
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
       irc2 = irc * irc,
       irc3 = irc2 * irc,
       irc4 = irc2 * irc2,
       irc6 = irc4 * irc2;

  epdata->epot_shift = 4 * irc6 * (irc6 - 1);
#ifdef MINMD_2D
  epdata->epot_tail = PI*rho*n*(.4*irc6 - 1)*irc4;
  epdata->pres_tail = PI*rho*rho*(2.4*irc6 - 3)*irc4;
#else
  epdata->epot_tail = 8*PI*rho*n/9*(irc6 - 3)*irc3;
  epdata->pres_tail = 32*PI*rho*rho/9*(irc6 - 1.5)*irc3;
#endif

  return epdata;
}


void lj_epot_data_free(lj_epot_data_t *epdata)
{
  free(epdata);
}


typedef struct {
  stat_accum_t *ekin_accum;
  stat_accum_t *epot_accum;
  stat_accum_t *etot_accum;
} md_lj_stat_t;


md_lj_stat_t *md_lj_stat_init(void)
{
  md_lj_stat_t *stat;
  XNEW(stat, 1);
  stat->ekin_accum = stat_accum_init();
  stat->epot_accum = stat_accum_init();
  stat->etot_accum = stat_accum_init();
  return stat;
}

void md_lj_stat_free(md_lj_stat_t *stat)
{
  stat_accum_free(stat->ekin_accum);
  stat_accum_free(stat->epot_accum);
  stat_accum_free(stat->etot_accum);
  free(stat);
}

void md_lj_stat_add(md_lj_stat_t *stat,
                 const lj_epot_data_t *epot_data,
                 real ekin)
{
  real epot = epot_data->epot,
       etot = ekin + epot;
  stat_accum_add(stat->ekin_accum, ekin);
  stat_accum_add(stat->epot_accum, epot);
  stat_accum_add(stat->etot_accum, etot);
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
  real etot; /* total energy */
  md_lj_stat_t *stat; /* statistical accumulators */

  rng_t *rng; /* random number generator */
  thermostat_t *thermostat;
} md_lj_t;


real md_lj_force(md_lj_t *);

md_lj_t *md_lj_open(md_lj_param_t *param)
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
  mdutils_init_face_centered_lattice(n, lj->l, lj->x, lj->rng);
  md_lj_force(lj);

  /* initialize velocities at the reference temperature */
  mdutils_init_velocities(n, lj->mass, lj->v, param->tp, lj->rng); 
  lj->ekin = mdutils_ekin(n, lj->mass, lj->v);

  /* compute the total energy */
  lj->etot = lj->ekin + lj->epot_data->epot;

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

  real epot6 = 0, epot12 = 0;
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
  lj_epot_data_t *epdata = lj->epot_data;
  epdata->epot_pr = epot12 - epot6;
  epdata->epot6 = epot6;
  epdata->epot12 = epot12;
  epdata->epot_shifted = epdata->epot_pr - npr * epdata->epot_shift;
  epdata->virial = 12 * epot12 - 6 * epot6;
  epdata->epot = epdata->epot_pr + epdata->epot_tail; /* unshifted with tail correction */
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
  lj->etot = lj->ekin + lj->epot_data->epot;
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
      fprintf(stderr, "step %lld: ekin %g + epot %g = etot %g\n",
          step, lj->ekin, lj->epot_data->epot, lj->etot);
    }
  }
}



void md_lj_print_stat(md_lj_t *lj)
{
  md_lj_stat_t *stat = lj->stat;
  int n = lj->n, n_dof = lj->n_dof;
  double ekin_mean, ekin_std;
  double epot_mean, epot_std;
  double etot_mean, etot_std;
  ekin_std = stat_accum_get_std(stat->ekin_accum, &ekin_mean);
  epot_std = stat_accum_get_std(stat->epot_accum, &epot_mean);
  etot_std = stat_accum_get_std(stat->etot_accum, &etot_mean);

  real epot_ref, pres_ref, fex_ref, muex_ref;
  epot_ref = ljeos3d_get(lj->rho, lj->tp, &pres_ref, &fex_ref, &muex_ref);

  printf("  nsteps %lld\n"
         "       n %d\n"
         "      tp %g+/-%g (cf: %g)\n"
         "  epot/n %g+/-%g (cf: %g)\n",
      stat->ekin_accum->count, n,
      2*ekin_mean/n_dof, 2*ekin_std/n_dof, lj->tp,
      epot_mean/n, epot_std/n, epot_ref);
}



#endif /* LJ_H__ */


