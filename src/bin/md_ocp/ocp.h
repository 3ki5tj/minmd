#ifndef OCP_H__
#define OCP_H__



#include "minmd_lib.h"



typedef struct {
  int n;
  real rho;
  real tp;
  real rc_def;
  real md_dt;
  int thermostat_type;
  real vr_dt;
  real ewald_tol;
} md_ocp_param_t;


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
} epot_vdw_data_t;

epot_vdw_data_t *epot_vdw_data_new(real rc, int n, real rho)
{
  epot_vdw_data_t *epdata;
  XNEW(epdata, 1);

  real irc = 1 / rc,
       irc3 = irc * irc * irc,
       irc6 = irc3 * irc3;

  epdata->epot_shift = 4 * irc6 * (irc6 - 1);
  epdata->epot_tail = 8*PI*rho*n/9*(irc6 - 3)*irc3;
  epdata->pres_tail = 32*PI*rho*rho/9*(irc6 - 1.5)*irc3;

  return epdata;
}

void epot_vdw_data_delete(epot_vdw_data_t *epdata)
{
  free(epdata);
}


typedef struct {
  real epot;
  real epot_shifted;
  real virial;
} epot_charge_data_t;

epot_charge_data_t *epot_charge_data_new(real rc, int n, real rho)
{
  epot_charge_data_t *epdata;
  XNEW(epdata, 1);

  epdata->epot = 0;
  epdata->epot_shifted = 0;
  epdata->virial = 0;

  return epdata;
}

void epot_charge_data_delete(epot_charge_data_t *epdata)
{
  free(epdata);
}



typedef struct {
  real epot;
  real epot_shifted;
  real virial;
  epot_vdw_data_t *vdw;
  epot_charge_data_t *charge;
} epot_data_t;

epot_data_t *epot_data_new(real rc, int n, real rho)
{
  epot_data_t *epdata;
  XNEW(epdata, 1);

  epdata->epot = 0;
  epdata->epot_shifted = 0;
  epdata->virial = 0;
  epdata->vdw = epot_vdw_data_new(rc, n, rho);
  epdata->charge = epot_charge_data_new(rc, n, rho);

  return epdata;
}

void epot_data_aggregate(epot_data_t *epdata)
{
  epdata->epot = epdata->vdw->epot + epdata->charge->epot;
  epdata->epot_shifted = epdata->vdw->epot_shifted + epdata->charge->epot_shifted;
  epdata->virial = epdata->vdw->virial + epdata->charge->virial;
}

void epot_data_delete(epot_data_t *epdata)
{
  epot_vdw_data_delete(epdata->vdw);
  epot_charge_data_delete(epdata->charge);
  free(epdata);
}



typedef struct {
  stat_accum_t *ekin_accum;
  stat_accum_t *epot_accum;
} md_ocp_stat_t;


md_ocp_stat_t *md_ocp_stat_new(void)
{
  md_ocp_stat_t *stat;
  XNEW(stat, 1);
  stat->ekin_accum = stat_accum_new();
  stat->epot_accum = stat_accum_new();
  return stat;
}

void md_ocp_stat_delete(md_ocp_stat_t *stat)
{
  stat_accum_delete(stat->ekin_accum);
  stat_accum_delete(stat->epot_accum);
  free(stat);
}

void md_ocp_stat_add(md_ocp_stat_t *stat,
                 const epot_data_t *epot_data,
                 real ekin)
{
  real epot = epot_data->epot;
  stat_accum_add(stat->ekin_accum, ekin);
  stat_accum_add(stat->epot_accum, epot);
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
  real *charge; /* charges */
  real (*x)[DIM]; /* positions */
  real (*v)[DIM]; /* velocities */
  real (*f)[DIM]; /* forces */

  ewald_t *ew;

  epot_data_t *epot_data; /* potential energy data */
  real ekin; /* kinetic energy */
  md_ocp_stat_t *stat; /* statistical accumulators */

  rng_t *rng; /* random number generator */
  thermostat_t *thermostat;
} md_ocp_t;


void md_ocp_init_face_centered_lattice(int, real, real (*)[DIM], rng_t *);
real md_ocp_force(md_ocp_t *);

md_ocp_t *md_ocp_new(md_ocp_param_t *param)
{
  md_ocp_t *ocp;
  int i, n = param->n;

  XNEW(ocp, 1);
  ocp->n = n;
  ocp->n_dof = n * DIM - DIM; /* - DIM*(DIM+1)/2; // with angular part */
  ocp->rho = param->rho;
  ocp->vol = n / ocp->rho;
  ocp->l = pow(ocp->vol, 1.0/DIM);
  ocp->tp = param->tp;
  ocp->md_dt = param->md_dt;

  /* initialize the cutoff distance */
  ocp->rc = (param->rc_def < ocp->l*0.5) ? param->rc_def : ocp->l*0.5;

  /* initialize the potential energy */
  ocp->epot_data = epot_data_new(ocp->rc, n, ocp->rho);

  ocp->rng = rng_new(0, 0);

  XNEW(ocp->mass, n);
  for (i = 0; i < n; i++) {
    ocp->mass[i] = 1;
  }

  XNEW(ocp->charge, n);
  for (i = 0; i < n; i++) {
    ocp->charge[i] = 1;
  }

  XNEW(ocp->x, n);
  XNEW(ocp->v, n);
  XNEW(ocp->f, n);

  /* initialize positions */
  md_ocp_init_face_centered_lattice(n, ocp->l, ocp->x, ocp->rng);

  /* initialize the ewald method for charge interaction */
  double ewald_tol = param->ewald_tol > 0 ? param->ewald_tol : 1e-6;
  ewald_param_t ewp = {
    .n = n,
    .box = {ocp->l, ocp->l, ocp->l},
    .rc = ocp->rc,
    .sigma = 0.0, /* choose alpha automatically */
    .tol = ewald_tol,
    .kee = 1.0,
    .charge = ocp->charge,
    .x = ocp->x,
    .f = ocp->f,
  };
  ocp->ew = ewald_new(EWALD_TYPE_DIRECT, &ewp);

  md_ocp_force(ocp);

  /* initialize velocities at the reference temperature */
  mdutils_init_velocities(n, ocp->mass, ocp->v, param->tp, ocp->rng); 
  ocp->ekin = mdutils_ekin(n, ocp->mass, ocp->v);

  /* initialize the thermostat */
  thermostat_vrescaling_param_t vrp = {
    .rng = ocp->rng
  };
  thermostat_param_t ts_param = {
    .n = n,
    .n_dof = n*DIM - DIM,
    .tp = param->tp,
    .dt = param->vr_dt,
    .boltz = 1,
    .mass = ocp->mass,
    .v = ocp->v,
    .algo_param = &vrp,
  };
  ocp->thermostat = thermostat_new(param->thermostat_type, &ts_param);

  /* initialize the statistical accumulators */
  ocp->stat = md_ocp_stat_new();

  return ocp;
}


void md_ocp_delete(md_ocp_t *ocp)
{
  ewald_delete(ocp->ew);
  thermostat_delete(ocp->thermostat);
  md_ocp_stat_delete(ocp->stat);
  free(ocp->mass);
  free(ocp->charge);
  free(ocp->x);
  free(ocp->v);
  free(ocp->f);
  epot_data_delete(ocp->epot_data);
  rng_delete(ocp->rng);
  free(ocp);
}



void md_ocp_init_face_centered_lattice(int n, real side_length, real (*x)[DIM], rng_t *r)
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


INLINE real *md_ocp_vec_pbc(real *x, real l, real invl)
{
  int d;
  for (d = 0; d < DIM; d++) {
    x[d] -= ((int)(x[d]*invl + 100.5) - 100.)*l;
  }
}


INLINE real md_ocp_pbc_dist2(real *dx, real *xi, real *xj, real l, real invl)
{
 vec_diff(dx, xi, xj);
 md_ocp_vec_pbc(dx, l, invl);
 return vec_sqr(dx);
}


real md_ocp_force_vdw(md_ocp_t *ocp, epot_vdw_data_t *epdata)
{
  real dx[DIM], fi[DIM], dr2, ir2, ir6, fs;
  real (*x)[DIM] = ocp->x;
  real (*f)[DIM] = ocp->f;
  real rc2 = ocp->rc * ocp->rc;
  real l = ocp->l, invl = 1/l;
  int i, j, npr = 0, n = ocp->n;

  for (i = 0; i < n; i++) {
    vec_zero(f[i]);
  }

  double epot6 = 0, epot12 = 0;
  for (i = 0; i < n - 1; i++) {
    vec_zero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = md_ocp_pbc_dist2(dx, x[i], x[j], l, invl);
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
  epdata->epot_shifted = epot_pr - npr * epdata->epot_shift;
  epdata->virial = (real) (12 * epot12 - 6 * epot6);
  epdata->epot = epot_pr + epdata->epot_tail; /* unshifted with tail correction */
  return epdata->epot;
}


real md_ocp_force_charge(md_ocp_t *ocp, epot_charge_data_t *epdata)
{
  /* call the ewald method to compute the electrostatic force */
  ewald_force_options_t ewf_opt = {
    /* electrostatic force is computed after the van der Waals forces
     * so there is no need to clear the force */
    .zero_forces = 0,
  };
  ewald_force(ocp->ew, &ewf_opt);

  epdata->epot = ocp->ew->data->energy;
  epdata->virial = ocp->ew->data->virial;
  return epdata->epot;
}


real md_ocp_force(md_ocp_t *ocp)
{
  md_ocp_force_vdw(ocp, ocp->epot_data->vdw);
  md_ocp_force_charge(ocp, ocp->epot_data->charge);

  epot_data_aggregate(ocp->epot_data);
  return ocp->epot_data->epot;
}


void md_ocp_vv(md_ocp_t *ocp)
{
  int i, n = ocp->n;
  double dt = ocp->md_dt,
         dth = dt*0.5;

  /* velocity verlet part 1 */
  for (i = 0; i < n; i++) {
    vec_sinc(ocp->v[i], ocp->f[i], dth/ocp->mass[i]);
    vec_sinc(ocp->x[i], ocp->v[i], dt);
  }

  md_ocp_force(ocp);
  
  /* velocity verlet part 2 */
  for (i = 0; i < n; i++) {
    vec_sinc(ocp->v[i], ocp->f[i], dth/ocp->mass[i]);
  }
}


void md_ocp_step(md_ocp_t *ocp, const md_running_param_t *mrp)
{
  ocp->ekin = thermostat_apply(ocp->thermostat);
  md_ocp_vv(ocp);
  ocp->ekin = thermostat_apply(ocp->thermostat);
}


void md_ocp_def_logger(md_ocp_t *ocp, const md_running_param_t *mrp, long long step)
{
  fprintf(stderr, "step %lld: tp %g, epot/n %g",
      step, ocp->ekin*2/ocp->n_dof, ocp->epot_data->epot/ocp->n);
  if (ocp->thermostat->type == THERMOSTAT_TYPE_NULL) {
    double epot_shifted = ocp->epot_data->vdw->epot_shifted
                        + ocp->epot_data->charge->epot;
    double etot_shifted = epot_shifted + ocp->ekin;
    fprintf(stderr, ", epot_shifted %g, etot_shifted %g",
        epot_shifted, etot_shifted);
  }
  fprintf(stderr, "\n");
}


void md_ocp_run(md_ocp_t *ocp, const md_running_param_t *mrp)
{
  long long step;

  for (step = 1; step <= mrp->nsteps; step++) {
    md_ocp_step(ocp, mrp);

    if (mrp->do_stat) {
      md_ocp_stat_add(ocp->stat, ocp->epot_data, ocp->ekin);
    }

    if (mrp->verbose > 1 && step % mrp->nst_print == 0) {
      md_ocp_def_logger(ocp, mrp, step);
    }
  }
}



void md_ocp_print_stat(md_ocp_t *ocp)
{
  md_ocp_stat_t *stat = ocp->stat;
  int n = ocp->n, n_dof = ocp->n_dof;
  double ekin_mean, ekin_std;
  double epot_mean, epot_std;
  ekin_std = stat_accum_get_std(stat->ekin_accum, &ekin_mean);
  epot_std = stat_accum_get_std(stat->epot_accum, &epot_mean);

  printf("  nsteps %lld\n"
         "       n %d\n"
         "      tp %g+/-%g (cf: %g)\n"
         "  epot/n %g+/-%g\n",
      stat->ekin_accum->count, n,
      2*ekin_mean/n_dof, 2*ekin_std/n_dof, ocp->tp,
      epot_mean/n, epot_std/n);
}



#endif /* OCP_H__ */


