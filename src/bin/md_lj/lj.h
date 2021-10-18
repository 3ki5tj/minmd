#ifndef LJ_H__
#define LJ_H__


#include "minmd_lib.h"



typedef struct {
  real rc_def;
  real rc;
} lj_cutoff_t;

void lj_cutoff_init(lj_cutoff_t *c, real rc_def, real l)
{
  c->rc_def = rc_def;
  c->rc = (rc_def < l*0.5) ? rc_def : l*0.5;
}



typedef struct {
  real epot, epot6, epot12;
  real epot_shift;
  real epot_tail;
  real vir;
  real pres;
  real pres_tail;
} lj_epot_data_t;


void lj_epot_data_init(lj_epot_data_t *epdata, lj_cutoff_t *cutoff, int n, real rho)
{
    real rc = cutoff->rc,
         irc = 1 / rc,
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
}


typedef struct {
  stat_accum_t *epot_acc;
  stat_accum_t *ekin_acc;
} lj_stat_t;


void lj_stat_init(lj_stat_t *stat_acc)
{
  stat_accum_init(stat_acc->epot_acc);
  stat_accum_init(stat_acc->ekin_acc);
}

void lj_stat_add(lj_stat_t *stat_acc,
                 const lj_epot_data_t *epot_data,
                 real ekin)
{
  stat_accum_add(stat_acc->epot_acc, epot_data->epot);
  stat_accum_add(stat_acc->ekin_acc, ekin);
}



typedef struct {
  int dim; /* dimension, 3 or 2 */
  int n; /* number of particles */
  int n_dof; /* number of degrees of freedom */
  real rho;
  real vol;
  real l; /* box size */
  lj_cutoff_t cutoff;

  real *mass; /* masses */
  real (*x)[DIM]; /* positions */
  real (*v)[DIM]; /* velocities */
  real (*f)[DIM]; /* forces */
  real *r2ij; /* cached pair distances */

  lj_epot_data_t epot_data; /* potential energy data */
  real ekin; /* kinetic energy */
  lj_stat_t stat_acc; /* statistical accumulators */

  rng_t *rng; /* random number generator */
} lj_md_t;



lj_md_t *lj_md_open(int n, real rho, real rc_def)
{
  lj_md_t *lj;
  int i;

  lj->dim = DIM;
  lj->n = n;
  lj->n_dof = n * DIM - DIM; /* - DIM*(DIM+1)/2; // with angular part */
  lj->rho = rho;
  lj->vol = n / rho;
  lj->l = pow(lj->vol, 1.0/DIM);

  /* initialize the cutoff parameters */
  lj_cutoff_init(&(lj->cutoff), rc_def, lj->l);

  /* initialize the potential energy */
  lj_epot_data_init(&(lj->epot_data), &(lj->cutoff), n, rho);

  lj->rng = rng_init(0, 0);

  XNEW(lj->mass, n);
  for (i = 0; i < n; i++) {
    lj->mass[i] = 1;
  }

  XNEW(lj->x, n);
  XNEW(lj->v, n);
  XNEW(lj->f, n);
  XNEW(lj->r2ij, n*n);

  /* initialize positions */
  mdutils_init_face_centered_lattice(n, lj->l, lj->x, lj->rng);

  /* initialize velocities at the reference temperature */
  real tp_ref = 1;
  mdutils_init_velocities(n, lj->mass, lj->v, tp_ref, lj->rng); 

  stat_accum_init(&(lj->stat_acc));

  return lj;
}



void lj_md_run(lj_md_t *lj, long long nsteps)
{
}



#endif /* LJ_H__ */


