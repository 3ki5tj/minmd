#ifndef CONSTRAINT_SHAKE_H__
#define CONSTRAINT_SHAKE_H__


#include "utils.h"
#include "vec.h"
#include "constraint_basic.h"


typedef struct {
  int citmax; /* maximum number of iterations for SHAKE */
  int vitmax; /* maximum number of iterations for RATTLE */
  double ctol; /* tolerance for distance constraints */
  double vtol; /* tolerance for velocity constraints */
  constraint_pairlist_t *prlist;
} constraint_shake_param_t;


constraint_shake_param_t *constraint_shake_param_clone(constraint_shake_param_t *src)
{
  constraint_shake_param_t *dest;

  XCLONE(dest, src, sizeof(*src)); /* allocate space and shallow copy */
  dest->prlist = constraint_pairlist_clone(src->prlist);
  return dest;
}

void constraint_shake_param_delete(constraint_shake_param_t *shp)
{
  constraint_pairlist_delete(shp->prlist);
  free(shp);
}



typedef struct {
  int i, j; /* copy of the indices of the pair */
  real dist_ref2; /* reference distance squared */
  real imass, jmass; /* 1/mass[i], 1/mass[j] */
  real massij; /* 1/(1/mass[i] + 1/mass[j]) */
  real halfmassij; /* 0.5/(1/mass[i] + 1/mass[j]) */
  real dx[DIM]; /* displacement vectors of constraint pairs */
  real dist2; /* |dx|^2, for RATTLE */
} constraint_shake_pair_data_t;

typedef struct {
  int npr;
  constraint_shake_pair_data_t *prd;
} constraint_shake_data_t;


constraint_shake_data_t *constraint_shake_data_new(constraint_param_t *csp)
{
  constraint_shake_param_t *shp = (constraint_shake_param_t *) csp->iparam;
  constraint_pairlist_t *prlist = shp->prlist;
  constraint_pair_param_t *prp = prlist->pr;

  constraint_shake_data_t *shd;
  constraint_shake_pair_data_t *prd;

  XNEW(shd, 1);
  shd->npr = prlist->n;
  XNEW(prd, prlist->n);
  shd->prd = prd;

  int ipr, i, j;
  real dist;

  for (ipr = 0; ipr < prlist->n; ipr++, prp++, prd++) {
    prd->i = i = prp->i;
    prd->j = j = prp->j;
    dist = prp->dist;
    prd->dist_ref2 = dist*dist;
    prd->imass = 1.0/csp->mass[i];
    prd->jmass = 1.0/csp->mass[j];
    prd->massij = 1.0/(1.0/csp->mass[i] + 1.0/csp->mass[j]);
    prd->halfmassij = 0.5/(1.0/csp->mass[i] + 1.0/csp->mass[j]);
  }

  return shd;
}


void constraint_shake_data_delete(constraint_shake_data_t *shd)
{
  free(shd->prd);
  free(shd);
}



void constraint_shake_compute_dx(constraint_shake_data_t *shd,
                                 real (*x)[DIM],
                                 int compute_dist2)
{
  constraint_shake_pair_data_t *prd = shd->prd;
  int ipr, i, j;

  for (ipr = 0; ipr < shd->npr; ipr++, prd++) {
    i = prd->i;
    j = prd->j;

    /* compute the displacement vectors of constraint pairs */
    vec_diff(prd->dx, x[i], x[j]);

    /* for RATTLE */
    if (compute_dist2) {
      prd->dist2 = vec_sqr(prd->dx);
    }
  }
}


constraint_t *constraint_shake_new(constraint_param_t *param)
{
  constraint_t *cs;

  XNEW(cs, 1);
  cs->type = CONSTRAINT_TYPE_SHAKE;

  XCLONE(cs->param, param, sizeof(*param));
  cs->param->iparam = constraint_shake_param_clone((constraint_shake_param_t *) param->iparam);

  constraint_shake_data_t *shd = constraint_shake_data_new(cs->param);
  cs->data = constraint_data_new(cs->param, shd);

  return cs;
}


void constraint_shake_delete(constraint_t *cs)
{
  constraint_shake_param_delete((constraint_shake_param_t *) cs->param->iparam);
  free(cs->param);

  constraint_shake_data_delete(cs->data->idata);
  constraint_data_delete(cs->data);

  free(cs);
}


/* SHAKE algorithm */
void constraint_shake_coordinates(constraint_t *cs)
{
  constraint_param_t *csp = cs->param;
  constraint_shake_param_t *shp = (constraint_shake_param_t *) csp->iparam;
  constraint_shake_data_t *shd = (constraint_shake_data_t *) cs->data->idata;
  
  /* compute dx0 for the current configuration
   * this won't change during iterations */
  constraint_shake_compute_dx(shd, csp->x0, 0);

  int it;
  for (it = 0; it < shp->citmax; it++) {
    double maxdev = 0;
    constraint_shake_pair_data_t *prd = shd->prd;
    int ipr, i, j;
    real (*x1)[DIM] = csp->x1;
    real dx[DIM], dist2, dev, dot, dist_ref2, *dx0;

    for (ipr = 0; ipr < shd->npr; ipr++, prd++) {
      i = prd->i;
      j = prd->j;
      dist_ref2 = prd->dist_ref2;
      dx0 = prd->dx;

      /* TODO: replace to the PBC version */
      vec_diff(dx, x1[i], x1[j]);
      dist2 = vec_sqr(dx);
      dev = dist_ref2 - dist2;
      if (fabs(dev) > maxdev) {
        maxdev = fabs(dev);
      }
      dot = vec_dot(dx, dx0);
      if (dot < 0.2*dist_ref2) { /* a very bad case */
        dot = 0.2*dist_ref2;
      }
      dev *= prd->halfmassij/dot;
      vec_sinc(x1[i], dx0,  dev*prd->imass);
      vec_sinc(x1[j], dx0, -dev*prd->jmass);
    }
    
    /* one iteration finished */
    if (maxdev < shp->ctol) {
      break;
    }
  }
}


/* RATTLE algorithm */
void constraint_shake_velocities(constraint_t *cs)
{
  constraint_param_t *csp = cs->param;
  constraint_shake_param_t *shp = (constraint_shake_param_t *) csp->iparam;
  constraint_shake_data_t *shd = (constraint_shake_data_t *) cs->data->idata;

  /* compute dx and dist2 for the current configuration
   * they won't change during iterations */
  constraint_shake_compute_dx(shd, csp->x1, 1);

  int it;
  real (*v)[DIM] = csp->v;

  for (it = 0; it < shp->vitmax; it++) {
    double maxdev = 0;
    int ipr, i, j;
    constraint_shake_pair_data_t *prd = shd->prd;
    real dev;

    for (ipr = 0; ipr < shd->npr; ipr++, prd++) {
      real *dx = prd->dx,
           dv[DIM];

      i = prd->i;
      j = prd->j;
      vec_diff(dv, v[i], v[j]);
      dev = vec_dot(dv, dx);
      //printf("ROUND %d: ipr %d %d %d: dev %g, maxdev %g\n", it, ipr, i, j, dev, maxdev);
      if (fabs(dev) > maxdev) {
        maxdev = fabs(dev);
      }
      dev *= -prd->massij/prd->dist2;
      vec_sinc(v[i], dx,  dev*prd->imass);
      vec_sinc(v[j], dx, -dev*prd->jmass);
    }

    /* one iteration finished */
    if (maxdev < shp->vtol) {
      break;
    }
  }
}


int constraint_shake_apply(constraint_t *cs, unsigned flags)
{
  if (flags & CONSTRAINT_APPLY_COORDINATES) {
    constraint_shake_coordinates(cs);
  }

  if (flags & CONSTRAINT_APPLY_VELOCITIES) {
    constraint_shake_velocities(cs);
  }

  return 0;
}


#endif /* CONSTRAINT_SHAKE_H__ */

