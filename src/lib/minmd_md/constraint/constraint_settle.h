#ifndef CONSTRAINT_SETTLE_H__
#define CONSTRAINT_SETTLE_H__

#include "utils.h"
#include "vec.h"
#include "constraint_basic.h"


/* the following parameters are only for the symmetric water molecule */
typedef struct {
  real dist_oh; /* in angstroms */
  real ang_hoh; /* in degrees */
  real mass_o, mass_h;
} constraint_settle_model_param_t;

enum {
  CONSTRAINT_SETTLE_WATER_SPC = 1,
  CONSTRAINT_SETTLE_WATER_TIP3P,
};

constraint_settle_model_param_t constraint_settle_water_spc_param = {
  .dist_oh = 1.0,
  .ang_hoh = 109.28,
  .mass_o = 12.01100,
  .mass_h =  1.00800,
};

constraint_settle_model_param_t constraint_settle_water_tip3p_param = {
  .dist_oh = 0.9572,
  .ang_hoh = 104.52,
  .mass_o = 12.01100,
  .mass_h =  1.00800,
};


/* for SETTLE, no need to fill the `n`, `mass` fields of constraint_param_t */
typedef struct {
  constraint_settle_model_param_t model_param;
  int nmol;
  int *id_o;
} constraint_settle_param_t;


enum {
  CONSTRAINT_SETTLE_NATOM = 3,
  CONSTRAINT_SETTLE_NPAIR = 3,
};

int constraint_settle_pair_ids_[CONSTRAINT_SETTLE_NPAIR][2] = {
  {1, 2},
  {2, 0},
  {0, 1},
};


typedef struct {
  real dist_oh, dist_hh, dist_on, dist_hn;
  real ang_hoh; /* in radians */
  real mass_o, mass_h;
  real mass_total;
  real frac_mass_h;
  int nmol;
  int *id_o;

  /* for the velocity algorithm */
  real invmass[CONSTRAINT_SETTLE_NATOM];
  real vinvmat[CONSTRAINT_SETTLE_NPAIR][CONSTRAINT_SETTLE_NPAIR];
} constraint_settle_data_t;


/* b = a^(-1) */
INLINE void constraint_settle_invert_matrix_(real b[DIM][DIM], real a[DIM][DIM])
{
   real as[DIM][DIM], det;
   int i, j;

   as[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
   as[0][1] = a[2][1]*a[0][2] - a[2][2]*a[0][1];
   as[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
   
   as[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
   as[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
   as[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];

   as[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
   as[2][1] = a[2][0]*a[0][1] - a[2][1]*a[0][0];
   as[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

   det = as[0][0] * a[0][0] + as[0][1] * a[1][0] + as[0][2] * a[2][0]; 
   for (i = 0; i < DIM; i++) {
     for (j = 0; j < DIM; j++) {
       b[i][j] = as[i][j]/det;
     }
   }
}

INLINE void constraint_settle_print_vec_(double *v, const char *name)
{
  printf("%s: %+8.5f %+8.5f %8.5f\n", name, v[0], v[1], v[2]);
}

INLINE void constraint_settle_print_mat_(double (*m)[3], const char *name)
{
  int i;
  printf("%s:\n", name);
  for (i = 0; i < 3; i++) {
    printf("    %+8.5f %+8.5f %8.5f\n", m[i][0], m[i][1], m[i][2]);
  }
  printf("\n");
}



void constraint_settle_data_new_velocities(constraint_settle_data_t *sed)
{
  real x[CONSTRAINT_SETTLE_NATOM][DIM] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  real dxpr[CONSTRAINT_SETTLE_NPAIR][DIM];
  real vmat[CONSTRAINT_SETTLE_NPAIR][CONSTRAINT_SETTLE_NPAIR], invmass;
  int ipr, jpr, i0, i1, j0;

  sed->invmass[0] = 1.0 / sed->mass_o;
  sed->invmass[1] = 1.0 / sed->mass_h;
  sed->invmass[2] = 1.0 / sed->mass_h;

  x[1][0] = -sed->dist_on;
  x[1][1] = +sed->dist_hn;
  x[2][0] = -sed->dist_on;
  x[2][1] = -sed->dist_hn;

  for (ipr = 0; ipr < CONSTRAINT_SETTLE_NPAIR; ipr++) {
    i0 = constraint_settle_pair_ids_[ipr][0];
    i1 = constraint_settle_pair_ids_[ipr][1];
    vec_diff(dxpr[ipr], x[i0], x[i1]);
  }

  for (ipr = 0; ipr < CONSTRAINT_SETTLE_NPAIR; ipr++) {
    i0 = constraint_settle_pair_ids_[ipr][0];
    i1 = constraint_settle_pair_ids_[ipr][1];
    for (jpr = 0; jpr < CONSTRAINT_SETTLE_NPAIR; jpr++) {
      if (ipr == jpr) {
        invmass = -(sed->invmass[i0] + sed->invmass[i1]);
      } else {
        j0 = constraint_settle_pair_ids_[jpr][0];
        // int j1 = constraint_settle_pair_ids_[jpr][1];
        if (i1 == j0) {
          invmass = sed->invmass[i1];
        } else { // i0 == j1
          invmass = sed->invmass[i0];
        }
      }
      vmat[ipr][jpr] = invmass * vec_dot(dxpr[ipr], dxpr[jpr]);
    }
  }

  constraint_settle_invert_matrix_(sed->vinvmat, vmat);
  //constraint_settle_print_mat_(sed->vinvmat, "Inverse velocity matrix");
}

constraint_settle_data_t *constraint_settle_data_new(constraint_param_t *csp)
{
  constraint_settle_param_t *sep = (constraint_settle_param_t *) csp->iparam;
  constraint_settle_model_param_t *model_param = &(sep->model_param);
  constraint_settle_data_t *sed;
  int i;

  XNEW(sed, 1);
  
  sed->dist_oh = model_param->dist_oh;
  sed->ang_hoh = model_param->ang_hoh * PI / 180;
  sed->dist_on = sed->dist_oh * cos(sed->ang_hoh/2);
  sed->dist_hn = sed->dist_oh * sin(sed->ang_hoh/2);
  sed->dist_hh = 2 * sed->dist_hn;
  //printf("OH %g, HH: %g; ON %g\n", sed->dist_oh, sed->dist_hh, sed->dist_on);

  sed->mass_o = model_param->mass_o;
  sed->mass_h = model_param->mass_h;
  sed->mass_total = sed->mass_o + sed->mass_h*2;
  sed->frac_mass_h = sed->mass_h / sed->mass_total;

  sed->nmol = sep->nmol;
  XNEW(sed->id_o, sed->nmol);
  for (i = 0; i < sed->nmol; i++) {
    sed->id_o[i] = sep->id_o[i];
  }

  constraint_settle_data_new_velocities(sed);

  return sed;
}



void constraint_settle_data_delete(constraint_settle_data_t *sed)
{
  free(sed->id_o);
  free(sed);
}



constraint_t *constraint_settle_new(constraint_param_t *param)
{
  constraint_t *cs;

  XNEW(cs, 1);
  cs->type = CONSTRAINT_TYPE_SETTLE;

  XCLONE(cs->param, param, sizeof(*param));
  XCLONE(cs->param->iparam, param->iparam, sizeof(constraint_settle_param_t));

  constraint_settle_data_t *sed = constraint_settle_data_new(cs->param);
  cs->data = constraint_data_new(cs->param, sed);

  return cs;
}


void constraint_settle_delete(constraint_t *cs)
{
  free(cs->param->iparam);
  free(cs->param);

  constraint_settle_data_delete(cs->data->idata);
  constraint_data_delete(cs->data);

  free(cs);
}


void constraint_settle_one_coordinates(constraint_settle_data_t *sed,
    real (*x0)[DIM], real (*x1)[DIM], real (*xout)[DIM],
    real frac_mass_h)
{
  enum {
    ID_A = 0,
    ID_B = 1,
    ID_C = 2,
  };

  real xab0[DIM], xac0[DIM], xz0[DIM], z0[DIM], dz0;
  real xab1[DIM], xac1[DIM], zab, zac;
  real vaxis[DIM], daxis2, daxis;
  real xab4[DIM], xac4[DIM];
  real sinth, sinth2, one_minus_costh;
  real xab0para[DIM], xab0perp[DIM], dxab0[DIM];
  real xac0para[DIM], xac0perp[DIM], dxac0[DIM];

  // parameter for the second rotation
  real xcom0[DIM], xbo0[DIM], xco0[DIM];
  real dxab[DIM], dxac[DIM];
  real tmpab[DIM], tmpac[DIM], tmpsum[DIM];
  real rhs, a, b, gam, a2b2;
  real xab4para[DIM], xab4perp[DIM], dxab4[DIM], yab4[DIM], xab5[DIM];
  real xac4para[DIM], xac4perp[DIM], dxac4[DIM], yac4[DIM], xac5[DIM];
  real x1sum[DIM], x5sum[DIM], dxa[DIM];

  // 1. apply a rotation around an in the x-y plane
  //    to satisfy the constraints for z-coordinates

  // build the coordinates
  vec_diff(xab0, x0[ID_A], x0[ID_B]);
  vec_diff(xac0, x0[ID_A], x0[ID_C]);
  vec_cross(xz0, xab0, xac0);
  dz0 = vec_norm(xz0);
  vec_vsmul(z0, xz0, 1.0/dz0);
  //printf("xz0 %g %g %g\n", xz0[0], xz0[1], xz0[2]);

  vec_diff(xab1, x1[ID_A], x1[ID_B]);
  zab = vec_dot(xab1, z0);
  vec_diff(xac1, x1[ID_A], x1[ID_C]);
  zac = vec_dot(xac1, z0);

  // vaxis: unnormalized axis of the first rotation
  vec_lincomb(vaxis, xab0, xac0, zac, -zab);
  //printf("vaxis %g %g %g; zab %g, zac %g\n", vaxis[0], vaxis[1], vaxis[2], zab, zac);
  daxis2 = vec_sqr(vaxis);
  if (daxis2 < 1e-12) { // zab == zac == 0
    // the first rotation is unnecessary
    vec_copy(xab4, xab0);
    vec_copy(xac4, xac0);
  } else {
    daxis = sqrt(daxis2);
    sinth = daxis/dz0;
    sinth2 = sinth * sinth;
    one_minus_costh = sinth2/(1 + sqrt(1-sinth2));
    
    // compute the rotated xab
    vec_vsmul(xab0para, vaxis, vec_dot(vaxis, xab0)/daxis2);
    vec_diff(xab0perp, xab0, xab0para);
  
    vec_lincomb(dxab0, xab0perp, z0, -one_minus_costh, zab);
    vec_add(xab4, xab0, dxab0);

    // compute the rotated xac
    vec_vsmul(xac0para, vaxis, vec_dot(vaxis, xac0)/daxis2);
    vec_diff(xac0perp, xac0, xac0para);
  
    vec_lincomb(dxac0, xac0perp, z0, -one_minus_costh, zac);
    vec_add(xac4, xac0, dxac0);
  }

  // 2. apply a rotation around the z-axis
  //
  // Try to find a rotation such that
  //   mBO xB0 x (R xBA4 - xBA1)
  // + mCO xC0 x (R xCA4 - xCA1) = 0
  // where O is the center of mass
  //
  // Assuming mB == mC
  //
  // xBO0 x (R - 1) xBA4 + xCO0 x (R - 1) xCA4 = xBO0 x (xBA1 - xBA4) + xCO0 x (xCA1 - xCA4)
  
  // compute the coordinates of the center of mass
  //    xa = xo + frac_mass_h * (xa - xb) + frac_mass_h * (xa - xc)
  // so
  //    xo = xa - frac_mass_h * (xab + xac)
  vec_add(tmpsum, xab0, xac0);
  vec_sadd(xcom0, x0[ID_A], tmpsum, -frac_mass_h);

  // compute the coordinates relative to the center of mass
  vec_diff(xbo0, x0[ID_B], xcom0);
  vec_diff(xco0, x0[ID_C], xcom0);

  // compute the right-hand side of the equation
  vec_diff(dxab, xab1, xab4);
  vec_cross(tmpab, xbo0, dxab);

  vec_diff(dxac, xac1, xac4);
  vec_cross(tmpac, xco0, dxac);

  vec_add(tmpsum, tmpab, tmpac);
  rhs = vec_dot(tmpsum, z0);

  // compute the left-hand side of the equation
  
  // compute the cosine part
  vec_cross(tmpab, xbo0, xab4);
  vec_cross(tmpac, xco0, xac4);
  vec_add(tmpsum, tmpab, tmpac);
  a = vec_dot(tmpsum, z0);

  // compute the sine part
  b = vec_dot(xbo0, xab4)
    + vec_dot(xco0, xac4);

  // solve the equation of angle
  //
  // a * [cos(theta) - 1] + b * sin(theta) = rhs
  //
  // solving the quadratic equation for sin(theta)
  gam = rhs + a;
  a2b2 = a*a + b*b;
  //printf("gam %g\n", gam);
  sinth = (b*gam - a*sqrt(a2b2 - gam*gam)) / a2b2;
  sinth2 = sinth * sinth;
  one_minus_costh = sinth2/(1 + sqrt(1 - sinth2));

  // apply the second the rotation to xab and xac
  vec_vsmul(xab4para, z0, vec_dot(z0, xab4));
  vec_diff(xab4perp, xab4, xab4para);
  vec_cross(yab4, z0, xab4perp);
  vec_lincomb(dxab4, xab4perp, yab4, -one_minus_costh, sinth);
  vec_add(xab5, xab4, dxab4);
  
  vec_vsmul(xac4para, z0, vec_dot(z0, xac4));
  vec_diff(xac4perp, xac4, xac4para);
  vec_cross(yac4, z0, xac4perp);
  vec_lincomb(dxac4, xac4perp, yac4, -one_minus_costh, sinth);
  vec_add(xac5, xac4, dxac4);

  // output the final coordinates
  // compute the coordinates relative to the oxygen atom
  // from x1 to x5, the center of mass is unchanged
  //
  // xo = xa1 - frac_mass_h * (xab1 + xac1)
  // xo = xa5 - frac_mass_h * (xab5 + xac5)
  //
  // xa5 = xa1 + frac_mass_h * ((xab5 + xac5) - (xab1 + xac1))
  vec_add(x5sum, xab5, xac5);
  vec_add(x1sum, xab1, xac1);
  vec_diff(dxa, x5sum, x1sum);
  vec_sadd(xout[ID_A], x1[ID_A], dxa, frac_mass_h);
  vec_diff(xout[ID_B], xout[ID_A], xab5);
  vec_diff(xout[ID_C], xout[ID_A], xac5);
}


/* SETTLE algorithm for coordinates */
void constraint_settle_coordinates(constraint_t *cs)
{
  constraint_param_t *csp = cs->param;
  constraint_settle_data_t *sed = (constraint_settle_data_t *) cs->data->idata;
  int i;
  real (*x0)[DIM], (*x1)[DIM];

  for (i = 0; i < sed->nmol; i++) {
    x0 = csp->x0 + sed->id_o[i];
    x1 = csp->x1 + sed->id_o[i];
    constraint_settle_one_coordinates(sed, x0, x1, x1, sed->frac_mass_h); 
  }
}


void constraint_settle_one_velocities(constraint_settle_data_t *sed,
    real (*x)[DIM], real (*v0)[DIM], real (*v)[DIM])
{
  int i, j, ipr;
  real dv[CONSTRAINT_SETTLE_NATOM][DIM] = {{0}};
  real dxpr[CONSTRAINT_SETTLE_NPAIR][DIM],
       dvpr[CONSTRAINT_SETTLE_NPAIR][DIM],
       xvpr[CONSTRAINT_SETTLE_NPAIR];

  for (ipr = 0; ipr < CONSTRAINT_SETTLE_NPAIR; ipr++) {
    i = constraint_settle_pair_ids_[ipr][0];
    j = constraint_settle_pair_ids_[ipr][1];
    vec_diff(dxpr[ipr], x[i], x[j]);
    vec_diff(dvpr[ipr], v0[i], v0[j]);
    xvpr[ipr] = vec_dot(dxpr[ipr], dvpr[ipr]);
  }

  for (ipr = 0; ipr < CONSTRAINT_SETTLE_NPAIR; ipr++) {
    real lambda = vec_dot(sed->vinvmat[ipr], xvpr);
    //printf("xvpr %g %g %g; lambda %g\n", xvpr[0], xvpr[1], xvpr[2], lambda);
    i = constraint_settle_pair_ids_[ipr][0];
    j = constraint_settle_pair_ids_[ipr][1];
    vec_sinc(dv[i], dxpr[ipr], +lambda);
    vec_sinc(dv[j], dxpr[ipr], -lambda);
    //printf("%g %g %g; lambda %g\n", dv[i][0], dv[i][1], dv[i][2], lambda);
  }

  for (i = 0; i < CONSTRAINT_SETTLE_NATOM; i++) {
    vec_sadd(v[i], v0[i], dv[i], sed->invmass[i]);
  }
}


/* SETTLE algorithm for velocities */
void constraint_settle_velocities(constraint_t *cs)
{
  constraint_param_t *csp = cs->param;
  constraint_settle_data_t *sed = (constraint_settle_data_t *) cs->data->idata;
  int i;
  real (*x1)[DIM], (*v)[DIM];

  for (i = 0; i < sed->nmol; i++) {
    x1 = csp->x1 + sed->id_o[i];
    v = csp->v + sed->id_o[i];
    constraint_settle_one_velocities(sed, x1, v, v);
  }
}


int constraint_settle_apply(constraint_t *cs, unsigned flags)
{
  if (flags & CONSTRAINT_APPLY_COORDINATES) {
    constraint_settle_coordinates(cs);
  }

  if (flags & CONSTRAINT_APPLY_VELOCITIES) {
    constraint_settle_velocities(cs);
  }

  return 0;
}



constraint_t *constraint_new(int, constraint_param_t *);

constraint_t *constraint_settle_water_new_ez(int water_model,
    int nmol, int *id_o,
    real (*x0)[DIM], real (*x1)[DIM], real (*v)[DIM])
{
  constraint_settle_model_param_t *model_param;
 
  if (water_model == CONSTRAINT_SETTLE_WATER_SPC) {
    model_param = &constraint_settle_water_spc_param;
  } else if (water_model == CONSTRAINT_SETTLE_WATER_TIP3P) {
    model_param = &constraint_settle_water_tip3p_param;
  }

  constraint_settle_param_t settle_param = {0};
  memcpy(&(settle_param.model_param), model_param, sizeof(*model_param));
  settle_param.nmol = nmol;
  settle_param.id_o = id_o;

  constraint_param_t csp = {
    .n = nmol * 3,
    .mass = NULL,
    .x0 = x0,
    .x1 = x1,
    .v = v,
    .iparam = &settle_param,
  };

  constraint_t *cs = constraint_new(CONSTRAINT_TYPE_SETTLE, &csp);

  return cs;
}



#endif /* CONSTRAINT_SETTLE_H__ */

