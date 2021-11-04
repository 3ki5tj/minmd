#ifndef EWALD_BASIC_H__
#define EWALD_BASIC_H__


#include "def.h"
#include "utils.h"
#include "vec.h"

/* generic Ewald parameters */
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
  void *iparam;
} ewald_param_t;

/* generic Ewald data */
#include "ewald_data.h"

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



#endif /* EWALD_BASIC_H__ */

