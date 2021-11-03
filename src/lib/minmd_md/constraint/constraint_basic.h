#ifndef CONSTRAINT_BASIC_H__
#define CONSTRAINT_BASIC_H__

#include "def.h"

typedef struct {
  int n; /* number of atoms */
  const real *mass; /* mass of atoms */
  real (*x0)[DIM]; /* original unperturbed coordinates that satisify constraints */
  real (*x1)[DIM]; /* perturbed coordinates that do not satisfy constraints */
  real (*xc)[DIM]; /* corrected coordinates, can be the same as x1 */
  real (*v1)[DIM]; /* original unperturbed velocities */
  real (*vc)[DIM]; /* corrected velocities */
  void *algo_param; /* algorithm-specific parameters; */
} constraint_param_t;


enum {
  CONSTRAINT_TYPE_NULL = 0,
  CONSTRAINT_TYPE_SHAKE = 1,
  /* CONSTRAINT_TYPE_SETTLE, */
  CONSTRAINT_TYPE_COUNT
};

typedef struct {
  void *algo_data; /* algorithm-specific data */
} constraint_data_t;

typedef struct {
  int type;
  constraint_param_t *param;
  constraint_data_t *data;
} constraint_t;


typedef struct {
  int n;
  int (*id)[2]; /* the indicies of the two involved atoms */
  double *dist; /* distance between the two atoms */
} constraint_pairlist_t;


INLINE constraint_pairlist_t *constraint_pairlist_new(int n)
{
  constraint_pairlist_t *pl;

  XNEW(pl, 1);
  pl->n = n;
  XNEW(pl->id, n);
  XNEW(pl->dist, n);
  return pl;
}


INLINE void constraint_pairlist_delete(constraint_pairlist_t *pr)
{
  free(pl->id);
  free(pl->dist);
  free(pl);
}



#endif /* CONSTRAINT_BASIC_H__ */
