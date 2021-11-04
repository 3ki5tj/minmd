#ifndef CONSTRAINT_BASIC_H__
#define CONSTRAINT_BASIC_H__


#include "def.h"
#include "constraint_pairlist.h"


typedef struct {
  int n; /* number of atoms */
  const real *mass; /* mass of atoms */
  real (*x0)[DIM]; /* original unperturbed coordinates that satisify constraints */
  real (*x1)[DIM]; /* perturbed coordinates that do not satisfy constraints */
  real (*v)[DIM]; /* original unperturbed velocities */
  void *iparam; /* algorithm-specific parameters */
} constraint_param_t;


typedef struct {
  /* currently no general data */
  void *idata; /* algorithm-specific data */
} constraint_data_t;


INLINE constraint_data_t *constraint_data_new(constraint_param_t *csp, void *idata)
{
  constraint_data_t *csd;

  XNEW(csd, 1);
  csd->idata = idata;

  return csd;
}


INLINE void constraint_data_delete(constraint_data_t *csd)
{
  free(csd);
}



enum {
  CONSTRAINT_TYPE_NULL = 0,
  CONSTRAINT_TYPE_SHAKE = 1,
  /* CONSTRAINT_TYPE_SETTLE, */
  CONSTRAINT_TYPE_COUNT
};


typedef struct {
  int type;
  constraint_param_t *param;
  constraint_data_t *data;
} constraint_t;


/* apply flags */
enum {
  CONSTRAINT_APPLY_COORDINATE = 1,
  CONSTRAINT_APPLY_COORDINATES = 1,
  CONSTRAINT_APPLY_POSITION = 1,
  CONSTRAINT_APPLY_POSITIONS = 1,
  CONSTRAINT_APPLY_VELOCITY = 2,
  CONSTRAINT_APPLY_VELOCITIES = 2,
};

#endif /* CONSTRAINT_BASIC_H__ */
