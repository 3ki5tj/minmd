#ifndef CONSTRAINT_BASIC_H__
#define CONSTRAINT_BASIC_H__


#include "def.h"
#include "constraint_pairlist.h"


typedef struct {
  int n; /* number of atoms */
  const real *mass; /* mass of atoms */
  real (*x0)[DIM]; /* original unperturbed coordinates that satisify constraints */
  real (*x1)[DIM]; /* perturbed coordinates that do not satisfy constraints */
  real (*xc)[DIM]; /* corrected coordinates, can be the same as x1 */
  real (*v1)[DIM]; /* original unperturbed velocities */
  real (*vc)[DIM]; /* corrected velocities */
  void *iparam; /* algorithm-specific parameters; */
} constraint_param_t;


enum {
  CONSTRAINT_TYPE_NULL = 0,
  CONSTRAINT_TYPE_SHAKE = 1,
  /* CONSTRAINT_TYPE_SETTLE, */
  CONSTRAINT_TYPE_COUNT
};


typedef struct {
  /* currently no general data */
  void *idata; /* algorithm-specific data */
} constraint_data_t;

typedef struct {
  int type;
  constraint_param_t *param;
  constraint_data_t *data;
} constraint_t;


typedef struct {
  int type;
  constraint_param_t *param;
  constraint_data_t *data;
} constraint_t;



#endif /* CONSTRAINT_BASIC_H__ */
