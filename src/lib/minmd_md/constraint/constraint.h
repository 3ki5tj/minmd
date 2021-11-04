#ifndef CONSTRAINT_H__
#define CONSTRAINT_H__


#include "constraint_shake.h"


constraint_t *constraint_new(int type, constraint_param_t *param)
{
  if (type == CONSTRAINT_TYPE_SHAKE) {
    return constraint_shake_new(param);
  } else {
    fprintf(stderr, "Error: unknown type of constraint solver %d\n", type);
    exit(1);
  }

  return NULL;
}


void constraint_delete(constraint_t *cs)
{
  if (cs->type == CONSTRAINT_TYPE_SHAKE) {
    constraint_shake_delete(cs);
  }
}


int constraint_apply(constraint_t *cs, unsigned flags)
{
  if (cs->type == CONSTRAINT_TYPE_SHAKE) {
    return constraint_shake_apply(cs, flags);
  }

  return 0;
}


#endif /* CONSTRAINT_H__ */

