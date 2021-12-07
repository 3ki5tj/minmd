#ifndef CONSTRAINT_H__
#define CONSTRAINT_H__


#include "constraint_shake.h"
#include "constraint_settle.h"


constraint_t *constraint_new(int type, constraint_param_t *param)
{
  if (type == CONSTRAINT_TYPE_SHAKE) {
    return constraint_shake_new(param);
  } else if (type == CONSTRAINT_TYPE_SETTLE) {
    return constraint_settle_new(param);
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
  } else if (cs->type == CONSTRAINT_TYPE_SETTLE) {
    constraint_settle_delete(cs);
  }
}


int constraint_apply(constraint_t *cs, unsigned flags)
{
  if (cs->type == CONSTRAINT_TYPE_SHAKE) {
    return constraint_shake_apply(cs, flags);
  } else if (cs->type == CONSTRAINT_TYPE_SETTLE) {
    return constraint_settle_apply(cs, flags);
  }

  return 0;
}



#endif /* CONSTRAINT_H__ */

