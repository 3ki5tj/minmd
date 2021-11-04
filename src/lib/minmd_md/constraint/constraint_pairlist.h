#ifndef CONSTRAINT_PAIRLIST_H__
#define CONSTRAINT_PAIRLIST_H__


#include "utils.h"


typedef struct {
  int i; /* index of the first atom */
  int j; /* index of the second atom */
  real dist; /* distance between the two atoms */
} constraint_pair_param_t;


typedef struct {
  int n;
  constraint_pair_param_t *pr;
} constraint_pairlist_t;


INLINE constraint_pairlist_t *constraint_pairlist_new(int n)
{
  constraint_pairlist_t *pl;

  XNEW(pl, 1);
  pl->n = n;
  XNEW(pl->pr, n);
  return pl;
}


INLINE constraint_pairlist_t *constraint_pairlist_clone(constraint_pairlist_t *src)
{
  constraint_pairlist_t *dest = constraint_pairlist_new(src->n);

  /* copy the pair lists */
  memmove(dest->pr, src->pr, src->n*sizeof(src->pr[0]));

  return dest;
}


INLINE void constraint_pairlist_delete(constraint_pairlist_t *pl)
{
  free(pl->pr);
  free(pl);
}


#endif /* CONSTRAINT_PAIRLIST_H__ */

