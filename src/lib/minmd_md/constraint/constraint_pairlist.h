#ifndef CONSTRAINT_PAIRLIST_H__
#define CONSTRAINT_PAIRLIST_H__


#include "def.h"


typedef struct {
  int n;
  int (*ids)[2]; /* the indicies of the two involved atoms */
  double *dist; /* distance between the two atoms */
} constraint_pairlist_t;


INLINE constraint_pairlist_t *constraint_pairlist_new(int n)
{
  constraint_pairlist_t *pl;

  XNEW(pl, 1);
  pl->n = n;
  XNEW(pl->ids, n);
  XNEW(pl->dist, n);
  return pl;
}


INLINE constraint_pairlist_t *constraint_pairlist_clone(constraint_pairlist_t *src)
{
  constraint_pairlist_t *dest = constraint_pairlist_new(src->n);

  /* copy the pair lists */
  memmove(dest->ids, src->ids, src->n*sizeof(src->id[0]));
  memmove(dest->dist, src->dist, src->n*sizeof(src->dist[0]));

  return dest;
}


INLINE void constraint_pairlist_delete(constraint_pairlist_t *pr)
{
  free(pl->id);
  free(pl->dist);
  free(pl);
}


#endif /* CONSTRAINT_PAIRLIST_H__ */
