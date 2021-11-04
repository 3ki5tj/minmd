#ifndef SHAKE_H__
#define SHAKE_H__


typedef struct {
  int itmax; /* maximum number of iterations */
  constraint_pairlist_t *pr;
} constraint_shake_param_t;

constraint_shake_param_clone(constraint_shake_param_t *dest, constraint_shake_param_t *src)
{
  constraint_shake_param_t *dest; 

  XCLONE(dest, src, sizeof(*src)); /* allocate space and shallow copy */
  dest->pr = constraint_pairlist_clone(src->pr);
  return dest;
}

constraint_shake_param_delete(constraint_shake_param_t *shp)
{
  constraint_pairlist_delete(shp->pr);
  free(shp);
}


constraint_t *constraint_shake_new(constraint_param_t *param)
{
  constraint_t *cs;

  XNEW(cs, 1);
  cs->type = CONSTRAINT_TYPE_SHAKE;

  XCLONE(cs->param, param, sizeof(*param));
  XCLONE(cs->param->iparam, param->iparam, sizeof(constraint_shake_param_t));

  return cs;
}


void constraint_shake_delete(constraint_t *cs)
{
  constraint_shake_param_delete((constraint_shake_param_t *) cs->param->iparam);
  free(cs->param);
  // constraint_shake_data_delete(cs->data->idata);
  
  free(cs);
}


void constraint_shake_apply(constaint_t *cs)
{
}


#endif /* SHAKE_H__ */

