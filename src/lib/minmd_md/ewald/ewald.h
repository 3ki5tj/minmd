#ifndef EWALD_H__
#define EWALD_H__


#include "ewald_direct.h"

ewald_t *ewald_init(int type, ewald_param_t *param)
{
  if (type == EWALD_TYPE_DIRECT) {
    return ewald_direct_init(param);
  } else {
    fprintf(stderr, "Error: ewald_init() does not support type %d\n", type);
    return NULL;
  }
}


real ewald_force(ewald_t *ew, ewald_force_options_t *ewf_opt)
{
  if (ew->type == EWALD_TYPE_DIRECT) {
    return ewald_direct_force(ew, ewf_opt);
  } else {
    fprintf(stderr, "Error: ewald_force() does not support type %d\n", ew->type);
    return 0;
  }
}


real ewald_free(ewald_t *ew)
{
  if (ew->type == EWALD_TYPE_DIRECT) {
    ewald_direct_free(ew);
  } else {
    fprintf(stderr, "Error: ewald_free() does not support type %d\n", ew->type);
  }
}


#endif /* EWALD_H__ */

