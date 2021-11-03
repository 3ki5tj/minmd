#include "mdutils.h"

int main(void)
{
  int n = 27;
  real x[100][DIM];
  rng_t *rng;

  rng = rng_new(0, 0);
  rng_delete(rng);

  return 0;
}
