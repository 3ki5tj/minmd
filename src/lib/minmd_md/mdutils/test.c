#include "mdutils.h"

int main(void)
{
  int n = 27;
  real x[100][DIM];
  rng_t *rng;

  rng = rng_init(0, 0);
  rng_free(rng);

  return 0;
}
