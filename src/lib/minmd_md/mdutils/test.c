#include "mdutils.h"

int main(void)
{
  int n = 27;
  real x[100][DIM];
  rng_t *rng;

  rng = rng_init(0, 0);
  mdutils_init_face_centered_lattice(n, 3.0, x, rng);
  rng_free(rng);

  return 0;
}
