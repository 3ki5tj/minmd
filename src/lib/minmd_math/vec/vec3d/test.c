#include "vec3d.h"

void test_rng_dir(rng_t *rng)
{
  real x[DIM];
  vec_rand_dir(x, rng);
  printf("norm = %g\n", x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

int main(void)
{
  rng_t *rng = rng_init(0, time(NULL));
  test_rng_dir(rng);
  rng_free(rng);
  return 0;
}
