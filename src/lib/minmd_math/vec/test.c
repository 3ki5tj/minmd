#include <stdio.h>
#include "vec.h"

void test_cross(void)
{
  real v1[DIM] = {1, 0, 0};
  real v2[DIM] = {0, 1, 0};
  real v3[DIM];

  vec_cross(v3, v1, v2);
  printf("v1 x v2 = (%g, %g, %g)\n", v3[0], v3[1], v3[2]);
}

void test_rng_dir(rng_t *rng)
{
  real x[DIM];
  vec_rand_dir(x, rng);
  printf("norm = %g\n", x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

int main(void)
{
  rng_t *rng = rng_init(0, time(NULL));

  test_cross();
  test_rng_dir(rng);

  rng_free(rng);
  return 0;
}
