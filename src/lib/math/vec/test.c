#include <stdio.h>
#include "vec.h"

int main(void)
{
  real v1[DIM] = {1, 0, 0};
  real v2[DIM] = {0, 1, 0};
  real v3[DIM];

  vec_cross(v3, v1, v2);
  printf("v1 x v2 = (%g, %g, %g)\n", v3[0], v3[1], v3[2]);
  return 0;
}
