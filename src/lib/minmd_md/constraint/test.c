#include "constraint.h"



int main(void)
{
  int n = 10, i;
  real *mass;
  real (*v)[DIM];

  /* allocate memory for masses and velocities */
  XNEW(mass, n);
  XNEW(v, n);
  for (i = 0; i < n; i++) {
    mass[i] = 1;
  }

  free(v);
  free(mass);
  return 0;
}
