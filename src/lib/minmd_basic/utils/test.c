#include <stdio.h>
#include "utils.h"

int main(void)
{
  real (*x)[2];
  int i;
  const int n = 10;

  XNEW(x, n);
  for (i = 0; i < n; i++) {
    x[i][0] = i;
    x[i][1] = 1.0*i*i;
  }

  for (i = 0; i < n; i++) {
    printf("%d %g %g\n", i, x[i][0], x[i][1]);
  }

  return 0;
}
