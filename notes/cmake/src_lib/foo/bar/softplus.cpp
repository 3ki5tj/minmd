#include "softplus.h"

/* log(exp(a) + exp(b)) */
real softplus(real a, real b)
{
  if (a < b) {
    real c = a; a = b; b = c;
  }
  return a + log(1 + exp(b-a));
}
