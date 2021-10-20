#include <time.h>
#include "rng.h"

void test_rng_gauss(rng_t *rng)
{
  int i;
  const int n = 10000000;
  double x, sx = 0, sxx = 0;

  for (i = 0; i < n; i++) {
    x = rng_gauss(rng);
    sx += x;
    sxx += x * x;
  }
  x = sx / n;
  printf("rng_gauss: %g +/- %g\n", x, sqrt(sxx/n - x*x));
}

void test_rng_chisqr(rng_t *rng)
{
  int i, j;
  const int n = 1000000, k = 10;
  double x, sx = 0, sxx = 0;
  double y, sy = 0, syy = 0;

  for (i = 0; i < n; i++) {
    x = rng_chisqr(rng, k);
    sx += x;
    sxx += x * x;

    for (y = 0, j = 0; j < k; j++) {
      x = rng_gauss(rng);
      y += x*x;
    }
    sy += y;
    syy += y*y;
  }
  x = sx / n;
  y = sy / n;
  printf("rng_chisqr: %g +/- %g\n"
         "   compare: %g +/- %g\n",
      x, sqrt(sxx/n - x*x),
      y, sqrt(syy/n - y*y));
}

int main(void)
{
  rng_t *rng = rng_init(0, time(NULL));

  test_rng_gauss(rng);
  test_rng_chisqr(rng);

  rng_free(rng);
  return 0;
}
