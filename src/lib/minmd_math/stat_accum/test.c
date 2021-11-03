#include "stat_accum.h"
#include "rng.h"

int main(void)
{
  stat_accum_t *acc;
  rng_t *rng;
  int n = 10000000, i;
  double x, xmean, xsig;

  rng = rng_new(0, 0);
  acc = stat_accum_new();
  for (i = 0; i < n; i++) {
    x = rng_gauss(rng);
    stat_accum_add(acc, x);
  }
  xsig = stat_accum_get_std(acc, &xmean);
  printf("%g samples, mean %g, std %g\n", 1.*acc->count, xmean, xsig);
  stat_accum_delete(acc);
  rng_delete(rng);
  return 0;
}
