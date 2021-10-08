#include <iostream>
#include <cmath>
#include <ctime>
#include "math.h"

int main(int argc, char **argv)
{
  real a = 1000.f, b = 0.f;
#ifdef MINMD_MYMATH
  real val = softplus(a, b);
#else
  real val = log(exp(a) + exp(b));
#endif
  std::cout << val << std::endl;

  rng_t *r = rng_init(0, time(NULL));
  std::cout << "random number: " << rng_gauss(r) << std::endl;
  rng_free(r);

  if (argc <= 1) {
    std::cout << argv[0] << " Version " << MINMD_VERSION_MAJOR << "."
      << MINMD_VERSION_MINOR << std::endl;
    std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
    return 1;
  }
  return 0;
}
