#include <iostream>
#include "softplus.h"

int main(void)
{
  real a = 3, b = 4;
  std::cout << softplus(a, b) << std::endl;
  return 0;
}
