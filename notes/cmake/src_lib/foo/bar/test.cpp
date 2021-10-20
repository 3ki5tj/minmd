#include <iostream>
#include "bar.h"

int main(void)
{
  real a = 3, b = 4;
  std::cout << softplus(a, b) << std::endl;
  return 0;
}
