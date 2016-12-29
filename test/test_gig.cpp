#include "gig/gig.h"
#include <iostream>
#include <iomanip>
#include <cassert>

using std::abs;

bool isclose(double x, double y) {
  double rtol = 1e-7;
  double atol = 0;
  return abs(x - y) <= atol + rtol * abs(y);
}

int main(int argc, char const *argv[]) {

  Random random(0);

  assert(isclose(random.gig(2.1, 0.1, 1.0), 1.30869321355819901));

  return 0;
}
