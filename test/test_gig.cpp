#include "gig/gig.h"
#include <cassert>

using std::abs;

bool isclose(double x, double y) {
  double rtol = 1e-7;
  double atol = 0;
  return abs(x - y) <= atol + rtol * abs(y);
}

void test_seed()
{
  Random random(0);

  assert(isclose(random.gig(2.1, 0.1, 1.0), 1.30869321355819901));
  assert(isclose(random.gig(2.1, 0.1, 10.4) , 0.156671928679679994811380084));
  assert(isclose(random.gig(2.1, 0.1, 1003.2), 0.007932940552363928338186481));
  assert(isclose(random.gig(0.1, 10.1, 0.001), 811.3064195112882543980958872));
  assert(isclose(random.gig(2.2, 0.001, 0.4), 16.09480769605515959597141773));
  assert(isclose(random.gig(-1, 1, 2), 0.79130783712979646526974875087));
}

void test_generator()
{
  std::mt19937_64 generator(0);
  Random random(generator);

  assert(isclose(random.gig(2.1, 0.1, 1.0), 1.30869321355819901));
  assert(isclose(random.gig(2.1, 0.1, 10.4) , 0.156671928679679994811380084));
  assert(isclose(random.gig(2.1, 0.1, 1003.2), 0.007932940552363928338186481));
  assert(isclose(random.gig(0.1, 10.1, 0.001), 811.3064195112882543980958872));
  assert(isclose(random.gig(2.2, 0.001, 0.4), 16.09480769605515959597141773));
  assert(isclose(random.gig(-1, 1, 2), 0.79130783712979646526974875087));
}

int main(int argc, char const *argv[]) {

  test_seed();
  test_generator();

  return 0;
}
