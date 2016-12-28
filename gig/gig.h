#ifndef GIG_H
#define GIG_H

#include <random>

class Random
{
private:
  std::default_random_engine generator;

  double normal(void)
  {
    return std::normal_distribution<double>()(generator);
  }

  double uniform(void)
  {
    return std::uniform_real_distribution<double>()(generator);
  }

  double gamma(double shape, double scale)
  {
    return std::gamma_distribution<double>(shape, scale)(generator);
  }

public:
  Random(unsigned int seed) : generator(seed) {}
  double gig(double lambda, double chi, double psi);
};

#endif
