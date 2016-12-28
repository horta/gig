#ifndef GIG_H
#define GIG_H

#include <random>

class Random
{
private:
  std::default_random_engine generator;

  double normal(void);
  double uniform(void);
  double gamma(double shape, double scale);

public:
  Random(unsigned int seed) : generator(seed) {}
  double gig(double lambda, double chi, double psi);
};

#endif
