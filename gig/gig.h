#ifndef GIG_H
#define GIG_H

#include <random>

using std::default_random_engine;

class Random {
private:
  default_random_engine generator;

public:
  Random(unsigned int seed) : generator(seed) {}
  Random(default_random_engine& generator) : generator(generator) {}
  double gig(double lambda, double chi, double psi);
};

#endif
