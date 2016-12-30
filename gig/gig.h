#ifndef GIG_H
#define GIG_H

#include <random>

template< class Generator >
class Random {
private:
  // std::default_random_engine generator;
  Generator generator;

  double normal(void) { return std::normal_distribution<double>()(generator); }

  double uniform(void) {
    return std::uniform_real_distribution<double>()(generator);
  }

  double gamma(double shape, double scale) {
    return std::gamma_distribution<double>(shape, scale)(generator);
  }

  double ratio_of_uniforms_noshift(double lambda, double lambda_old,
                                   double omega, double alpha);
  double unnamed_approach(double lambda, double lambda_old, double omega,
                         double alpha);
  double ratio_of_uniforms_mode(double lambda, double lambda_old, double omega,
                                double alpha);

public:
  // template< class Generator >
  // result_type
  // Generator& g
  Random(Generator& g) : generator(g) {}
  Random(unsigned int seed);
  double gig(double lambda, double chi, double psi);
};

template< >
Random<std::default_random_engine>::Random(unsigned int seed)
  : generator(seed)
{
}

#endif
