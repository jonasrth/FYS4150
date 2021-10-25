#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle
{
public:
  double q_;    // Charge of particle in elementary charge units
  double m_;    // Mass of particle in atomic mass units
  arma::vec r_; // Position vector of particle in micrometers
  arma::vec v_; // Velocity vector of particle in micrometers/microsecond

  // Constructor
  Particle(double q, double m, arma::vec r, arma::vec v);
};

#endif
