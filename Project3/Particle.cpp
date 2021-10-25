#include "Particle.hpp"

// Definition of the constructor
Particle::Particle(double q, double m, arma::vec r, arma::vec v)
{
  // Assign input values (q,m,r,v) to the class member variables (q_, m_, r_, v_)
  q_ = q;
  m_ = m;
  r_ = r;
  v_ = v;
}
