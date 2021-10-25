#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <armadillo>
#include "Particle.hpp"

class PenningTrap
{
public:
  double B0_;  // Magnetic field strength
  double V0_;  // Electric potential
  double d_;   // Characteristic dimension
  double f_;        // E_field oscillation amplitude
  double omegaV_;   // E-field oscillation frequency
  std::vector<Particle> Particles_;

  // Constructor
  // B0 given in Tesla, V0 given in Volts and d given in cm.
  PenningTrap(double B0, double V0, double d, double f = 0, double omegaV = 0);

  double t_ = 0;  // time in microseconds

  double k_e = 1.38935333e5;    // Coulomb constant in base units
  double T = 9.64852558e1;      // Tesla in base units
  double V = 9.64852558e7;      // Volt in base units
  double cm = 1e4;              // Centimeter in base units

  // Method for adding particle to PenningTrap
  void add_particle(Particle p);

  // Method for finding external electric field:
  arma::vec E_field(arma::vec r);

  // Method for finding external magnetic field:
  arma::vec B_field(arma::vec r);

  // Method for finding force from external field on particle i:
  arma::vec F_external(int i);

  // Method for finding force from particle j on particle i
  arma::vec F_particle(int i, int j);

  // Method for finding total force from other particles on particle i
  arma::vec F_particle_tot(int i);

  // Method for finding total force on particle i
  arma::vec F_tot(int i, bool interaction);

  // Method for evolving the system one time step (dt) using Runge-Kutta 4th order
  void evolve_RK4(double dt, bool interaction);

  // Method for evolving the system one time step (dt) using Forward Euler
  void evolve_forward_Euler(double dt, bool interaction);
};

#endif
