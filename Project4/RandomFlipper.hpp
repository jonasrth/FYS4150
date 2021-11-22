#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <armadillo>
#include <chrono>

class RandomFlipper
{
public:
  // Member variables that hold a <random> generator
  // and distribution (uniform distribution [0,1)).
  std::mt19937 generator;
  std::uniform_real_distribution<double> uniform_dist = std::uniform_real_distribution<double>(0.0 ,1.0);

  // Current state
  arma::mat S_;  // Square matrix representing lattice with spins +-1
  int L_;     // Length of matrix S

  double E_;      // Total energy of lattice in coupling constants J
  double M_;      // Total magnetisation of lattice
  double T_;      // Temperature of lattice in J/k_B
  arma::vec w_;   // Stores Boltzmann factors to prevent solving exp(...) at each step

  int current_state = 0;

  // Constructor
  RandomFlipper(arma::mat S, double T);

  double find_energy();

  void update_state();

  void monte_carlo_cycle();

};

#endif
