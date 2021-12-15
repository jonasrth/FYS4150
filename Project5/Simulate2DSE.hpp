#ifndef __Simulate2DSE_hpp__
#define __Simulate2DSE_hpp__

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <armadillo>

class Simulate2DSE
{
public:

  // Parameters for simulation
  int M_;        // Size of simulated matrix
  double h_;     // Step size in x- and y direction
  double dt_;    // Step size in time

  // Parameters defining gaussian wave packet
  double x_c_, y_c_;          // x- and y position
  double sigma_x_, sigma_y_;  // x- and y standard deviation
  double p_x_, p_y_;          // x- and y momentum

  // Parameters for potential
  std::string slit_;   // "single", "double" or "triple" slit
  double v0_;         // Barrier potential size
  arma::mat v_;       // Matrix holding the potential

  // For solving PDEs
  arma::cx_double r_;     // Constants involved in PDE matrix equation
  arma::sp_cx_mat A_,B_;  // Matrices involved in PDE matrix equations

  // States of system
  arma::cx_mat U_init_; // Complex matrix representing initial state of 2D SE
  arma::cx_vec u_;      // Complex vector representing current state of 2D SE
  arma::cx_cube U_;     // Complex cube representing all states of 2D SE



  // Constructor
  // -------------
  // Parameters:
  // h                - [1] Spacial step size
  // dt               - [1] Time step size
  // x_c, y_c         - [1] x and y position of initial wave packet
  // sigma_x, sigma_y - [1] x and y standard deviation of initial wave packet
  // p_x, p_y         - [1] x and y momentum of initial wave packet
  // v0               - [1] Size of barrier potential
  // slit             - [1] "single", "double" or "triple",
  //                        defines number of slits in the barrier potential
  Simulate2DSE(double h, double dt, double x_c, double y_c,
                double sigma_x, double sigma_y, double p_x, double p_y,
                double v0, std::string slit = "double");

  // Converts matrix indices to vector indices
  int ij_to_k(int i, int j);

  // Initialises inital state of 2D SE
  void initialise_u(arma::cx_mat& u);

  // Initialises barrier potential
  void initialise_v(arma::mat& v);

  // Converts vector a to matrix A
  void vec_to_mat(arma::cx_mat& A, arma::cx_vec a);

  // Converts matrix A to vector a
  void mat_to_vec(arma::cx_vec& a, arma::cx_mat A);

  // Fills matrix for solving the 2D SE, with diagonal a and remaining elements r
  void fill_matrix(arma::sp_cx_mat& A, arma::cx_vec a, arma::cx_double r);

  // Makes matrices A and B for solving 2D SE, using fill_matrix function'.
  void make_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B);

  // Simulates and saves all states of the 2D SE for a time period T.
  void simulate(double T);

};

#endif
