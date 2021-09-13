//
// Including headers and declaring useful functions
//

// Ensuring header file is only included once for each compilation
#ifndef __functions_hpp__
#define __functions_hpp__

// Headers:
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>


// Function declarations:


// Returns u(x) for a given x
double u_func(double x);

// Returns f(x) for a given x
double f_func(double x);

// Function estimating the values of u(x) for n equally spaced x-values between
// x_init and x_end using an algorithm for solving the Poisson equation
// as a matrix equation with a general tridiagonal matrix
std::vector<std::vector<double> > approximate_poisson_general(double (func)(double),
                                  double u_init, double u_end,
                                  double x_init, double x_end, int steps);

// Same as function above, but for a special tridiagonal matrix with signature
// (-1,2,-1)
std::vector<std::vector<double> > approximate_poisson_special(double (func)(double),
                                  double u_init, double u_end,
                                  double x_init, double x_end, int steps);

// Estimates absolute error between estimated u(x) and approximated u(x),
// i.e. v(x)
std::vector<std::vector<double> > abs_error(std::vector<double> x,
                                  std::vector<double> v, double (u_func)(double));

// Estimates relative error between estimated u(x) and approximated u(x),
// i.e. v(x)
std::vector<std::vector<double> > rel_error(std::vector<double> x,
                                  std::vector<double> v, double (u_func)(double));

// Loads values of two vectors to a text file so it can be read by
// a python program
void load_to_file(std::string filename, std::vector<double> x, std::vector<double> y);




#endif
