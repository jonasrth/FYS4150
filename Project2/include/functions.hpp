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
#include <armadillo>


// Function declarations:

// Create a tridiagonal matrix tridiag(a,d,e) of size N*N,
// from scalar input a, d, and e. That is, create a matrix where
// - all N-1 elements on the subdiagonal have value a
// - all N elements on the diagonal have value d
// - all N-1 elements on the superdiagonal have value e
arma::mat create_tridiagonal(int N, double a, double d, double e);

// Returns the analytical eigenvalues to the solution of our symmetric
// tridiagonal matrix equation describing the second order differential equation
// for a special case of a buckling beam. The matrix has size (N*N) and
// signature tridiag(a,d,a).
arma::vec analytical_eigenvalues(int N, double a, double d);

// Returns the analytical eigenvectors to the solution of our symmetric
// tridiagonal matrix equation describing the second order differential equation
// for a special case of a buckling beam. The matrix has size (N*N) and
// signature tridiag(a,d,a).
arma::mat analytical_eigenvectors(int N, double a, double d);

// Determine the max off-diagonal element of a symmetric matrix A
// - Saves the matrix element indicies to k and l
// - Returns absolute value of A(k,l) as the function return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues,
                        arma::mat& eigenvectors, const int maxiter,
                        int& iterations, bool& converged);

// Writes vector x and matrix y to file with name "filename"
// - x must have same amount of elements as y has rows.
void write_to_file(std::string filename, arma::vec x, arma::mat y);

// Writes x as well as the first three eigenvectors to file for the
// analytical- and numerical solution of the buckling beam second order
// differential equation.
// - N is the length of the solution
// - eps and maxiter are the same as in the jacobi_eigensolver function
void compare_jacobi_analytical(int N, double eps, const int maxiter);



#endif
