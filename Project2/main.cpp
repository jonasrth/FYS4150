
#include "functions.hpp"


int main(){

  // matrix size N, step size h, diagonal d, sub/super-diagonal a,
  int N;
  double h;
  double a,d;

  // matrix A, eigenvalues and eigenvectors
  arma::mat A;
  arma::vec eigval;
  arma::mat eigvec;

  // analytical eigenvalues and eigenvectors
  vec analytical_eigval;
  arma::mat analytical_eigvec;

  //
  // Problem 3:
  //

  std::cout << std::endl << "PROBLEM 3:" << std::endl << std::endl;

  N = 6;
  h = 1.0/(N+1);

  a = -1/(h*h);
  d = 2/(h*h);

  A = create_tridiagonal(N,a,d,a);

  std::cout << "Matrix A:" << std::endl << A << std::endl;

  // Finding analytical eigenvalues and eigenvectors:
  analytical_eigval = analytical_eigenvalues(N,a,d);
  analytical_eigvec = analytical_eigenvectors(N,a,d);

  // Finding eigenvalues and eigenvectors through armadillo:
  arma::eig_sym(eigval,eigvec,A);

  std::cout << "Analytical eigenvalues:" << std::endl << analytical_eigval << std::endl;
  std::cout << "Armadillo eigenvalues:" << std::endl << eigval << std::endl;
  std::cout << "Analytical eigenvectors (columns):" << std::endl << analytical_eigvec << std::endl;
  std::cout << "Armadillo eigenvectors (columns):" << std::endl << eigvec << std::endl;

  //
  // Problem 4b:
  //

  std::cout << std::endl << "PROBLEM 4b:" << std::endl << std::endl;

  arma::mat X(4,4,arma::fill::eye);
  A = X;

  A(0,3) = 0.5; A(3,0) = 0.5; A(1,2) = -0.7; A(2,1) = -0.7;

  std::cout << "Matrix A:" << std::endl << A << std::endl;

  int k,l;

  double A_max = max_offdiag_symmetric(A,k,l);

  std::cout << "A has max off-diagonal element of value " << A_max
            << " at indeces k = " << k << ", l = " << l << "." << std::endl;

  //
  // Problem 5:
  //

  std::cout << std::endl << std::endl << "PROBLEM 5b:" << std::endl << std::endl;

  A = create_tridiagonal(N,a,d,a);

  double eps = pow(10,-16);
  int maxiter = 100000;
  int iterations;
  bool converged;

  jacobi_eigensolver(A,eps,eigval,eigvec,maxiter,iterations,converged);

  std::cout << "Analytical eigenvalues:" << std::endl << analytical_eigval << std::endl;
  std::cout << "Jacobi eigenvalues:" << std::endl << eigval << std::endl;
  std::cout << "Analytical eigenvectors (columns):" << std::endl << analytical_eigvec << std::endl;
  std::cout << "Jacobi eigenvectors (columns):" << std::endl << eigvec << std::endl;

  //
  // Problem 6a:
  //

  int N_max = 50;

  arma::vec N_vals(N_max-1),iter_vals(N_max-1);

  for(int i=2; i<N_max+1; i++){
    A = create_tridiagonal(i,a,d,a);
    jacobi_eigensolver(A,eps,eigval,eigvec,maxiter,iterations,converged);

    N_vals(i-2) = i;
    iter_vals(i-2) = iterations;
  }

  write_to_file("transformations_tridiag.txt",N_vals,iter_vals);

  for(int i=2; i<N_max+1; i++){
    // Generate random N*N matrix
    arma::mat A = arma::mat(i, i).randn();

    // Symmetrize the matrix by reflecting the upper triangle to lower triangle
    A = arma::symmatu(A);

    jacobi_eigensolver(A,eps,eigval,eigvec,maxiter,iterations,converged);

    N_vals(i-2) = i;
    iter_vals(i-2) = iterations;
  }

  write_to_file("transformations_dense.txt",N_vals,iter_vals);

  //
  // Problem 7:
  //

  compare_jacobi_analytical(10,eps,maxiter);

  compare_jacobi_analytical(100,eps,maxiter);


  return 0;

}
