#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <armadillo>
//#include <time.h>


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

void compare_jacobi_analytical(int N, double eps, const int maxiter);


int main(){


  // Problem 3:

  int N;
  double h;
  double a,d;

  N = 6;
  h = 1.0/(N+1);

  a = -1/(h*h);
  d = 2/(h*h);


  arma::mat A;
  A = create_tridiagonal(N,a,d,a);

  arma::vec eigval;
  arma::mat eigvec;

  arma::vec eigval_ = analytical_eigenvalues(N,a,d);
  arma::mat eigvec_ = analytical_eigenvectors(N,a,d);

  arma::eig_sym(eigval,eigvec,A);

  //std::cout << B << std::endl;

  //std::cout << eigvec << std::endl;
  //std::cout << eigvec_ << std::endl;

  /*
  // Problem 4a:

  arma::mat A(4,4,arma::fill::eye);
  arma::mat R(4,4,arma::fill::eye);

  A(0,3) = 0.5; A(3,0) = 0.5; A(1,2) = -0.7; A(2,1) = -0.7;

  std::cout << A << std::endl;

  int k,l;

  double A_max = max_offdiag_symmetric(A,k,l);

  std::cout << A_max << std::endl;
  std::cout << k << l << std::endl;

  jacobi_rotate(A,R,k,l);

  std::cout << A << std::endl;
  std::cout << R << std::endl;
  */

  // Problem 5:

  A = create_tridiagonal(N,a,d,a);

  double eps = pow(10,-16);
  int maxiter = 100000;
  int iterations;
  bool converged;

  arma::vec eigenvalues;
  arma::mat eigenvectors;

  jacobi_eigensolver(A,eps,eigenvalues,eigenvectors,maxiter,iterations,converged);

  std::cout << "Iterations: " << iterations << std::endl;

  std::cout << eigenvalues << std::endl;
  std::cout << eigenvectors << std::endl;

  //std::cout << eigenvalues(sort_index(eigenvalues)) << std::endl;

  /*
  // Problem 6a:


  for(int i=2; i<21; i++){
    A = create_tridiagonal(i,a,d,a);
    jacobi_eigensolver(A,eps,eigenvalues,eigenvectors,maxiter,iterations,converged);
    std::cout << iterations << std::endl;
  }

  for(int i=2; i<21; i++){
    // Generate random N*N matrix
    arma::mat A = arma::mat(i, i).randn();

    // Symmetrize the matrix by reflecting the upper triangle to lower triangle
    A = arma::symmatu(A);

    jacobi_eigensolver(A,eps,eigenvalues,eigenvectors,maxiter,iterations,converged);
    std::cout << iterations << std::endl;

  }

  */

  // Problem 7:

  compare_jacobi_analytical(10,eps,maxiter);

  compare_jacobi_analytical(100,eps,maxiter);


  return 0;

}

void compare_jacobi_analytical(int N, double eps, const int maxiter){

  double h = 1.0/(N+1);
  double a = -1/(h*h), d = 2/(h*h);

  arma::vec x = arma::linspace(0,1,N+2);

  // Finding analytical solutions and writing to file:

  arma::mat eigvec;

  eigvec = analytical_eigenvectors(N,a,d);

  eigvec.shed_cols(3,N-1);
  eigvec.insert_rows(0,1);
  eigvec.insert_rows(N+1,1);

  write_to_file("solution_N"+std::to_string(N)+"_analytical.txt",x,eigvec);

  // Finding solutions via jacobi rotation method:

  arma::vec eigval;
  int iterations;
  bool converged;

  arma::mat A = create_tridiagonal(N,a,d,a);

  jacobi_eigensolver(A,eps,eigval,eigvec,maxiter,iterations,converged);

  eigvec.shed_cols(3,N-1);
  eigvec.insert_rows(0,1);
  eigvec.insert_rows(N+1,1);

  std::cout << iterations << std::endl;

  write_to_file("solution_N"+std::to_string(N)+".txt",x,eigvec);

}

void write_to_file(std::string filename, arma::vec x, arma::mat y){
  /*
  Prints vector x and y to file with name "filename"
  */

  // Create file:
  std::ofstream outfile ("text_files/"+filename);

  outfile << std::setw(12) << "x";
  for(int i=0; i<size(y)[1]; i++){
    outfile << std::setw(12) << "y"+std::to_string(i)+"(x)";
  }
  outfile << std::endl;

  for(int i=0; i<size(x)[0]; i++){
    outfile << std::setw(12) << std::setprecision(4) << std::scientific << x(i);
    for(int j=0; j<size(y)[1]; j++){
      outfile << std::setw(12) << std::setprecision(4) << std::scientific << y(i,j);
    }
    outfile << std::endl;
  }

  outfile.close();
}

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues,
                        arma::mat& eigenvectors, const int maxiter,
                        int& iterations, bool& converged){

  int N = size(A)[0];

  arma::mat D = A;
  arma::mat R(N, N, arma::fill::eye);

  int k,l;
  double maxval = max_offdiag_symmetric(D,k,l);

  iterations = 0;
  while((abs(maxval) > eps) && (iterations<maxiter)){

    jacobi_rotate(D,R,k,l);

    //std::cout << maxval << std::endl;

    maxval = max_offdiag_symmetric(D,k,l);
    iterations += 1;
  }

  arma::vec eigval(N);
  arma::mat eigvec(N,N);
  arma::uvec sort = sort_index(diagvec(D));

  eigval = diagvec(D);
  eigval = eigval(sort);

  for(int i=0; i<N; i++){
    eigvec.col(i) = arma::normalise(R.col(sort(i)));
  }

  eigenvalues = eigval;
  eigenvectors = eigvec;


  // Checking if convergence criteria was reached before maxiter:
  if(iterations<maxiter){
    converged = true;
  }
  else{
    converged = false;
  }
}

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){

  int N = size(A)[0];

  double tau = (A(l,l)-A(k,k))/(2*A(k,l));

  double t,c,s;

  // Finding tan
  if(tau>=0){
    t = -tau + sqrt(1 + tau*tau);
  }
  else{
    t = -tau - sqrt(1 + tau*tau);
  }

  // Finding cos and sin:
  c = 1/sqrt(1 + t*t);
  s = c*t;

  // Stores updated values so we keep using the old ones for a bit:
  double a_kk,a_ik,r_ik;


  // Updating matrix A:

  a_kk = A(k,k)*c*c - 2*A(k,l)*c*s + A(l,l)*s*s;
  A(l,l) = A(l,l)*c*c + 2*A(k,l)*c*s + A(k,k)*s*s;
  A(k,l) = 0;
  A(l,k) = 0;

  A(k,k) = a_kk;

  for(int i=0; i<N; i++){
    // Updating the i =/= k,l elements of matrix A
    if((i != k) && (i!= l)){
      a_ik = A(i,k)*c - A(i,l)*s;
      A(k,i) = a_ik;
      A(i,l) = A(i,l)*c + A(i,k)*s;
      A(l,i) = A(i,l);

      A(i,k) = a_ik;
    }
  }

  // Updating matrix R:

  for(int i=0; i<N; i++){
    r_ik = R(i,k)*c - R(i,l)*s;
    R(i,l) = R(i,l)*c + R(i,k)*s;

    R(i,k) = r_ik;
  }

}

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l){

  double maxval = 0;
  int N = size(A)[0];

  for(int j=1; j<N; j++){
    for(int i=0; i<j; i++){
      if(abs(A(i,j)) > abs(maxval)){
        maxval = A(i,j);
        k = i;
        l = j;
      }
    }
  }

  return maxval;
}

arma::vec analytical_eigenvalues(int N, double a, double d){

  arma::vec eigval(N);

  double pi = 3.14159265359;

  for(int i=1; i<N+1; i++){
    eigval(i-1) = d + 2*a*cos(i*pi/(N+1));
  }

  return eigval;
}

arma::mat analytical_eigenvectors(int N, double a, double d){

  arma::mat eigvec(N,N);

  double pi = 3.14159265359;

  for(int j=1; j<N+1; j++){
    // Runs through all columns (eigenvectors)
    for(int i=1; i<N+1; i++){
      // Runs through all rows in column (elements in eigenvectors)
      eigvec(i-1,j-1) = sin(i*j*pi/(N+1));
    }
  }

  return arma::normalise(eigvec);
}

arma::mat create_tridiagonal(int N, double a, double d, double e){

  // Start from Diagonal matrix
  arma::mat A = arma::mat(N, N, arma::fill::eye)*d;

  // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
  for(int i = 0; i<N-1; i++){
    A(i+1,i) = a;   // Filling subdiagonal
    A(i,i+1) = e;   // Filling superdiagonal
  }

  // Fill last row (row index n-1)

  return A;
}
