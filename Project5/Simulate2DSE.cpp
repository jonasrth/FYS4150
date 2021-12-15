#include "Simulate2DSE.hpp"

Simulate2DSE::Simulate2DSE(double h, double dt, double x_c, double y_c,
                            double sigma_x, double sigma_y, double p_x, double p_y,
                            double v0, std::string slit)
{
  // Parameters for simulation:
  M_ = (int)(1 + 1/h);
  h_ = h;
  dt_ = dt;

  // Parameters for initial wave packet
  x_c_ = x_c;
  y_c_ = y_c;
  sigma_x_ = sigma_x;
  sigma_y_ = sigma_y;
  p_x_ = p_x;
  p_y_ = p_y;

  // Potential
  v0_ = v0;
  slit_ = slit;
  Simulate2DSE::initialise_v(v_);

  // For solving PDEs
  r_ = arma::cx_double(0,dt/2/h/h);
  Simulate2DSE::make_matrices(A_,B_);

  // Initial states of system
  Simulate2DSE::initialise_u(U_init_);
  Simulate2DSE::mat_to_vec(u_,U_init_);
}

int Simulate2DSE::ij_to_k(int i, int j)
{
  // Takes matrix indexes i,j and return vector index k
  return (i-1) + (j-1)*(M_-2);
}

void Simulate2DSE::initialise_u(arma::cx_mat& u)
{

  //
  // Initialises wave function matrix u as wavepacket defined by class parameters
  // (x_c,y_c,sigma_x,sigma_y,p_x,p_y)
  //

  u = arma::cx_mat(M_,M_,arma::fill::zeros);

  // Saving exp(...) values to use when calculating initial state
  double x, y;
  arma::cx_vec exp_x(M_,arma::fill::zeros),exp_y(M_,arma::fill::zeros);
  for(int k=1; k<M_-1; k++){
    x = k*h_;
    y = k*h_;
    exp_x(k) = exp(arma::cx_double(-pow((x - x_c_)/sigma_x_,2)/2, p_x_*(x - x_c_)));
    exp_y(k) = exp(arma::cx_double(-pow((y - y_c_)/sigma_y_,2)/2, p_y_*(y - y_c_)));
  }

  // Using saved exp(...) values to calculate initial state
  for(int i=1; i<M_-1; i++){
    for(int j=1; j<M_-1; j++){
      u(i,j) = exp_x(i)*exp_y(j);
    }
  }

  // Normalising initial state:
  double sum = arma::accu(arma::real(arma::conj(u)%u));
  u = u/sqrt(sum);
}

void Simulate2DSE::initialise_v(arma::mat& v)
{

  //
  // Initialises potential matrix v with potential shape and magnitude defined
  // by class parameters (v0, slit)
  //

  v = arma::mat(M_,M_,arma::fill::zeros);

  for(int i=1; i<M_-1; i++){

    // Setting up barrier
    if((i*h_ >= 0.49) && (i*h_ <= 0.51)){
      v.row(i) = arma::mat(1,M_,arma::fill::value(v0_));

      // Setting up single slit
      if(slit_ == "single"){
        for(int j=1; j<M_-1; j++){

          // Creating slit
          if((j*h_ >= 0.475) && (j*h_ <= 0.525)){
            v(i,j) = 0;
          }
        }

      // Setting up double slit
      }else if(slit_ == "double"){
        for(int j=1; j<M_-1; j++){

          // Creating slits
          if((j*h_ >= 0.425) && (j*h_ <= 0.475)){
            v(i,j) = 0;
          }else if((j*h_ >= 0.525) && (j*h_ <= 0.575)){
            v(i,j) = 0;
          }
        }

      // Setting of triple slit
      }else if(slit_ == "triple"){
        for(int j=1; j<M_-1; j++){

          // Creating slits
          if((j*h_ >= 0.375) && (j*h_ <= 0.425)){
            v(i,j) = 0;
          }else if((j*h_ >= 0.475) && (j*h_ <= 0.525)){
            v(i,j) = 0;
          }else if((j*h_ >= 0.575) && (j*h_ <= 0.625)){
            v(i,j) = 0;
          }
        }
      }
    }
  }
}

void Simulate2DSE::vec_to_mat(arma::cx_mat& A, arma::cx_vec a)
{

  //
  // Turns vector representing internal points of the box into matrix
  //

  A = arma::cx_mat(M_,M_,arma::fill::zeros);

  int k;
  for(int i=1; i<M_-1; i++){
    for(int j=1; j<M_-1; j++){
      k = ij_to_k(i,j);
      A(i,j) = a(k);
    }
  }
}

void Simulate2DSE::mat_to_vec(arma::cx_vec& a, arma::cx_mat A)
{

  //
  // Turns matrix representing the box into vector containing the internal points
  //

  int N = (M_-2)*(M_-2);
  a = arma::cx_vec(N);

  int k;
  for(int i=1; i<M_-1; i++){
    for(int j=1; j<M_-1; j++){
      k = ij_to_k(i,j);
      a(k) = A(i,j);
    }
  }
}

void Simulate2DSE::fill_matrix(arma::sp_cx_mat& A, arma::cx_vec a, arma::cx_double r)
{

  //
  // Fills matrix A in the shape defined by the Crank-Nicolson discretisation
  // of the Schrödinger equation. The innter diagonal is set as the vector a,
  // while the outer diagonals are set at the double r.
  //

  int m = M_-2;
  int N = m*m;
  A = arma::sp_cx_mat(N,N);

  // Filling diagonal
  A.diag() = a;

  // Filling outer r diagonals
  for(int i=0; i<N-m; i++)
  {
    A(i,i+m) = r;
    A(i+m,i) = r;
  }

  // Filling inner r diagonals
  for(int i=0; i<N-1; i++)
  {
    A(i,i+1) = r;
    A(i+1,i) = r;
  }

  // Removing appropriate elements in inner r diagonal
  for(int i=0; i<m-1; i++)
  {
    A(m*(i+1)-1,m*(i+1)) = 0;
    A(m*(i+1),m*(i+1)-1) = 0;
  }
}

void Simulate2DSE::make_matrices(arma::sp_cx_mat& A, arma::sp_cx_mat& B)
{

  //
  // Makes matrices A and B to be used in the matrix equation representing the
  // Schrödinger equation.
  //

  int N = (M_-2)*(M_-2);
  arma::cx_vec a(N), b(N);

  int k;
  arma::cx_double c;

  // Setting up diagonals a and b for matrices A and B
  for(int i=1; i<M_-1; i++){
    for(int j=1; j<M_-1; j++){
      k = ij_to_k(i,j);
      c = 4.0*r_ + arma::cx_double(0,dt_/2*v_(i,j));
      a(k) = 1.0 + c;
      b(k) = 1.0 - c;
    }
  }

  // Setting up matrices A and B
  fill_matrix(A,a,-r_);
  fill_matrix(B,b,+r_);
}

void Simulate2DSE::simulate(double T)
{

  //
  // Solves the next time step after a time dt_, until reaching a time T.
  //

  // Number of time steps
  int N = (int)(T/dt_ + 1);

  // Setting up complex cube holding all states
  U_ = arma::cx_cube(M_,M_,N);
  vec_to_mat(U_.slice(0),u_);

  arma::cx_vec b;

  // Solving as symmetric matrix equation since A is symmetric
  arma::superlu_opts opts;
  opts.symmetric = true;

  // Performing simulation:
  for(int n=1; n<N; n++){
    b = B_*u_;

    // Solving matrix equation A_*u_ = b, using LU decomposition
    u_ = spsolve(A_, b, "superlu");

    // Saving current state
    vec_to_mat(U_.slice(n),u_);
  }
}
