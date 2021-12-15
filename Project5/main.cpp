#include "Simulate2DSE.hpp"

// Finds total probability of each slice of wave function cube (each state)
void total_probability(Simulate2DSE SE, std::string filename);

int main()
{

  // Parameters:

  double h,dt,x_c,y_c,sigma_x,sigma_y,p_x,p_y,v0,T;
  std::string slit;

  h = 0.005;
  dt = 2.5e-5;
  x_c = 0.25;
  y_c = 0.5;
  sigma_x = 0.05;
  sigma_y = 0.05;
  p_x = 200;
  p_y = 0;

  //
  // Problem 7
  //

  T = 0.008;

  // Finding total probability over time for box with no barrier:

  v0 = 0;

  Simulate2DSE SE = Simulate2DSE(h, dt, x_c, y_c, sigma_x, sigma_y, p_x, p_y, v0);
  SE.simulate(T);
  total_probability(SE, "no_slit_total_prob.txt");


  // Finding total probability over time for box with double slit barrier:

  v0 = 1e10;
  sigma_y = 0.1;
  slit = "double";

  SE = Simulate2DSE(h, dt, x_c, y_c, sigma_x, sigma_y, p_x, p_y, v0, slit);
  SE.simulate(T);
  total_probability(SE, "double_slit_total_prob.txt");



  //
  // Problem 8
  //

  T = 0.002;

  v0 = 1e10;
  sigma_y = 0.2;
  slit = "double";

  SE = Simulate2DSE(h, dt, x_c, y_c, sigma_x, sigma_y, p_x, p_y, v0, slit);
  SE.simulate(T);

  // Cube to hold snapshots:
  arma::cx_cube snapshots(SE.M_,SE.M_,3);
  // x and y dimension of U_ (same as M-2):
  int N = arma::size(SE.U_)[2];

  // Saving snapshots of t = 0, t = 0.001, t = 0.002:
  snapshots.slice(0) = SE.U_.slice(0);
  snapshots.slice(1) = SE.U_.slice((int)((N-1)/2));
  snapshots.slice(2) = SE.U_.slice(N-1);

  snapshots.save("text_files/double_slit_snapshots.txt", arma::arma_ascii);



  //
  // Problem 9
  //


  // Saving wave function at t=0.002 for single slit barrier potential
  slit = "single";
  SE = Simulate2DSE(h, dt, x_c, y_c, sigma_x, sigma_y, p_x, p_y, v0, slit);
  SE.simulate(T);
  SE.U_.slice(N-1).save("text_files/p_single.txt",arma::arma_ascii);


  // Saving wave function at t=0.002 for double slit barrier potential
  slit = "double";
  SE = Simulate2DSE(h, dt, x_c, y_c, sigma_x, sigma_y, p_x, p_y, v0, slit);
  SE.simulate(T);
  SE.U_.slice(N-1).save("text_files/p_double.txt",arma::arma_ascii);

  // Saving wave function at t=0.002 for triple slit barrier potential
  slit = "triple";
  SE = Simulate2DSE(h, dt, x_c, y_c, sigma_x, sigma_y, p_x, p_y, v0, slit);
  SE.simulate(T);
  SE.U_.slice(N-1).save("text_files/p_triple.txt",arma::arma_ascii);



  return 0;
}

void total_probability(Simulate2DSE SE, std::string filename)
{

  //
  // Finds total probability of each slice of wave function cube (each state).
  // Writes results to file with name filename
  //

  arma::cube U_prob = arma::real(arma::conj(SE.U_)%SE.U_);

  int N = arma::size(U_prob)[2];
  arma::vec prob(N);
  for(int n = 0; n<N; n++){
    prob(n) = arma::accu(U_prob.slice(n));
  }

  prob.save("text_files/"+filename,arma::arma_ascii);
}
