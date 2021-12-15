#include "Simulate2DSE.hpp"

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
  sigma_y = 0.1;
  p_x = 200;
  p_y = 0;


  // Finding total probability for box with no barrier:

  Simulate2DSE SE = Simulate2DSE(h, dt, x_c, y_c, sigma_x, sigma_y, p_x, p_y,
                    v0 = 0, slit = "double");
  SE.simulate(T = 0.001);

  total_probability(SE, "no_slit_total_prob.txt");


  // Finding total probability for box with double slit barrier:

  SE = Simulate2DSE(h, dt, x_c, y_c, sigma_x, sigma_y, p_x, p_y, v0 = 1e10, slit = "double");
  SE.simulate(T = 0.001);

  total_probability(SE, "double_slit_total_prob.txt");

  //SE.U_.save("text_files/cube.txt",arma::arma_ascii);

  //U_prob.save("text_files/cube.txt",arma::arma_ascii);



  return 0;
}

void total_probability(Simulate2DSE SE, std::string filename)
{
  arma::cube U_prob = arma::real(arma::conj(SE.U_)%SE.U_);

  int N = arma::size(U_prob)[2];
  arma::vec prob(N);
  for(int n = 0; n<N; n++){
    prob(n) = arma::accu(U_prob.slice(n));
  }

  prob.save("text_files/"+filename,arma::arma_ascii);
}
