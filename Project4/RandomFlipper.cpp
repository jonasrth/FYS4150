#include "RandomFlipper.hpp"

RandomFlipper::RandomFlipper(arma::mat S, double T)
{
  S_ = S;
  L_ = arma::size(S)[0];

  E_ = RandomFlipper::find_energy();
  M_ = arma::accu(S);
  T_ = T;
  w_ = arma::vec(17);
  w_(0) = exp(8.0/T); w_(4) = exp(4.0/T); w_(8) = 1;
  w_(12) = exp(-4.0/T); w_(16) = exp(-8.0/T);

  //w_ = {exp(8.0/T), exp(4.0/T), 1.0, exp(-4.0/T), exp(-8.0/T)};
}

double RandomFlipper::find_energy()
{

  // Finds total energy by multiplying S elementwise with matrices
  // rolled in the x- and y- direction.

  arma::mat S_roll_x = S_, S_roll_y = S_;

  S_roll_x.insert_cols(0,S_roll_x.col(L_-1));
  S_roll_x.shed_col(L_);

  S_roll_y.insert_rows(0,S_roll_y.row(L_-1));
  S_roll_y.shed_row(L_);

  return -(arma::accu(S_%S_roll_x) + arma::accu(S_%S_roll_y));
}

void RandomFlipper::update_state()
{

  int k = (int)(L_*uniform_dist(generator));
  int l = (int)(L_*uniform_dist(generator));
  double r = uniform_dist(generator);

  int dE = 2*S_(k,l)*(S_((k+1)%L_,l) + S_((L_+k-1)%L_,l)
                      + S_(k,(l+1)%L_) + S_(k,(L_+l-1)%L_));

  // Transform dE into index to compare Boltzmann factor to r
  if(r < w_(dE + 8)){
    S_(k,l) *= -1;   // Update spin state
    E_ += dE;        // Update energy
    M_ += 2*S_(k,l);  // Update magnetisation
  }
}

void RandomFlipper::monte_carlo_cycle()
{
  int N = L_*L_;
  for(int i = 0; i<N; i++){
    RandomFlipper::update_state();
  }
}
