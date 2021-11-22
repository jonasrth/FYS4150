#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <armadillo>
#include <omp.h>
#include <random>

class RandomFlipper
{

public:
  // Member variables that hold a <random> generator
  // and distribution (uniform distribution [0,1)).
  mt19937 generator;
  uniform_real_distribution<double> uniform_dist = uniform_real_distribution<double>(0.0 ,1.0);

  // Current state
  arma::mat S_;  // Matrix representing lattice with spins +-1

  double E_;     // Total energy of lattice in coupling constants J
  double M_;     // Total magnetisation of lattice
  double T_;     // Temperature of lattice in J/k_B

  int current_state = 0;


  // Constructor
  RandomWalker(arma::mat S, double T)
  {
    E_ = find_energy(S);
    M_ = arma::accu(S);
    T_ = T;
  }

  // Finds total energy of matrix S
  double find_energy(arma::mat S)
  {

    // Finds total energy by multiplying S elementwise with matrices
    // rolled in the x- and y- direction.

    arma::mat S_roll_x = S, S_roll_y = S;
    int L = arma::size(S)[0];

    S_roll_x.insert_cols(0,S_roll_x.col(L-1));
    S_roll_x.shed_col(L);

    S_roll_y.insert_rows(0,S_roll_y.row(L-1));
    S_roll_y.shed_row(L);

    return -(arma::accu(S%S_roll_x) + arma::accu(S%S_roll_y));
  }

  // Generate random state by flipping one spin
  void update_state()
  {
    double r = uniform_dist(generator);
    if (r < 0.5)
    {
      current_state--;
    }
    else
    {
      current_state++;;
    }
  }
};


double find_energy(arma::mat S);

void initialise(arma::mat& S, int seed);

void update_state(arma::mat& S, int k, int l, double r, int& E, int& M);

int main()
{

  int L = 2;
  int E0 = 0;

  arma::mat S(L,L);

  initialise(S,1583026);

  //std::cout << S << std::endl;
  //std::cout << find_energy(S) << std::endl;

  // Finding energy of lattice:

  int k,l;
  double r;

  // Construct a Mersenne Twister 19937 random number generator with a given seed
  std::mt19937 generator(1583034);

  std::uniform_int_distribution<int> index_distribution(0,L-1);
  std::uniform_int_distribution<int> seed_distribution(0,100000);
  std::uniform_real_distribution<double> real01_distribution(0.0,1.0);


  int N = L*L;
  int n = 1000;

  int E, M;
  double E_sum = 0, M_sum = 0, E2_sum = 0, M2_sum = 0;

  /*
  E = find_energy(S);
  M = arma::accu(S);

  for(int i = 0; i<n; i++){
    for(int j = 0; j<N; j++){
      k = index_distribution(generator);
      l = index_distribution(generator);
      r = real01_distribution(generator);

      update_state(S,k,l,r,E,M);
    }

    E_sum += E;
    E2_sum += E*E;
    M_sum += abs(M);
    M2_sum += M*M;
  }

  std::cout << "  <eps> : " << 1.0*E_sum/n/N << std::endl;
  std::cout << "  <|m|> : " << 1.0*M_sum/n/N << std::endl;
  std::cout << "<eps^2> : " << 1.0*E2_sum/n/N/N << std::endl;
  std::cout << "  <m^2> : " << 1.0*M2_sum/n/N/N << std::endl;
  std::cout << "    c_V : " << (1.0*E2_sum/n - pow(1.0*E_sum/n,2))/(1.0*N) << std::endl;
  std::cout << "    chi : " << (1.0*M2_sum/n - pow(1.0*M_sum/n,2))/(1.0*N) << std::endl;
  */

  L = 20;

  S = arma::mat(L,L,arma::fill::ones);

  initialise(S,12350);

  N = L*L;
  n = 10000;

  E_sum = 0; M_sum = 0; E2_sum = 0; M2_sum = 0;

  E = find_energy(S);
  M = arma::accu(S);

  //std::cout << S << std::endl;
  std::cout << E << std::endl;

  arma::vec E_values(n),M_values(n);

  for(int i = 0; i<n; i++){
    for(int j = 0; j<N; j++){
      k = (int)(L*real01_distribution(generator));
      l = (int)(L*real01_distribution(generator));
      r = real01_distribution(generator);

      update_state(S,k,l,r,E,M);
    }


    E_sum += E;

    /*
    E2_sum += E*E;
    M_sum += abs(M);
    M2_sum += M*M;
    */

    //std::cout << E << std::endl;
    //std::cout << S << std::endl;

    E_values(i) = 1.0*E/N;
    M_values(i) = 1.0*abs(M)/N;
  }

  std::cout << 1.0*E_sum/N/n << std::endl;
  E_values.save("text_files/E_values.txt",arma::arma_ascii);
  M_values.save("text_files/M_values.txt",arma::arma_ascii);
  //std::cout << E_values << std::endl;

  return 0;
}

double find_energy(arma::mat S)
{
  arma::mat S_rx = S, S_ry = S;
  int L = arma::size(S)[0];

  S_rx.insert_cols(0,S_rx.col(L-1));
  S_rx.shed_col(L);

  S_ry.insert_rows(0,S_ry.row(L-1));
  S_ry.shed_row(L);

  return -(arma::accu(S%S_rx) + arma::accu(S%S_ry));
}

void initialise(arma::mat& S, int seed)
{
  std::mt19937 generator(seed);
  std::uniform_int_distribution<int> index_distribution(0,1);

  for(int k = 0; k<arma::size(S)[0]; k++){
    for(int l = 0; l<arma::size(S)[0]; l++){
      S(k,l) = 1 - 2*index_distribution(generator);
    }
  }
}

void update_state(arma::mat& S, int k, int l, double r, int& E, int& M)
{

  int L = arma::size(S)[0];

  int dE = 2*S(k,l)*(S((k+1)%L,l) + S((L+k-1)%L,l) + S(k,(l+1)%L) + S(k,(L+l-1)%L));

  double beta = 1.0;

  if(dE <= 0){
    S(k,l) *= -1;
    E += dE;
    M += 2*S(k,l);
  }
  else if (r < exp(-beta*dE)){
    S(k,l) *= -1;
    E += dE;
    M += 2*S(k,l);
  }
}
