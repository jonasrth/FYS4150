#include "RandomFlipper.hpp"

// initialises S randomly
void initialise(arma::mat& S, int seed);

// loops over n MC cycles at a temperature T for a lattic of size L
// unordered - bool to tell is initialisation of lattice should choose
//             random spin direction or point them all in the same direction
// save - bool to tell if you want to save to file
// print - bool to tell if you want to print expectation values, etc.
void MC_cycles(std::string filename, int L, double T, int n,
                bool unordered = true, bool save = true, bool print = false);

// loops over temperatures and saves quantities to file "filename"
// L - size of lattice
// Tmin - min temperature value
// Tmax - max temperature value
// T_steps - numper of steps to take
// n - number of monte carlo cycles
// n - burn-in time to use
void temperature_loop(std::string filename, int L, double Tmin, double Tmax,
                      int T_steps, int n, int n_burnin);

int main()
{

  // Printing quantities for 2x2 lattice at T = 1.0 for 10^6 MC cycles:
  std::cout << "For L = 2, T = 1.0 J/k_B and n = 10^6:" << std::endl;
  MC_cycles("", 2, 1.0, 1000000, true, false, true);
  std::cout << std::endl;

  // Saving 20x20 lattice
  MC_cycles("EM_T1_unordered_L20", 20, 1.0, 10000);
  MC_cycles("EM_T1_ordered_L20", 20, 1.0, 10000, false);

  MC_cycles("EM_T2_unordered_L20", 20, 2.4, 10000);
  MC_cycles("EM_T2_ordered_L20", 20, 2.4, 10000, false);

  std::cout << "For L = 20, T = 1.0 J/k_B and n = 10^6:" << std::endl;
  MC_cycles("EM_T1_L20", 20, 1.0, 1000000, true, true, true);
  std::cout << std::endl;
  std::cout << "For L = 20, T = 2.4 J/k_B and n = 10^6:" << std::endl;
  MC_cycles("EM_T2_L20", 20, 2.4, 1000000, true, true, true);
  std::cout << std::endl;


  std::vector<int> L_values = {40,60,80,100};
  std::string filename;

  // Finding quantities as function of temperature
  // WARNING: takes hours to run.

  for(int i = 0; i < L_values.size(); i++){
    filename = "init_temps_L"+std::to_string(L_values[i]);
    std::cout << filename << std::endl;
    temperature_loop(filename, L_values[i], 2.1, 2.4, 10, 1000000, 2000);
  }

  for(int i = 0; i < L_values.size(); i++){
    filename = "zoom_peak_temps_2_L"+std::to_string(L_values[i]);
    std::cout << filename << std::endl;
    temperature_loop(filename, L_values[i], 2.26, 2.30, 40, 1000000, 2000);
  }

  return 0;
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

void MC_cycles(std::string filename, int L, double T, int n, bool unordered,
                bool save, bool print)
{

  int N = L*L;
  arma::mat S(L,L);

  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();

  if (unordered == true){
    initialise(S,seed);
  }
  else{
    S.ones();
  }

  RandomFlipper RF(S,T);
  RF.generator.seed(seed+1);

  // To find expectation values:
  double E_sum = 0, E2_sum = 0, M_sum = 0, M2_sum = 0;

  // Matrix to hold mean energy per spin and mean magnetisation per spin
  // at the state at the end of each MC cycle
  arma::mat EM(n,2);

  for(int i = 0; i<n; i++){

    RF.monte_carlo_cycle();

    E_sum += RF.E_;
    E2_sum += RF.E_*RF.E_;
    M_sum += abs(RF.M_);
    M2_sum += RF.M_*RF.M_;

    //std::cout << E_sum << std::endl;

    EM(i,0) = 1.0*RF.E_/N;
    EM(i,1) = 1.0*abs(RF.M_)/N;
  }

  if (save == true){
    EM.save("text_files/"+filename+".txt",arma::arma_ascii);
  }
  if (print == true){
    std::cout << "<eps>   [J]   : " << (double)E_sum/N/n << std::endl;
    std::cout << "<|m|>   [1]   : " << (double)M_sum/N/n << std::endl;
    std::cout << "<eps^2> [J^2] : " << (double)E2_sum/N/N/n << std::endl;
    std::cout << "<m^2>   [1]   : " << (double)M2_sum/N/N/n << std::endl;
    std::cout << "C_V     [k_B] : " << (double)((E2_sum/n - E_sum*E_sum/n/n)/N/T/T) << std::endl;
    std::cout << "chi     [1/J] : " << (double)((M2_sum/n - M_sum*M_sum/n/n)/N/T) << std::endl;
  }
}

void temperature_loop(std::string filename, int L, double Tmin, double Tmax,
                      int T_steps, int n, int n_burnin)
{
  int N = L*L;

  arma::vec T = arma::linspace(Tmin,Tmax,T_steps+1);
  arma::mat data(T_steps+1,5);

  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();

  int n_eff = n - n_burnin;

  #pragma omp parallel
  {
    const int my_thread = omp_get_thread_num();

    #pragma omp for
    for(int i = 0; i<(T_steps+1); i++){

      arma::mat S(L,L);
      initialise(S,seed + my_thread + 1);

      RandomFlipper RF(S,Tmin);
      RF = RandomFlipper(S,T(i));
      RF.generator.seed(seed - my_thread - 1);


      // Not saving burn-in values
      for(int j = 0; j<n_burnin; j++){
        RF.monte_carlo_cycle();
      }

      // Saving values after burn-in time n_burnin
      for(int j = n_burnin; j<n; j++){

        RF.monte_carlo_cycle();

        data(i,1) += RF.E_;
        data(i,2) += abs(RF.M_);
        data(i,3) += RF.E_*RF.E_;
        data(i,4) += RF.M_*RF.M_;
      }

      data(i,0) = T(i);
      data(i,3) = (data(i,3)/n_eff - pow(data(i,1)/n_eff,2))/T(i)/T(i)/N;
      data(i,4) = (data(i,4)/n_eff - pow(data(i,2)/n_eff,2))/T(i)/N;
      data(i,1) = data(i,1)/N/n_eff;
      data(i,2) = data(i,2)/N/n_eff;

    }
  }
  data.save("text_files/"+filename+".txt",arma::arma_ascii);
}
