#include "PenningTrap.hpp"
#include <omp.h>

// Writes to matrix the solution for the movement of particles in penning trap PT
// for N steps over time T. Can turn on/off interactions between particles
// with bool interaction. particles variable denotes how many particles' motion
// to save.
arma::mat simulate(PenningTrap PT, std::string algo, int particles, double T, int N, bool interaction);

// Calculates analytical solution for single particle in penning trap, with
// initial particle state on form r = {x0,0,z0} and v = {0,v0,0}.
arma::mat analytical_solution(PenningTrap PT, double T, int N);

arma::mat relative_error(PenningTrap PT, std::string algo, double T, double h);

double error_convergence_rate(PenningTrap PT, std::string algo, double T, arma::vec h);

double count_particles_outside(PenningTrap PT);

arma::mat particles_outside(PenningTrap PT, double T, int N, arma::vec omegaV, bool interaction);

int main()
{

  // Initialising PenningTrap
  double B0 = 1, V0 = 10, d = 1;
  PenningTrap PT(B0,V0,d);

  // Initialising and adding particles to PenningTrap
  double q = 1, m = 40.078;
  Particle p(q,m,{1,0,1},{0,1,0});
  PT.add_particle(p);
  p = Particle(q,m,{-1,0,-1},{0,-1,0});
  PT.add_particle(p);

  // Setting simulation time T and number of steps N
  double T = 100;
  int N = 10000;

  // Matrix for holding simulation results:
  arma::mat S;

  /*
  // Simulating two particles without interactions:
  S = simulate(PT,"RK4",2,T,N,false);
  S.save("text_files/2p_positions.txt",arma::arma_ascii);

  // Finding analytical solution:
  S = analytical_solution(PT,T,N);
  S.save("text_files/analytical_positions.txt",arma::arma_ascii);

  // Simulating two particles with interactions:
  S = simulate(PT,"RK4",2,T,N,true);
  S.save("text_files/2p_interaction_positions.txt",arma::arma_ascii);


  PT = PenningTrap(B0,V0,d);
  p = Particle(q,m,{1,0,1},{0,1,0});
  PT.add_particle(p);

  // Step sizes to test error for:
  arma::vec h = arma::logspace(0,-4,5);

  // Finding relative error for different h, using Forward Euler:
  for(int i = 0; i<arma::size(h)[0]; i++){
    arma::mat rel_error = relative_error(PT,"FE",T,h[i]);
    rel_error.save("text_files/FE_rel_error_h"+std::to_string(i)+".txt",arma::arma_ascii);
  }

  // Finding relative error for different h, using Runge-Kutta 4:
  for(int i = 0; i<arma::size(h)[0]; i++){
    arma::mat rel_error = relative_error(PT,"RK4",T,h[i]);
    rel_error.save("text_files/RK4_rel_error_h"+std::to_string(i)+".txt",arma::arma_ascii);
  }

  // Defining error convergence variable
  double r_err;

  // Error convergence for Forward Euler:
  r_err = error_convergence_rate(PT,"FE",T,h);
  std::cout << "r_err for FE: " << r_err << std::endl;

  // Error convergence for Runge-Kutta 4:
  r_err = error_convergence_rate(PT,"RK4",T,h);
  std::cout << "r_err for RK4: " << r_err << std::endl;

  */

  // Initialising new PenningTrap

  V0 = 0.0025;
  d = 0.05;
  PT = PenningTrap(B0,V0,d,0.7);

  int N_p = 100;  // Number of particles
  arma::vec r,v;  // Initialising random position and velocity vectors

  arma::arma_rng::set_seed(12345);  // Setting a seed for the RNG

  for(int i = 0; i < N_p; i++){
    r = arma::vec(3).randn() * 0.1 * PT.d_;   // Random initial position
    v = arma::vec(3).randn() * 0.1 * PT.d_;   // Random initial velocity

    // Add particle with random initial position and velocity
    p = Particle(q,m,r,v);
    PT.add_particle(p);
  }

  arma::vec omegaV_values = arma::linspace(0.2,2.5,150);

  // Finding number of escaped particles as function of induced frequency

  /*
  PT.f_ = 0.1;
  S = particles_outside(PT,500,10000,omegaV_values,false);
  S.save("text_files/NW_f0.1.txt",arma::arma_ascii);

  PT.f_ = 0.4;
  S = particles_outside(PT,500,10000,omegaV_values,false);
  S.save("text_files/NW_f0.4.txt",arma::arma_ascii);

  PT.f_ = 0.7;
  S = particles_outside(PT,500,10000,omegaV_values,false);
  S.save("text_files/NW_f0.7.txt",arma::arma_ascii);
  */

  // Finds particles ejected for induced frequency range (0.42,0.46)
  // with interactions off:
  PT.f_ = 0.1;
  omegaV_values = arma::linspace(0.42,0.46,50);
  S = particles_outside(PT,500,10000,omegaV_values,false);
  S.save("text_files/NW_f0.1_zoom.txt",arma::arma_ascii);

  // Finds particles ejected for induced frequency range (0.42,0.46)
  // with interaction on:
  PT.f_ = 0.1;
  omegaV_values = arma::linspace(0.42,0.46,50);
  S = particles_outside(PT,500,10000,omegaV_values,true);
  S.save("text_files/NW_f0.1_zoom_interaction.txt",arma::arma_ascii);

  return 0;
}

arma::mat particles_outside(PenningTrap PT, double T, int N, arma::vec omegaV, bool interaction)
{

  double dt = T/(N-1);
  arma::mat NW(arma::size(omegaV)[0],2);

  #pragma omp parallel for
  for(int k = 0; k < arma::size(omegaV)[0]; k++){
    PenningTrap PT_ = PT;
    PT_.omegaV_ = omegaV[k];
    for(int i = 0; i < N-1; i++){
      PT_.evolve_RK4(dt, interaction);
    }
    NW.row(k)[0] = omegaV[k];
    NW.row(k)[1] = count_particles_outside(PT_)/PT_.Particles_.size();
    std::cout << std::setw(12) << NW.row(k)[0] << std::setw(12) <<  NW.row(k)[1] << std::endl;
  }
  return NW;
}

double count_particles_outside(PenningTrap PT)
{
  double n = 0;
  for(int i = 0; i < PT.Particles_.size(); i++){
    if(arma::norm(PT.Particles_[i].r_) > PT.d_){
      n += 1;
    }
  }
  return n;
}

double error_convergence_rate(PenningTrap PT, std::string algo, double T, arma::vec h)
{
  arma::vec max_error = h;

  double abs_error;

  // Finding max error for step sizes in h:
  for(int k = 0; k < arma::size(h)[0]; k++){

    int N = (int)T/h[k] + 1;

    arma::mat r_exact = analytical_solution(PT,T,N).cols(1,3);
    arma::mat r_approx = simulate(PT,algo,1,T,N,false).cols(1,3);

    max_error[k] = 0;

    for(int i = 0; i < N; i++){
      abs_error = abs(arma::norm(r_exact.row(i) - r_approx.row(i)));
      if(abs_error > max_error[k]){
        max_error[k] = abs_error;
      }
    }
  }

  // Finding error convergence rate:
  double r_err = 0;
  for(int k = 1; k < arma::size(h)[0]; k++){
    r_err += log(max_error[k]/max_error[k-1])/log(h[k]/h[k-1]);
  }

  return 0.25*r_err;
}

arma::mat relative_error(PenningTrap PT, std::string algo, double T, double h)
{
  int N = (int)T/h + 1;
  arma::mat rel_error(N,2);

  rel_error.col(0) = arma::linspace(0,T,N);

  arma::mat r_exact = analytical_solution(PT,T,N).cols(1,3);
  arma::mat r_approx = simulate(PT,algo,1,T,N,false).cols(1,3);

  for(int i = 0; i < N; i++){
    rel_error.row(i)[1] = abs(arma::norm(r_exact.row(i) - r_approx.row(i))
                          /arma::norm(r_exact.row(i)));
  }

  return rel_error;
}

arma::mat simulate(PenningTrap PT, std::string algo, int particles, double T, int N, bool interaction)
{

  double dt = T/(N-1);

  double t = 0;
  arma::vec r,v;
  arma::mat trv(N,1 + 6*particles);

  if(algo == "RK4"){
    for(int i = 0; i < N; i++){
        trv.row(i)[0] = t;
        for(int p = 0; p < particles; p++){
          r = PT.Particles_[p].r_;
          v = PT.Particles_[p].v_;
          for(int j = 0; j < 3; j++){
            trv.row(i)[1 + 6*p + j] = r[j];
            trv.row(i)[4 + 6*p + j] = v[j];
          }
        }

        t = t + dt;
        PT.evolve_RK4(dt,interaction);
    }
  }
  else{
    for(int i = 0; i < N; i++){
        trv.row(i)[0] = t;
        for(int p = 0; p < particles; p++){
          r = PT.Particles_[p].r_;
          v = PT.Particles_[p].v_;
          for(int j = 0; j < 3; j++){
            trv.row(i)[1 + 6*p + j] = r[j];
            trv.row(i)[4 + 6*p + j] = v[j];
          }
        }

        t = t + dt;
        PT.evolve_forward_Euler(dt,interaction);
      }
  }

  return trv;
}

arma::mat analytical_solution(PenningTrap PT, double T, int N)
{
  double dt = T/N;

  Particle p = PT.Particles_[0];

  double x0 = p.r_[0], z0 = p.r_[2], v0 = p.v_[1];

  double omega_0,omega_z,omega_p,omega_m;

  omega_0 = p.q_*PT.B0_/p.m_;
  omega_z = sqrt(2*p.q_*PT.V0_/(PT.d_*PT.d_)/p.m_);

  omega_p = (omega_0 + sqrt(omega_0*omega_0 - 2*omega_z*omega_z))/2;
  omega_m = (omega_0 - sqrt(omega_0*omega_0 - 2*omega_z*omega_z))/2;

  double A_p, A_m;

  A_p = (v0 + omega_m*x0)/(omega_m - omega_p);
  A_m = -(v0 + omega_p*x0)/(omega_m - omega_p);

  double x,y,z;
  arma::mat sol(N,4);
  arma::vec t = arma::linspace(0,T,N);

  for(int i = 0; i < N; i++){

    x = A_p*cos(omega_p*t[i]) + A_m*cos(omega_m*t[i]);
    y = -A_p*sin(omega_p*t[i]) - A_m*sin(omega_m*t[i]);
    z = z0*cos(omega_z*t[i]);

    sol.row(i) = {t[i],x,y,z};
  }

  return sol;
}
