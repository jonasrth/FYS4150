#include "PenningTrap.hpp"

// Definition of the constructor
PenningTrap::PenningTrap(double B0, double V0, double d, double f, double omegaV)
{
  // Assign input values (B0,V0,d) to the class member variables (B0_,V0_,d_,)
  B0_ = B0*T;
  V0_ = V0*V;
  d_ = d*cm;
  f_ = f;
  omegaV_ = omegaV;
}

void PenningTrap::add_particle(Particle p)
{
  Particles_.push_back(p);
}

arma::vec PenningTrap::E_field(arma::vec r)
{
  arma::vec scale = {1,1,-2};
  return V0_*(1 + f_*cos(omegaV_*t_))/(d_*d_)*(r%scale);
}

arma::vec PenningTrap::B_field(arma::vec r)
{
  arma::vec e_z = {0,0,1};
  return B0_*e_z;
}

arma::vec PenningTrap::F_external(int i)
{
  Particle p = Particles_[i];

  arma::vec r = p.r_, v = p.v_;
  double q = p.q_;

  arma::vec E = PenningTrap::E_field(r), B = PenningTrap::B_field(r);

  arma::vec F;
  if(arma::norm(r) > d_){
    // E- and B-field set to zero when outside penningtrap (|r| > d)
    F = {0,0,0};
  }
  else{
    F = q*(E + arma::cross(v,B));
  }

  return F;
}

arma::vec PenningTrap::F_particle(int i, int j)
{
  Particle p_i = Particles_[i], p_j = Particles_[j];

  double q_i = p_i.q_, q_j = p_j.q_;
  arma::vec r_i = p_i.r_, r_j = p_j.r_;

  return k_e*q_i*q_j*(r_i - r_j)/pow(arma::norm(r_i - r_j),3);
}

arma::vec PenningTrap::F_particle_tot(int i)
{
  arma::vec F = {0,0,0};
  for(int j = 0; j < Particles_.size(); j++){
    if(i != j){
      F += F_particle(i,j);
    }
  }
  return F;
}

arma::vec PenningTrap::F_tot(int i, bool interaction)
{
  arma::vec F;
  if(interaction == true){
    F = F_particle_tot(i) + F_external(i);
  }
  else{
    F = F_external(i);
  }
  return F;
}

void PenningTrap::evolve_forward_Euler(double dt, bool interaction)
{
  t_ += dt;
  for(int i = 0; i < Particles_.size(); i++){
    Particle p_i = Particles_[i];
    arma::vec a_i = F_tot(i,interaction)/p_i.m_;
    Particles_[i].v_ = p_i.v_ + a_i*dt;
    Particles_[i].r_ = p_i.r_ + p_i.v_*dt;

  }
}

void PenningTrap::evolve_RK4(double dt, bool interaction)
{
  double h = dt;
  Particle p_i = Particles_[0];
  int N = Particles_.size();

  arma::cube u(3,2,N), u_new(3,2,N);
  arma::cube k1(3,2,N),k2(3,2,N),k3(3,2,N),k4(3,2,N);

  for(int i = 0; i < N; i++){
    p_i = Particles_[i];

    // Saving state at time t:
    u.slice(i).col(0) = p_i.r_;
    u.slice(i).col(1) = p_i.v_;

    // Finding k1:
    k1.slice(i).col(0) = h*p_i.v_;
    k1.slice(i).col(1) = h*F_tot(i,interaction)/p_i.m_;
  }

  // k1 and k2 found at t + dt/2
  t_ += dt/2;

  for(int i = 0; i < N; i++){
    // Guess for state at t + dt/2:
    Particles_[i].r_ = u.slice(i).col(0) + k1.slice(i).col(0)/2;
    Particles_[i].v_ = u.slice(i).col(1) + k1.slice(i).col(1)/2;
  }

  for(int i = 0; i < N; i++){
    p_i = Particles_[i];

    // Finding k2:
    k2.slice(i).col(0) = h*p_i.v_;
    k2.slice(i).col(1) = h*F_tot(i,interaction)/p_i.m_;

    // Other guess for state at t + dt/2
    Particles_[i].r_ = u.slice(i).col(0) + k2.slice(i).col(0)/2;
    Particles_[i].v_ = u.slice(i).col(1) + k2.slice(i).col(1)/2;
  }

  for(int i = 0; i < N; i++){
    p_i = Particles_[i];

    // Finding k3:
    k3.slice(i).col(0) = h*p_i.v_;
    k3.slice(i).col(1) = h*F_tot(i,interaction)/p_i.m_;
  }

  // k4 found at t + dt
  t_ += dt/2;

  for(int i = 0; i < N; i++){
    // Guess for state at t + dt:
    Particles_[i].r_ = u.slice(i).col(0) + k3.slice(i).col(0);
    Particles_[i].v_ = u.slice(i).col(1) + k3.slice(i).col(1);
  }

  for(int i = 0; i < N; i++){
    p_i = Particles_[i];

    // Finding k4:
    k4.slice(i).col(0) = h*p_i.v_;
    k4.slice(i).col(1) = h*F_tot(i,interaction)/p_i.m_;
  }

  // New positions and velocities:
  u_new = u + (1./6)*(k1 + 2*k2 + 2*k3 + k4);

  for(int i = 0; i < N; i++){
    Particles_[i].r_ = u_new.slice(i).col(0);
    Particles_[i].v_ = u_new.slice(i).col(1);
  }


}
