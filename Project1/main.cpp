
#include "functions.hpp"

int main()
{
  // Problem 2:
  // Writing exact solution of u(x) to file

  int n = 10000;            // Number of points in x
  double a = 0,b = 1;       // Start- and end point of x
  double h = (b-a)/(n-1);   // Step size

  std::vector<double> x(n); // Vector to hold x-values
  std::vector<double> u(n); // Vector to hold u-values

  // Setting the values of x and u:
  for (int i = 0; i<n; i++){
    x[i] = a + i*h;
    u[i] = u_func(x[i]);
  }

  load_to_file("u_exact.txt",x,u); // Load to text file



  std::vector<std::vector<double> > xv;
  std::vector<std::vector<double> > x_abs;
  std::vector<std::vector<double> > x_rel;

  for (int i = 1; i<6; i++){
    // Problem 7:
    // Writing approximate solution of u(x) to file using general algorithm,
    // for n = 10^i
    n = (int)pow(10,i);
    xv = approximate_poisson_general(f_func,0,0,0,1,n);
    load_to_file("v_n"+std::to_string(n)+".txt",xv[0],xv[1]);

    // Problem 8a:
    // Writing absolute error of approximate solution to file, for n = 10^i
    x_abs = abs_error(xv[0],xv[1],u_func);
    load_to_file("abs_error_n"+std::to_string(n)+".txt",x_abs[0],x_abs[1]);

    // Problem 8b:
    // Writing relative error of approximate solution to file, for n = 10^i
    x_rel = rel_error(xv[0],xv[1],u_func);
    load_to_file("rel_error_n"+std::to_string(n)+".txt",x_rel[0],x_rel[1]);
  }

  // Problem 9:
  // Writing max error as a function of number of steps n to file,
  // up to n = 10^7
  int n_ = 61;              // Number of points to calculate
  double dn_ = 6.0/(n_-1);  // Step size
  // Vectors to write to file:
  std::vector<double> points(n_),max_error(n_);
  for (int i = 0; i<n_; i++){
    int n = (int)pow(10,1+i*dn_);
    xv = approximate_poisson_general(f_func,0,0,0,1,n);
    x_rel = rel_error(xv[0],xv[1],u_func);

    points[i] = (double)n;
    // Finds element with maximum value:
    max_error[i] = *max_element(x_rel[1].begin(), x_rel[1].end());
  }

  load_to_file("max_error.txt",points,max_error);

  // Problem 10:
  std::vector<std::vector<double> > test;
  std::vector<double> n_points(6),duration_general(6,0),duration_special(6,0);
  clock_t t1,t2;
  for (int i = 1; i<7; i++){
    int n = (int)pow(10,i);
    n_points[i-1] = (double)n;

    for (int j = 0; j<100; j++){
      t1 = clock();
      test = approximate_poisson_general(f_func,0,0,0,1,n);
      t2 = clock();
      duration_general[i-1] += ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    }
    duration_general[i-1] *= 0.01;

    for (int j = 0; j<100; j++){
      t1 = clock();
      test = approximate_poisson_special(f_func,0,0,0,1,n);
      t2 = clock();
      duration_special[i-1] += ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    }
    duration_special[i-1] *= 0.01;
  }

  load_to_file("duration_general.txt",n_points,duration_general);
  load_to_file("duration_special.txt",n_points,duration_special);

  return 0;
}
