//
// Definitions of functions declared in "functions.hpp"
//

#include "functions.hpp"


double u_func(double x){
  // Returns u(x) for a given x
  return 1 - (1-exp(-10))*x - exp(-10*x);
}

double f_func(double x){
  // Returns f(x) for a given x
  return 100*exp(-10*x);
}

std::vector<std::vector<double> > approximate_poisson_general(double (func)(double),
                                  double u_init, double u_end,
                                  double x_init, double x_end, int steps){
  /*
  General algorithm for solving a matrix equation with a general tridiagonal
  matrix.

  Parameters:
        func - RHS of Poisson equation
        u_init - u(x_init)
        u_end - u(x_end)
        x_init - initial value of x
        x_end - final value of x

  Returns 2D vector containing x and approximation to u, v(x)
  */
  int N = steps;
  int n = N-2;
  double h = (x_end - x_init)/(N - 1);

  std::vector<double> a(n,-1),b(n,2),c(n,-1),g(n);
  std::vector<double> x(N),v(N);

  // Setting all values of x
  x[0] = x_init;
  x[N-1] = x_end;
  for(int i = 1; i<N-1; i++){
    x[i] = x_init + i*h;
  }

  // Setting endpoints of v:
  v[0] = u_init;
  v[N-1] = u_end;

  // Setting all values of solution vector g
  g[0] = func(x[1])*h*h - a[0]*u_init;
  g[n-1] = func(x[n])*h*h - c[0]*u_end;
  for(int i = 1; i<n-1; i++){
    g[i] = func(x[i+1])*h*h;
  }

  // Forward substition:
  for(int i = 1; i<n; i++){
    b[i] = b[i] - (a[i]/b[i-1])*c[i-1];
    g[i] = g[i] - (a[i]/b[i-1])*g[i-1];
  }

  // Backward substition:
  v[n] = g[n-1]/b[n-1];
  for(int i = n-2; i>-1; i--){
    v[i+1] = (g[i] - c[i]*v[i+2])/b[i];
  }

  std::vector<std::vector<double> > xv(2);
  xv[0] = x;
  xv[1] = v;

  return xv;
}

std::vector<std::vector<double> > approximate_poisson_special(double (func)(double),
                                  double u_init, double u_end,
                                  double x_init, double x_end, int steps){
  /*
  Special algorithm for solving the tridiagonal matrix equation with
  matrix signature (-1,2,-1).

  Parameters:
        func - RHS of Poisson equation
        u_init - u(x_init)
        u_end - u(x_end)
        x_init - initial value of x
        x_end - final value of x

  Returns 2D vector containing x and approximation to u, v(x)
  */

  int N = steps;
  int n = N-2;
  double h = (x_end - x_init)/(N - 1);

  std::vector<double> a(n,-1),b(n,2),c(n,-1),g(n);
  std::vector<double> x(N),v(N);

  // Setting all values of x:
  x[0] = x_init;
  x[N-1] = x_end;
  for(int i = 1; i<N-1; i++){
    x[i] = x_init + i*h;
  }

  // Setting endpoints of v:
  v[0] = u_init;
  v[N-1] = u_end;

  // Setting all values of solution vector g:
  g[0] = func(x[1])*h*h + u_init;
  g[n-1] = func(x[n])*h*h + u_end;
  for(int i = 1; i<n-1; i++){
    g[i] = func(x[i+1])*h*h;
  }

  // Forward substitiution:
  for(int i = 1; i<n; i++){
    b[i] = b[i] - 1/b[i-1];
    g[i] = g[i] + g[i-1]/b[i-1];
  }

  // Backward substitution:
  v[n] = g[n-1]/b[n-1];
  for(int i = n-2; i>-1; i--){
    v[i+1] = (g[i] + v[i+2])/b[i];
  }

  // Storing x and v in and array so we can return it:
  std::vector<std::vector<double> > xv(2);
  xv[0] = x;
  xv[1] = v;

  return xv;
}

std::vector<std::vector<double> > abs_error(std::vector<double> x,
                                  std::vector<double> v, double (u_func)(double)){
  /*
  Takes: vector x, vector v and function u
  Returns: 2D array with x and the absolute difference between u and v
  */
  std::vector<double> abs_error(x.size()-2);
  std::vector<double> new_x(x.size()-2);

  for (int i = 0; i<x.size()-2; i++){
    abs_error[i] = abs(v[i+1] - u_func(x[i+1]));
    new_x[i] = x[i+1];
  }
  std::vector<std::vector<double> > x_abs(2);

  x_abs[0] = new_x;
  x_abs[1] = abs_error;

  return x_abs;
}

std::vector<std::vector<double> > rel_error(std::vector<double> x,
                                  std::vector<double> v, double (u_func)(double)){

  std::vector<std::vector<double> > x_rel = abs_error(x,v,u_func);

  /*
  Takes: vector x, vector v and function u
  Returns: 2D array with x and the relative difference between u and v
  */

  for (int i=0; i<x_rel[1].size(); i++){
    x_rel[1][i] /= abs(u_func(x_rel[0][i]));
  }

  return x_rel;
}

void load_to_file(std::string filename, std::vector<double> x, std::vector<double> y){
  /*
  Prints vector x and y to file with name "filename" 
  */

  // Create file:
  std::ofstream outfile ("text_files/"+filename);

  outfile << std::setw(12) << "x" << std::setw(12) << "y(x)" << std::endl;

  for (int i = 0; i<x.size(); i++){
    outfile << std::setw(12) << std::setprecision(4) << std::scientific << x[i]
            << std::setw(12) << std::setprecision(4) << std::scientific << y[i]
            << std::endl;
  }

  outfile.close();
}
