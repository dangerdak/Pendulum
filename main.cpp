#include "Pendulum.h"
#include <cmath>

using namespace std;

int main() {

	//initialise theta and omega to values at t=0	
	double theta = 0.1, omega = 0.0;
	
	//value for damping constant, gamma and step size, dt
	const double gamma = 0.0;
	double dt = 0.500;
	

	//find number of iterations from total time wanted. Add 0.5 so conversion to int rounds properly.
	const int t_max = 1000;
	double n_max_db = 0.5 + t_max / dt; 
	int n_max = (int) n_max_db;

	const int l = 1;
       	const double g = 1;

	//m is mass of upper pendulum, M is mass of lower pendulum
	//const int M = 1;
	const int m = 1;
	double alpha = gamma / (m * sqrt(l * g));
	
	euler(theta, omega, alpha, dt, n_max, g, l, m);
	rk4(theta, omega, alpha, dt, n_max, g, l, m);
	leapfrog(theta, omega, alpha, dt, n_max, g, l, m);
	//step_energy(theta, omega, m, l, g, alpha, dt, t_max);
	//euler_range(theta, omega, alpha, dt, n_max, g, l, m);
	
	//dbpend(theta, alpha, dt, n_max, M, m); 
	


	return 0;
	
}
