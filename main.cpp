#include "Pendulum.h"
#include <cmath>

using namespace std;

int main() {
	
	const double theta = 0.1, omega = 0;
	const double gamma = 0.0, m = 1, l = 1, g = 9.8;
	double alpha = gamma / (m * sqrt(l * g));
	double dt = 0.1;

	const int t_max = 1000;
	double n_max_db = 0.5 + t_max / dt; //find number of iterations from total time wanted. Add 0.5 so conversion to int rounds properly.
	const int n_max = (int) n_max_db;
	
	euler(theta, omega, alpha, dt, n_max, g, l, m);
	rk4(theta, omega, alpha, dt, n_max, g, l, m);
	leapfrog(theta, omega, alpha, dt, n_max, g, l, m);

	average_range(theta, omega, alpha, dt, n_max, g, l, m);


	return 0;
	
}
