#include "Pendulum.h"
#include <cmath>

using namespace std;

int main() {
	
	const int n_max = 400;
	const double theta = 0.1, omega = 0;
	const double lambda = 0.2, m = 1, l = 1, g = 9.8;
	double alpha = lambda / (m * sqrt(l * g));
	double dt = alpha;

	rk4(theta, omega, alpha, dt, n_max, g, l, m);
	euler(theta, omega, alpha, dt, n_max, g, l, m);
	leapfrog(theta, omega, alpha, dt, n_max, g, l, m);

	return 0;
}
