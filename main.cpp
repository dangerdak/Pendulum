#include "Pendulum.h"
#include <cmath>

using namespace std;

int main() {
	
	const int n_max = 5000;
	const double theta = 0.1, omega = 0;
	const double lambda = 0.3, m = 1, l = 1, g = 9.8;
	double alpha = lambda / (m * sqrt(l * g));
	double dt = 0.8 * alpha;

	euler(theta, omega, alpha, dt, n_max, g, l, m);
	rk4(theta, omega, alpha, dt, n_max, g, l, m);
	leapfrog(theta, omega, alpha, dt, n_max, g, l, m);
	return 0;
}
