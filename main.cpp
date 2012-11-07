#include "Pendulum.h"
#include <cmath>

using namespace std;

int main() {
	
	const int n_max = 500;
	const double theta = 0.1, omega = 0;
	const double gamma = 0.2, m = 1, l = 1, g = 9.8;
	double alpha = gamma / (m * sqrt(l * g));
	double dt = alpha;

	euler(theta, omega, alpha, dt, n_max, g, l, m);
	rk4(theta, omega, alpha, dt, n_max, g, l, m);
	leapfrog(theta, omega, alpha, dt, n_max, g, l, m);
	return 0;
}
