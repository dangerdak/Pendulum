#include "Pendulum.h"

using namespace std;

int main() {
	
	const int n_max = 6000;
	const double theta = 0.1, omega = 0;
	const double dt = 0.01, lambda = 0, m = 1, l = 1, g = 9.8;
	rk4(theta, omega, n_max, dt, lambda, m, l, g);

	return 0;
}
