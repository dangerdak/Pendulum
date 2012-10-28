#include "Pendulum.h"
#include <fstream>
#include <cmath>
using namespace std;

//Single Pendulum Euler Method
void euler(const double theta_0, const double omega_0, const int n_max,
	const double dt, const double lambda, const double m, const double l, 
	const double g) {

	const double alpha = (lambda / m) * sqrt(1 / (l * g));
	double theta = theta_0, omega = omega_0;
	double t = 0;
	ofstream f("euler.csv");
		
	for(int n = 0; n < n_max; n++) {
		
		f << t << "," << theta << "," << omega << endl;
		t += dt;
		theta += omega * dt;
		omega -= (alpha + theta) * dt;

	}	
}

//Single Pendulum Leapfrog Method

