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

	f.close();

}

//Single Pendulum Leapfrog Method
void leapfrog(const double theta_0, const double omega_0, const int n_max,
	const double dt, const double lambda, const double m, const double l, 
	const double g) {

	const double alpha = (lambda / m) * sqrt(1 / (l * g));
	double theta_plus, theta_minus, omega_plus, omega_minus,
		theta_n = theta_0, omega_n = omega_0;
	double t = 0;
	ofstream f("leapfrog.csv");
	
	for(int n = 0; n < n_max; n++) {
		
		f << t << "," << theta_n << "," << omega_n << endl;
		t += dt;
		
		if (n==0) {
			theta_plus = theta_n + omega_n * dt;
			omega_plus = omega_n - (alpha + theta_n) * dt;
		}
		else {
			theta_minus = theta_n;
			omega_minus = omega_n;
			theta_n = theta_plus;
			omega_n = omega_plus;

			theta_plus = theta_minus + 2 * omega_n * dt;
			omega_plus = omega_minus - 2 * (alpha * omega_n + theta_n) * dt; 	
		}
	
	}

	f.close();

}
