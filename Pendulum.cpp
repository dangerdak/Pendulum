#include "Pendulum.h"
#include <fstream>
#include <cmath>
using namespace std;

//Single Pendulum Euler Method
void euler(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m) {

	double theta_plus = theta_0, omega = omega_0;
	double t = 0;
	ofstream f("euler.csv");
		
	for(int n = 0; n < n_max; n++) {
		
		f << t << "," << theta_plus << "," << omega << "," << get_energy(theta_plus, omega, m, l, g) << endl;

		t += dt;
		double theta = theta_plus;
		theta_plus = theta + omega * dt;
		omega -= (alpha * omega + theta) * dt;
	
	}	


	f.close();

}

//Single Pendulum Leapfrog Method
void leapfrog(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m) {

	double theta_plus, omega_plus,
		theta_n = theta_0, omega_n = omega_0;
	double t = 0;
	ofstream f("leapfrog.csv");
	
	for(int n = 0; n < n_max; n++) {
		
		f << t << "," << theta_n << "," << omega_n << "," << get_energy(theta_n, omega_n, m, l, g) << endl;
		t += dt;
		
		if (n==0) {
			theta_plus = theta_n + omega_n * dt;
			omega_plus = omega_n - (alpha + theta_n) * dt;
		}
		else {
			double theta_minus = theta_n;
			double omega_minus = omega_n;
			theta_n = theta_plus;
			omega_n = omega_plus;

			theta_plus = theta_minus + 2 * omega_n * dt;
			omega_plus = omega_minus - 2 * (alpha * omega_n + theta_n) * dt; 	
		}
	
	}

	f.close();

}

//Single Pendulum RK4 Method
void rk4(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m) {

	double theta = theta_0, omega = omega_0;
	double t = 0;
	ofstream f("rk4.csv");

	for(int n = 0; n < n_max; n++) {

		f << t << "," << theta << "," << omega << endl;
		
		double k1 = omega * dt;
		double l1 = -(alpha * omega + theta) * dt;
		
		double k2 = (omega + 0.5 * l1) * dt;
		double l2 = -(alpha * (omega + 0.5 * l1) + theta + 0.5 * k1) * dt;
		
		double k3 = (omega + 0.5 * l2) * dt;
		double l3 = -(alpha * (omega + 0.5 * l2) + theta + 0.5 * k2) * dt;
		
		double k4 = (omega + l3) * dt;
		double l4 = -(alpha * (omega + l3) + theta + k3) * dt;

		double k = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		double l = (l1 + 2 * l2 + 2 * l3 + l4) / 6;	

		t += dt;	
		theta += k;
		omega += l;

	}

	f.close();

}

//Energy check 
double get_energy(double theta, double omega, const double g, const double l, const double m) {

	double pe = m * g * l * (1 - cos(theta));
	double ke = 0.5 * g * l * l * m * omega * omega;
	double energy = ke + pe;

	return energy;

} 
