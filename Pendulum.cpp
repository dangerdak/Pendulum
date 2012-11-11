#include "Pendulum.h"
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

//Single Pendulum Euler Method
void euler(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m) {

	double theta_plus = theta_0, omega = omega_0;
	double t = 0;
	ofstream f("euler.csv");
		
	for(int n = 0; n < n_max; n++) {
		
		f << t << "," << theta_plus << "," << omega << "," << get_energy(theta_plus, omega, m, l, g) << "," << theta_0 * cos(t) << endl;

		euler_update(theta_plus, omega, alpha, t, dt);
	}	


	f.close();

}
//Euler update
void euler_update(double &theta_plus, double &omega, double alpha, double &t, double dt) {
	t +=dt;
	double theta = theta_plus;
	theta_plus = theta + omega * dt;
	omega -= (alpha * omega + theta) * dt;
}

//Allows energy to be plotted as a function of step size and time, but in damped case energy range is too big to see detail of oscillations
void step_energy(double theta, double omega, const int m, const int l, const double g, double alpha, double dt, const int t_max) {
	int h_max = 7;
	ofstream stepfile("euler_step.csv");
	
	for(int h = 0; h <h_max; h++) {
		double omega_plus = omega;
		double theta_plus = theta;
		double t = 0.0;

		double n_max_db = 0.5 + t_max / dt; 
		int n_max = (int) n_max_db;

		for(int i = 0; i < n_max; i++) {
			euler_update(theta_plus, omega_plus, alpha, t, dt);
			stepfile << t << "," << dt << "," << get_energy(theta_plus, omega_plus, g, l, m) << endl;
		}
	dt += 0.001;	
	}
	stepfile.close();

}

//Single Pendulum Leapfrog Method
void leapfrog(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const int l, const int m) {

	double theta_plus, omega_plus,
		theta_n = theta_0, omega_n = omega_0;
	double t = 0;
	ofstream f("leapfrog.csv");
	
	for(int n = 0; n < n_max; n++) {
		
		f << t << "," << theta_n << "," << omega_n << "," << get_energy(theta_n, omega_n, g, l, m) << endl;
		
		leapfrog_update(theta_n, theta_plus, omega_n, omega_plus, alpha, t, dt, n);

	}

	f.close();

}

//Leapfrog update
void leapfrog_update(double &theta_n, double &theta_plus, double &omega_n, double &omega_plus, double alpha, double &t, double dt, int n) {
	
	t += dt;
		
	if (n==0) {
		theta_plus = theta_n + omega_n * dt;
		omega_plus = omega_n - (alpha * omega_n + theta_n) * dt;
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

//Single Pendulum RK4 Method
void rk4(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m) {

	double theta = theta_0, omega = omega_0;
	double t = 0;
	ofstream f("rk4.csv");

	for(int n = 0; n < n_max; n++) {

		f << t << "," << theta << "," << omega <<  "," << get_energy(theta, omega, g, l, m) << "," << theta_0 * cos(t) << endl;
		
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


double get_energy(double theta, double omega, const double g, const double l, const double m) {

	double pe = m * g * l * (1 - (1 - theta*theta/2));
	double ke = 0.5 * g * l * l * m * omega * omega;
	double energy = ke + pe;

	return energy;

} 

//Find maximum in section of array
double maximum(double *array, int start, int finish) {
	double max = array[start];

	for(int n = start + 1; n < finish; n++)
		if (max < array[n])
			max = array[n];
	return max;
}

//Find minimum in section of array
double minimum(double *array, int start, int finish) {
	double min = array[start];

	for(int n = start + 1; n < finish; n++)
		if (min > array[n])
			min = array[n];
	return min;
}

//Find moving energy range for eulers method
void euler_range(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m) {

	double theta_plus = theta_0, omega = omega_0;
	double t = 0;
	double *energy = new double[n_max];	

	for(int n = 0; n < n_max; n++) {
		
		energy[n] = get_energy(theta_plus, omega, m, l, g);

		t += dt;
		double theta = theta_plus;
		theta_plus = theta + omega * dt;
		omega -= (alpha * omega + theta) * dt;
	
	}	

	ofstream f("energy_range.csv");
	int window = n_max / 100;
	double *range = new double[n_max - window];

	for(int i = 0; i < n_max - window; i++) {

		range[i] = maximum(energy, i, i + window - 1) - minimum(energy, i, i + window - 1);
		
		f << range[i] << endl;
	}
	
	f.close();

}

//Double Pendulum
void dbpend(const double theta_0, double alpha, const double dt, const double n_max, const double M, const double m) {
	
	const double phi_0 = 0, omega_0 = 0, v_0 = 0;
	const double R = M/m;
	double k[4];
	double theta = theta_0;
	double phi = phi_0;
	double omega = omega_0;
	double v = v_0;
	double t = 0;

	ofstream f("dbpend.csv");
	//ofstream file("kvalues.csv");

	for(int n = 1; n < n_max; n++) {	
		
		f << t << "," << theta << "," << phi << "," << omega << "," << v << endl;

		fill_k(theta, phi, omega, v, k, alpha, dt, R);
	
		//debug - checking k values
		//file << k[0] << "," << k[1] << "," << k[2] << "," << k[3] << endl;

		theta += k[0];
		phi += k[1];
		omega += k[2];
		v += k[3];
		t += dt;

	}

	f.close();
	//file.close();

}

	

//Funtion to return k values (k_theta, k_phi, k_omega, k_v) ie update values
void fill_k(double theta, double phi, double omega, double v, double k[], double alpha, const double dt, const double R) {

	double k1[4], k2[4], k3[4], k4[4];

	k1[0] = dt * omega;
	k1[1] = dt * v;
	k1[2] = -dt * ( (R + 1) * theta - R * phi + alpha * omega);
	k1[3] = dt * ( (R + 1) * theta - (R + 1) * phi + 
			alpha * (1 - 1 / R) * omega - alpha * v / R);

	k2[0] = dt * (omega + 0.5 * k1[2]);
	k2[1] = dt * (v + 0.5 * k1[3]);
	k2[2] = -dt * ( (R + 1) * (theta + 0.5 * k1[0]) - R * (phi + 0.5 * k1[1]) + 
			alpha * (omega + 0.5 * k1[2]) );
	k2[3] = dt * ( (R + 1) * (theta + 0.5 * k1[0]) - (R + 1) * (phi + 0.5 * k1[1]) + 
			alpha * (1 - 1 / R) * (omega + 0.5 * k1[2]) - alpha * (v + 0.5 * k1[3]) / R);

	k3[0] = dt * (omega + 0.5 * k2[2]);
	k3[1] = dt * (v + 0.5 * k2[3]);
	k3[2] = -dt * ( (R + 1) * (theta + 0.5 * k2[0]) - R * (phi + 0.5 * k2[1]) + 
			alpha * (omega + 0.5 * k2[2]) );
	k3[3] = dt * ( (R + 1) * (theta + 0.5 * k2[0]) - (R + 1) * (phi + 0.5 * k2[1]) + 
			alpha * (1 - 1 / R) * (omega + 0.5 * k2[2]) - alpha * (v + 0.5 * k2[3]) / R);
	
	k4[0] = dt * (omega + k3[0]);
	k4[1] = dt * (v + k3[3]);
	k4[2] = -dt * ( (R + 1) * (theta + k3[0]) - R * (phi + k3[1]) + 
			alpha * (omega + k3[2]) );
	k4[3] = dt * ( (R + 1) * (theta + k3[0]) - (R + 1) * (phi + k3[1]) + 
			alpha * (1 - 1 / R) * (omega + k3[2]) - alpha * (v + k3[3]) / R);

	for (int n=0; n < 4; n++)
		k[n] = (k1[n] + 2 * k2[n] + 2 * k3[n] + k4[n]) / 6;
}






