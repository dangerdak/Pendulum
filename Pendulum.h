//implement euler method for single pendulum
void euler(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m);

//implement leapfrog method for single pendulum
void leapfrog(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m);

//implement rk4 method for single pendulum
void rk4(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m);

//find change in energy of single pendulum over time
//(as predicted by any of the above methods)
double get_energy(double theta, double omega, const double g, const double l, const double m);

//find range of oscillations in energy (doesn't work for damped case)
void euler_range(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m);

//implement rk4 for a double pendulum
void dbpend(const double theta_0, double alpha, const double dt, const double n_max, const double M, 
		const double m);

//find update (k) values for double pendulum rk4 method
void fill_k(double theta, double phi, double omega, double v, double k[], 
		const double alpha, const double dt, const double R);


