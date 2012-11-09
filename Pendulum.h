void euler(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m);

void leapfrog(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m);

void rk4(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m);

double get_energy(double theta, double omega, const double g, const double l, const double m);

void average_range(const double theta_0, const double omega_0, double alpha, double dt,
	const int n_max, const double g, const double l, const double m);


