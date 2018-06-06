#ifndef BELLE_EFFICIENCY
#define BELLE_EFFICIENCY
#include<vector>

// BELLE_D0ToK0Spipi_efficiency_correction

double compute_legendre_polynomial(double x, int legendre_order);

double Efficiency(double *x);

double Efficiency(std::vector<double> x) {
	return Efficiency(&x[0]);
}

#endif//BELLE_EFFICIENCY
