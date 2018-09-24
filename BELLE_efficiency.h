#ifndef BELLE_EFFICIENCY__
#define BELLE_EFFICIENCY__
#include<vector>

// BELLE_D0ToK0Spipi_efficiency_correction

double compute_legendre_polynomial(double x, int legendre_order);

double Efficiency(double *x);

double Efficiency(std::vector<double> x) {
	return Efficiency(&x[1]); // Start at one, since the 0th entry is s = m_D^2
}

#endif//BELLE_EFFICIENCY__
