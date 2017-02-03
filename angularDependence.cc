#include"angularDependence.h"
#include<iostream>
angularDependence::angularDependence(kinematicSignature kinSignature, std::string name):
	_kinSignature(kinSignature), _name(name) {}

std::complex<double> angularDependence::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature.nKin()) {
		std::cerr << "angularDependence::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature.nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	std::cerr << "angularDependence::eval(...): ERROR: Calling base-class method. Returning zero.d" << std::endl;
	return std::complex<double>(0.,0.);
}

sameMassZeroS::sameMassZeroS(double fsMass) : angularDependence(kinematicSignature(1), "sameMassZeroS"), _fsMass(fsMass) {}

std::complex<double> sameMassZeroS::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature.nKin()) {
		std::cerr << "sameMassZeroS::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature.nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	return std::complex<double>(1.,0.);
}

sameMassOneP::sameMassOneP(double fsMass) : angularDependence(kinematicSignature(1), "sameMassOneP"), _fsMass(fsMass) {}

std::complex<double> sameMassOneP::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature.nKin()) {
		std::cerr << "sameMassOneP::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature.nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	double retVal = (1 + (kin[1] - _fsMass*_fsMass)/kin[0]) * (kin[0] - kin[1] - 2*kin[2] + 3*_fsMass*_fsMass);
	return std::complex<double>(retVal, 0.);
}

sameMassTwoD::sameMassTwoD(double fsMass) : angularDependence(kinematicSignature(1), "sameMassTwoD"), _fsMass(fsMass) {}

std::complex<double> sameMassTwoD::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature.nKin()) {
		std::cerr << "sameMassTwoD::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature.nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
//	std::cerr << "sameMassTwoD::eval(...): ERROR: Calculations for this wave are not done yet. Returning zero" << std::endl;
	const double s   = kin[0];
	const double s12 = kin[1];
	const double m2  = _fsMass*_fsMass;
	double retVal    = (2*m2*pow(m2-2,2) + (-3* m2*m2 - 6 *m2 * s + s*s)*s12 + 2*(3*m2 - s)*s12*s12 + 7*pow(s12,3))*(m2*m2 - 2*s*s - 2*s*s12 + s12*s12 - 2*m2 * (s+s12))/(648 * s * s * s12);
	return std::complex<double> (retVal, 0.);
}
	
