#include"angularDependence.h"
#include<iostream>
angularDependence::angularDependence(std::shared_ptr<kinematicSignature> kinSignature, std::string name):
	_kinSignature(kinSignature), _name(name) {}

std::complex<double> angularDependence::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "angularDependence::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	std::cerr << "angularDependence::eval(...): ERROR: Calling base-class method. Returning zero.d" << std::endl;
	return std::complex<double>(0.,0.);
}

sameMassZeroS::sameMassZeroS(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassZeroS"), _fsMass(fsMass) {}

std::complex<double> sameMassZeroS::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "sameMassZeroS::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
//	std::cout << "Called sameMassZeroS with " << kin[0] << " " << kin[1] << " " << kin[2] << " giving 1." << std::endl;
	return std::complex<double>(1.,0.);
}

sameMassOneP::sameMassOneP(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassOneP"), _fsMass(fsMass) {}

std::complex<double> sameMassOneP::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "sameMassOneP::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	double retVal = (1 + (kin[1] - _fsMass*_fsMass)/kin[0]) * (kin[0] - kin[1] - 2*kin[2] + 3*_fsMass*_fsMass);
//	std::cout << "Called sameMassOneP with " << kin[0] << " " << kin[1] << " " << kin[2] << " giving " << retVal << std::endl;
	return std::complex<double>(retVal, 0.);
}

sameMassTwoD::sameMassTwoD(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassTwoD"), _fsMass(fsMass) {}

std::complex<double> sameMassTwoD::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "sameMassTwoD::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
//	std::cerr << "sameMassTwoD::eval(...): ERROR: Calculations for this wave are not done yet. Returning zero" << std::endl;
	const double s   = kin[0];
	const double s12 = kin[1];
	const double s13 = kin[2];
	const double m2  = _fsMass*_fsMass;
//	double retVal    = (2*m2*pow(m2-s,2) + (-3* m2*m2 - 6 *m2 * s + s*s)*s12 + 2*(3*m2 - s)*s12*s12 + 7*pow(s12,3))*(m2*m2 - 2*s*s - 2*s*s12 + s12*s12 - 2*m2 * (s+s12))/(648 * s * s * s12);
	double retVal = - (m2*m2 - 2*s*s - 2*s*s12 + s12*s12 - 2*m2*(s+s12));
	retVal *= 2*pow(m2,3) + m2*m2*(-4*s + 9*s12) + s12*(pow(s-s12,2) + 6*(-s+s12)*s13 + 6*s13*s13) + 2*m2*(s*s + 3*s*s12 - 3*s12*(s12 + 3*s13));
	retVal /= 648*s*s*s12;
	return std::complex<double> (retVal, 0.);
}

sameMassZeroSnonRelativistic::sameMassZeroSnonRelativistic(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassZeroSnonRelativistic"), _fsMass(fsMass) {}

std::complex<double> sameMassZeroSnonRelativistic::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "sameMassZeroSnonRelativistic::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
//	std::cout << "Called sameMassZeroS with " << kin[0] << " " << kin[1] << " " << kin[2] << " giving 1." << std::endl;
	return std::complex<double>(1.,0.);
}

sameMassOnePnonRelativistic::sameMassOnePnonRelativistic(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassOnePnonRelativistic"), _fsMass(fsMass) {}

std::complex<double> sameMassOnePnonRelativistic::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "sameMassOnePnonRelativistic::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
//	std::cerr << "sameMassTwoD::eval(...): ERROR: Calculations for this wave are not done yet. Returning zero" << std::endl;
	const double s   = kin[0];
	const double s12 = kin[1];
	const double s13 = kin[2];
	const double m2  = _fsMass*_fsMass;

	double retVal    = -1./4.; // Use definitions as in the paper
//	double retVal    = - pow(s12/s,.5)/4.; //  Use definitions as by Dima
	retVal *= s - s12 - 2*s13 + 3*m2;
	return std::complex<double> (retVal, 0.);
}

sameMassTwoDnonRelativistic::sameMassTwoDnonRelativistic(double fsMass) : angularDependence(std::make_shared<kinematicSignature>(1), "sameMassTwoDnonRelativistic"), _fsMass(fsMass) {}

std::complex<double> sameMassTwoDnonRelativistic::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "sameMassTwoDnonRelativistic::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
//	std::cerr << "sameMassTwoD::eval(...): ERROR: Calculations for this wave are not done yet. Returning zero" << std::endl;
	const double s   = kin[0];
	const double s12 = kin[1];
	const double s13 = kin[2];
	const double m2  = _fsMass*_fsMass;
	double retVal = 0.;
	retVal += -(-4*m2 + s12)*(m2*m2 + pow(s-s12,2) - 2* m2 * (s+s12));
	retVal += 3*s12*pow(3*m2 + s - s12 - 2*s13,2);
	retVal /= 48*s;
	return std::complex<double> (retVal, 0.);
}

ratioOfDependences::ratioOfDependences(std::shared_ptr<angularDependence> numerator, std::shared_ptr<angularDependence> denominator) : 
	angularDependence(numerator->kinSignature(), "ratioOfDependences"), _numerator(numerator), _denominator(denominator) 
{
		if (not (*(numerator->kinSignature()) == *(denominator->kinSignature()))) {
			std::cout << "ratioOfDependences::ratioOfDependences(...): ERROR: Different kinematic signatures found" << std::endl;
		}
}

std::complex<double> ratioOfDependences::eval(const std::vector<double>& kin) const {
	std::complex<double> den = _denominator->eval(kin);
	if (den == std::complex<double>(0.,0.) ) {
		return std::complex<double> (0.,0.);
	}
	return _numerator->eval(kin)/den;
}

arbitraryMass_S_nonRelativistic::arbitraryMass_S_nonRelativistic(size_t isobarIndex, std::vector<double> fsMasses) : angularDependence(std::make_shared<kinematicSignature>(2), "arbitraryMass_S_nonRelativistic"), _isobarIndex(isobarIndex), _fsMasses(fsMasses) {
	if (_fsMasses.size() != 3) {
		std::cout << "arbitraryMass_S_nonRelativistic(...): ERROR: Number of final state paricles masses has to be 3" << std::endl;
		throw;
	}
	if (isobarIndex != 12 and isobarIndex != 13 and isobarIndex != 23) {
		std::cout <<  "arbitraryMass_S_nonRelativistic(...): ERROR: None of the three possible isobar masses match (12, 13 and 23) the given value:" << isobarIndex << std::endl;
		throw;
	}

};

std::complex<double> arbitraryMass_S_nonRelativistic::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "arbitraryMass_S_nonRelativistic::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	return std::complex<double>(1.,0.);
}

arbitraryMass_P_nonRelativistic::arbitraryMass_P_nonRelativistic(size_t isobarIndex, std::vector<double> fsMasses) : angularDependence(std::make_shared<kinematicSignature>(2), "arbitraryMass_S_nonRelativistic"), _isobarIndex(isobarIndex), _fsMasses(fsMasses) {
	if (_fsMasses.size() != 3) {
		std::cout << "arbitraryMass_P_nonRelativistic(...): ERROR: Number of final state paricles masses has to be 3" << std::endl;
		throw;
	}
	if (isobarIndex != 12 and isobarIndex != 13 and isobarIndex != 23) {
		std::cout <<  "arbitraryMass_P_nonRelativistic(...): ERROR: None of the three possible isobar masses match (12, 13 and 23) the given value:" << isobarIndex << std::endl;
		throw;
	}

};

std::complex<double> arbitraryMass_P_nonRelativistic::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "arbitraryMass_P_nonRelativistic::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	const double& s   = kin[0];
	const double& s12 = kin[1];
	const double& s13 = kin[2];

	const double& m1 = _fsMasses[0];
	const double& m2 = _fsMasses[1];
	const double& m3 = _fsMasses[2];

	const double s23 = s + m1*m1 + m2*m2 + m3*m3 - s12 - s13;
	if (_isobarIndex == 12) {
		return (s23 - s13 - (s - m3*m3) * (m1*m1 - m2*m2)/s12)/4;
	}
	if (_isobarIndex == 13) {
		return (s12 - s23 - (s - m2*m2) * (m3*m3 - m1*m1)/s13)/4;
	}
	if (_isobarIndex == 23) {
		return (s13 - s12 - (s - m1*m1) * (m2*m2 - m3*m3)/s23)/4;
	}
	std::cout << "arbitraryMass_P_nonRelativistic::eval(...): ERROR: None of the isobarCombinations matched. Return 0." << std::endl;
	return std::complex<double>(0.,0.);
}


