#include"efficiencyFunction.h"
#include<iostream>
efficiencyFunction::efficiencyFunction() : _kinSignature(std::make_shared<kinematicSignature>(0)) {}

double efficiencyFunction::call(const std::vector<double>& kin) const {
	(void) kin;
	std::cerr << "efficiencyFunction::operator(): ERROR: Called base class method. Returning zero." << std::endl;
	return 0.;
}

threeParticlPerfectEfficiency::threeParticlPerfectEfficiency (std::shared_ptr<kinematicSignature> kinSig) {
	_kinSignature = kinSig;
}

double threeParticlPerfectEfficiency::call(const std::vector<double>& kin) const {
	(void) kin;
	return 1.;
}
