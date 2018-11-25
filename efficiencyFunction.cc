#include"efficiencyFunction.h"
#include"BELLE_efficiency.h"
#include<iostream>
#include<limits>
#include"constants.h"
#include"math.h"
efficiencyFunction::efficiencyFunction() :_kin1max(std::numeric_limits<double>::infinity()) ,_kinSignature(std::make_shared<kinematicSignature>(0)) {}

double efficiencyFunction::eval(const std::vector<double>& kin) const {
	(void) kin;
	std::cerr << "efficiencyFunction::operator(): ERROR: Called base class method. Returning zero." << std::endl;
	return 0.;
}

threeParticlPerfectEfficiency::threeParticlPerfectEfficiency (std::shared_ptr<kinematicSignature> kinSig) {
	_kinSignature = kinSig;
}

double threeParticlPerfectEfficiency::eval(const std::vector<double>& kin) const {
	(void) kin;
	return 1.;
}

BELLE_DtoKpipi_efficiency::BELLE_DtoKpipi_efficiency() 
{
	_kinSignature = std::make_shared<kinematicSignature>(2);
}

double BELLE_DtoKpipi_efficiency::eval(const std::vector<double>& kin) const {
	if (kin[1] > _kin1max) {
		return 0.;
	}
	return Efficiency(kin);
}
