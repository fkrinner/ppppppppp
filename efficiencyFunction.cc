#include"efficiencyFunction.h"
#include"BELLE_efficiency.h"
#include<iostream>
#include"constants.h"
#include"math.h"
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

BELLE_DtoKpipi_efficiency::BELLE_DtoKpipi_efficiency(): _minM2Pisquared(0.), _maxAbsCosT(1.) {
	_kinSignature = std::make_shared<kinematicSignature>(2);
}

double BELLE_DtoKpipi_efficiency::call(const std::vector<double>& kin) const {
	if (kin[2] < _minM2Pisquared) {
		double p1pK = (kin[1] - mPi*mPi - mKs*mKs)/2;
		double Epi  = pow(kin[2],.5)/2.; // Half the isobar mass in the isobar rest frame
		double EKs  = (mD0*mD0 - kin[2] - mKs*mKs)/4/Epi; // 4*Epi = 2*mPiPi = 2*pow(kin[2],.5)
		double pPi  = pow(Epi*Epi - mPi*mPi, .5);
		double pKs  = pow(EKs*EKs - mKs*mKs, .5);

		double cosT = (Epi*EKs - p1pK)/pPi/pKs;

		if (abs(cosT) > _maxAbsCosT) {
			return 0.;
		}
	}
	return Efficiency(kin);
}
