#include"efficiencyFunction.h"
#include"BELLE_efficiency.h"
#include<iostream>
#include<limits>
#include"constants.h"
#include"math.h"
efficiencyFunction::efficiencyFunction() : _kinSignature(std::make_shared<kinematicSignature>(0)) {}

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

BELLE_DtoKpipi_efficiency::BELLE_DtoKpipi_efficiency() {
	_kinSignature = std::make_shared<kinematicSignature>(2);
}

double BELLE_DtoKpipi_efficiency::eval(const std::vector<double>& kin) const {
	return Efficiency(kin);
}

BELLE_DtoKpipi_efficiency_CP::BELLE_DtoKpipi_efficiency_CP(const std::vector<double>& fs_masses) {
	_kinSignature = std::make_shared<kinematicSignature>(2);
	if (fs_masses.size() != 3) {
		std::cout << "BELLE_DtoKpipi_efficiency_CP::BELLE_DtoKpipi_efficiency_CP(...): ERROR: Number of final state masses needs to be three. Given:";
		for (const double& mass : fs_masses) {
			std::cout << " " << mass;
		}
		std::cout << std::endl;
		throw;
	}
	_fs_masses_square_sum = 0.; 
	for (const double& mass : fs_masses) {
		_fs_masses_square_sum += mass*mass;
	}
}

double BELLE_DtoKpipi_efficiency_CP::eval(const std::vector<double>& kin) const {
	std::vector<double> CP_kin(3,0.);
	CP_kin[0] = kin[0];
	CP_kin[1] = kin[0] + _fs_masses_square_sum - kin[1] - kin[2];
	CP_kin[2] = kin[2];

	return Efficiency(CP_kin);
}
