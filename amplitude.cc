#include"amplitude.h"
#include<iostream>

amplitude::amplitude(kinematicSignature kinSignature, std::string name):
	_kinSignature(kinSignature), _name(name) {}

std::complex<double> amplitude::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature.nKin()) {
		std::cerr << "amplitude::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature.nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	std::cerr << " amplitude::eval(...): ERROR: Calling eval of the base class. Returning zero." << std::endl;
	return std::complex<double>(0.,0.);
}

threeParticleIsobaricAmplitude::threeParticleIsobaricAmplitude(bool boseSymmetrize, std::string name, std::shared_ptr<massShape> shape, std::shared_ptr<angularDependence> angDep):
	amplitude(angDep->kinSignature(), name), _bose(boseSymmetrize), _massShape(shape), _angDep(angDep) {}


std::complex<double> threeParticleIsobaricAmplitude::evalSingleSymmTerm(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature.nKin()) {
		std::cerr << "threeParticleIsobaricAmplitude::evalSingleSymmTerm(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature.nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	std::vector<size_t> isobarMassIndices = _kinSignature.isobarMassIndices();
	if (isobarMassIndices.size() != 1) {
		std::cerr << "threeParticleIsobaricAmplitude::evalSingleSymmTerm(...): ERROR: Number of isobar masses is not one (" << isobarMassIndices.size() << "). Returning zero" << std::endl;
		return std::complex<double>(0.,0.);
	}
	std::complex<double> retVal = _massShape->eval(kin[isobarMassIndices[0]]) * _angDep->eval(kin);
	return retVal;
}

std::complex<double> threeParticleIsobaricAmplitude::eval(const std::vector<double>& kin) const {
	std::complex<double> retVal(0.,0.);
	if (_bose) {
		std::vector<std::vector<double> > symmTerms = _kinSignature.getBoseSymmetrizedKinematics(kin);
		for (std::vector<double>& boseKin : symmTerms) {
			retVal += evalSingleSymmTerm(boseKin);
		}
		retVal /= symmTerms.size();
	} else {
		retVal += evalSingleSymmTerm(kin);
	}
	return retVal;
}
