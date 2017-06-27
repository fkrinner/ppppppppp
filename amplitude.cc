#include"amplitude.h"
#include<iostream>

amplitude::amplitude(std::shared_ptr<kinematicSignature> kinSignature, std::string name):
	_kinSignature(kinSignature), _name(name) {}

std::complex<double> amplitude::eval(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "amplitude::eval(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	std::cerr << " amplitude::eval(...): ERROR: Calling eval of the base class. Returning zero." << std::endl;
	return std::complex<double>(0.,0.);
}

threeParticleIsobaricAmplitude::threeParticleIsobaricAmplitude(bool boseSymmetrize, std::string name, std::shared_ptr<massShape> shape, std::shared_ptr<angularDependence> angDep):
	amplitude(angDep->kinSignature(), name), _bose(boseSymmetrize), _massShape(shape), _angDep(angDep) {}


std::complex<double> threeParticleIsobaricAmplitude::evalSingleSymmTerm(const std::vector<double>& kin) const {
//	std::cout << "Called single term with " << kin[0] << " " << kin[1] << " " << kin[2] << std::endl;

	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "threeParticleIsobaricAmplitude::evalSingleSymmTerm(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}
	std::vector<size_t> isobarMassIndices = _kinSignature->isobarMassIndices();
	if (isobarMassIndices.size() != 1) {
		std::cerr << "threeParticleIsobaricAmplitude::evalSingleSymmTerm(...): ERROR: Number of isobar masses is not one (" << isobarMassIndices.size() << "). Returning zero" << std::endl;
		return std::complex<double>(0.,0.);
	}
	std::complex<double> retVal = _massShape->eval(kin[isobarMassIndices[0]]) * _angDep->eval(kin);
	return retVal;
}

std::complex<double> threeParticleIsobaricAmplitude::eval(const std::vector<double>& kin) const {
//	std::cout << "------------------------------------------------------------------" << std::endl;
//	std::cout << "Called amplitude with " << kin[0] << " " << kin[1] << " " << kin[2] << std::endl;
	std::complex<double> retVal(0.,0.);
	if (_bose) {
		std::vector<std::vector<double> > symmTerms = _kinSignature->getBoseSymmetrizedKinematics(kin);
		for (std::vector<double>& boseKin : symmTerms) {
			retVal += evalSingleSymmTerm(boseKin);
		}
		retVal /= symmTerms.size();
	} else {
		retVal += evalSingleSymmTerm(kin);
	}
	return retVal;
}

threeParticlaIsobaricAmplitudeNoBose::threeParticlaIsobaricAmplitudeNoBose(size_t isobarIndex, std::string name, std::shared_ptr<massShape> shape, std::shared_ptr<angularDependence> angDep, std::vector<double> fsMasses):
	amplitude(angDep->kinSignature(), name), _isobarIndex(isobarIndex), _sumFSmasses(0.),  _massShape(shape), _angDep(angDep) {
	
	if (fsMasses.size() != 3) {
		std::cout << "threeParticlaIsobaricAmplitudeNoBose::threeParticlaIsobaricAmplitudeNoBose(...): ERROR: Number of final state masses differs from three" << std::endl;
		throw;
	}
	for (double& fsMass : fsMasses) {
		_sumFSmasses += fsMass*fsMass;
	}
	if (isobarIndex != 12 and isobarIndex != 13 and isobarIndex != 23) {
		std::cout <<  "threeParticlaIsobaricAmplitudeNoBose::threeParticlaIsobaricAmplitudeNoBose(...): ERROR: None of the three possible isobar masses match (12, 13 and 23) the given value:" << isobarIndex << std::endl;
		throw;
	}
}

std::complex<double> threeParticlaIsobaricAmplitudeNoBose::eval(const std::vector<double>& kin) const {
	std::complex<double> retVal = _angDep->eval(kin);
	if (_isobarIndex == 12) {

	}
	double mIsob = 0.;
	if (_isobarIndex == 12) {
		mIsob = kin[1];
	} else if (_isobarIndex == 13) {
		mIsob = kin[2];
	} else if (_isobarIndex == 23) {
		mIsob = kin[0] + _sumFSmasses - kin[1] - kin[2];
	} else {
		std::cout << "threeParticlaIsobaricAmplitudeNoBose::eval(...): ERROR: Invalid isobarIndex: " << _isobarIndex << ". Returning 0" << std::endl;
		return std::complex<double>(0.,0.);
	} 
	retVal *= _massShape->eval(mIsob);
	return retVal;
}





