#include"amplitude.h"
#include<iostream>
#include<fstream>
#include<math.h>

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
	double sIsob = 0.;
	if (_isobarIndex == 12) {
		sIsob = kin[1];
	} else if (_isobarIndex == 13) {
		sIsob = kin[2];
	} else if (_isobarIndex == 23) {
		sIsob = kin[0] + _sumFSmasses - kin[1] - kin[2];
	} else {
		std::cout << "threeParticlaIsobaricAmplitudeNoBose::eval(...): ERROR: Invalid isobarIndex: " << _isobarIndex << ". Returning 0" << std::endl;
		return std::complex<double>(0.,0.);
	} 
//	if (sIsob < 0.) {
//		std::cout << _isobarIndex << " " << kin[0] << " " << kin[1] << " " << kin[2] << " " << _sumFSmasses << std::endl;
//	}
	retVal *= _massShape->eval(sIsob);
	return retVal;
}


dalitzMonomialAmplitude::dalitzMonomialAmplitude(std::shared_ptr<kinematicSignature> kinSignature, double exponent1, double exponent2):
	amplitude(kinSignature, std::string("monomial_")+std::to_string(exponent1)+std::string("_")+std::to_string(exponent2)), _exponent1(exponent1), _exponent2(exponent2) {}

std::complex<double> dalitzMonomialAmplitude::eval(const std::vector<double>& kin) const {
	return std::complex<double>(pow(kin[1],_exponent1) * pow(kin[2], _exponent2));
}

dalitzPolynomialAmplitude::dalitzPolynomialAmplitude(std::shared_ptr<kinematicSignature> kinSignature, const std::string configurationFile, double xMin, double xMax, double yMin, double yMax) :
	amplitude(kinSignature, std::string("dalitzPolynomialAmplitude_from_")+configurationFile), 
	_nTerms(0),
	_xMin(xMin), _xMax(xMax), _xWidth(xMax-xMin), _yMin(yMin), _yMax(yMax), _yWidth(yMax-yMin),
	_xExponents(),
	_yExponents(),
	_coefficients() {

	std::ifstream fin(configurationFile.c_str(), std::ifstream::in);
	size_t degX;
	size_t degY;
	double coeff;
	while (true) {
		if (!(fin>>degX)) {
			break;
		}
		if (!(fin>>degY)) {
			std::cout << "mCosTintensPolynomial(...): ERROR: Invalid configuration file '" << configurationFile << "' (degY)" << std::endl;
			throw;
		}
		if (!(fin>>coeff)) {
			std::cout << "mCosTintensPolynomial(...): ERROR: Invalid configuration file '" << configurationFile << "' (coeff)" << std::endl;
			throw;
		}
		_xExponents.push_back(degX);
		_yExponents.push_back(degY);
		_coefficients.push_back(coeff);
	}
	_nTerms = _coefficients.size();
}

std::complex<double> dalitzPolynomialAmplitude::eval(const std::vector<double>& kin) const {
	double retVal = 0.;
	std::pair<double, double> XY = getXY(kin);
	for (size_t c = 0; c < _nTerms; ++c) {
		retVal += _coefficients[c] *pow(XY.first, _xExponents[c]) * pow(XY.second, _yExponents[c]);
	}
	return std::complex<double>(pow(retVal*retVal, .25), 0.);
}

std::pair<double, double> dalitzPolynomialAmplitude::getXY(const std::vector<double>& kin) const {
	double X = 2*(kin[1]-_xMin)/_xWidth-1.;
	double Y = 2*(kin[2]-_yMin)/_yWidth-1.;
	return std::pair<double, double>(X,Y);
}

mCosTintensPolynomial::mCosTintensPolynomial(std::shared_ptr<kinematicSignature> kinSignature, const std::string configurationFile, double motherMass, std::vector<double> fsMasses, size_t isobarCombination) :
	amplitude(kinSignature, std::string("mCosTintensPolynomial_from_")+configurationFile), 
	_isobarCombination(isobarCombination), 
	_nTerms(0), 
	_mWidth(0.),
	_cosTwidth(2.),
	_mLimits(0.,0.),
	_cosTlimits(-1.,1.),
	_xExponents(),
	_yExponents(),
	_coefficients(),
	_fsMasses(fsMasses) {
	if (_isobarCombination == 12) {
		_mLimits = std::pair<double, double>(_fsMasses[0] + _fsMasses[1], motherMass - _fsMasses[2]);
	} else if (_isobarCombination == 13) {
		_mLimits = std::pair<double, double>(_fsMasses[0] + _fsMasses[2], motherMass - _fsMasses[1]);
	} else if (_isobarCombination == 23) {
		std::cout << "mCosTintensPolynomial(...): NOT IMPLEMENTED ERROR: _isobarCombination == 23 not implemented yet" << std::endl;
		throw;
	} else {
		std::cout << "mCosTintensPolynomial(...): ERROR: Invalid _isobarCombination: " << _isobarCombination << std::endl;
		throw;
	}
	std::ifstream fin(configurationFile.c_str(), std::ifstream::in);
	size_t degX;
	size_t degY;
	double coeff;
	while (true) {
		if (!(fin>>degX)) {
			break;
		}
		if (!(fin>>degY)) {
			std::cout << "mCosTintensPolynomial(...): ERROR: Invalid configuration file '" << configurationFile << "'" << std::endl;
			throw;
		}
		if (!(fin>>coeff)) {
			std::cout << "mCosTintensPolynomial(...): ERROR: Invalid configuration file '" << configurationFile << "'" << std::endl;
			throw;
		}
		_xExponents.push_back(degX);
		_yExponents.push_back(degY);
		_coefficients.push_back(coeff);
	}
	_nTerms = _coefficients.size();
	_mWidth = _mLimits.second - _mLimits.first;
}

std::complex<double> mCosTintensPolynomial::eval(const std::vector<double>& kin) const {
	double retVal = 0.;
	std::pair<double, double> XY = getXY(getMcosT(kin));
	if (isnan(XY.first) || isnan(XY.second)) {
		std::cout << "mCosTintensPolynomial::eval(...): WARNING: NaN encountered for kin = {" << kin[0];
		for (size_t k = 1; k < kin.size(); ++k) {
			std::cout << ", " << kin[k];
		}
		std::cout << "}: return (0.,0.)" << std::endl;
		return std::complex<double>(0.,0.);
		
	}
	for (size_t c = 0; c < _nTerms; ++c) {
		retVal += _coefficients[c] *pow(XY.first, _xExponents[c]) * pow(XY.second, _yExponents[c]);
	}
	return std::complex<double>(pow(retVal*retVal, .25), 0.);
}

std::pair<double, double> mCosTintensPolynomial::getMcosT(const std::vector<double>& kin) const {
	double isobarMass  = 0.;
	double isobarMass2 = 0.;
	double fsm1        = 0.;
	double fsm2        = 0.;
	if (_isobarCombination == 12) { 
		isobarMass  = pow(kin[1],.5);
		isobarMass2 = pow(kin[2],.5);
		fsm1        = _fsMasses[0];
		fsm2        = _fsMasses[1];
	} else if (_isobarCombination == 13) {
		isobarMass  = pow(kin[2],.5);
		isobarMass2 = pow(kin[1],.5);
		fsm1        = _fsMasses[0];
		fsm2        = _fsMasses[2];
	}

	double p1p2 = (isobarMass2*isobarMass2 - fsm1*fsm1 - fsm2*fsm2)/2;
	double E1   = isobarMass/2.;
	double E2   = (kin[0] - isobarMass*isobarMass - fsm2*fsm2)/4/E1;
	double p1   = pow(E1*E1 - _fsMasses[0]*_fsMasses[0],.5);
	double p2   = pow(E2*E2 - _fsMasses[1]*_fsMasses[1],.5);
	double cosT = (E1*E2 - p1p2)/p1/p2;


	return std::pair<double, double>(isobarMass, cosT);
}

std::pair<double, double> mCosTintensPolynomial::getXY(const std::pair<double, double> mCosT) const {
	return std::pair<double,double>(2*(mCosT.second - _cosTlimits.first)/_cosTwidth-1.,2*(mCosT.first-_mLimits.first)/_mWidth-1.);
}

bool mCosTintensPolynomial::setCosTLimits(const std::pair<double,double> newLimits) {
	if (newLimits.first >= newLimits.second) {
		std::cout << "mCosTintensPolynomial::setCosTLimits(...): ERROR: Limits not ordered" << std::endl;
		return false;
	}
	_cosTlimits = newLimits;
	_cosTwidth  = _cosTlimits.second - _cosTlimits.first;
	return true;
}

bool mCosTintensPolynomial::setMlimits(const std::pair<double,double> newLimits) {
	if (newLimits.first >= newLimits.second) {
		std::cout << "mCosTintensPolynomial::setMlimits(...): ERROR: Limits not ordered" << std::endl;
		return false;
	}
	_mLimits = newLimits;
	_mWidth = _mLimits.second - _mLimits.first;
	return true;
}






