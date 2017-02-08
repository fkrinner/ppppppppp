#include"massShape.h"

#include<iostream>

massShape::massShape(std::string name, std::vector<double> pars, std::vector<std::string> parNames):
	_nPar(pars.size()), _name(name), _parameters(pars) {
	if (parNames.size() == _nPar) {
		_parameterNames = parNames;
	} else {
		if (parNames.size() > 0) {
			std::cerr << "massShape::massShape(...): ERROR: Number of parameter names does not match and is not zero." << std::endl;
		}
		_parameterNames = std::vector<std::string>();
		for (size_t i = 0; i < _nPar; ++i) {
			_parameterNames.push_back(std::string("unname_parameter_")+std::to_string(i));
		}
	}
}

std::complex<double> massShape::eval(double m) const {
	std::cerr << "massShape::eval(double): ERROR: Evaluating the bas class at mass m = " << m << ", returning zero." << std::endl;
	return std::complex<double>(0.,0.);
}

bool massShape::setParameter(size_t n, double val) {
	if (n < _nPar) {
		_parameters[n] = val;
		return true;
	}
	return false;
}
	
std::pair<bool, double> massShape::getParameter(size_t n) const {
	if (n<_nPar){
		return std::pair<bool, double>(true, _parameters[n]);	
	}
	return std::pair<bool, double>(false, 0.);
}

bool massShape::setParName(size_t n, std::string name) {
	if(n<_nPar) {
		_parameterNames[n] = name;
		return true;
	}
	return false;
}

std::pair<bool, std::string> massShape::getParName(size_t n) const {
	if (n<_nPar) {
		return std::pair<bool, std::string>(true, _parameterNames[n]);
	}
	return std::pair<bool, std::string>(false, "");
}
//------------------------------------------------------------------------------
simpleBW::simpleBW(double mass, double width):
	massShape("simpleBW", {mass, width}, {"m0","G0"}) {}

std::complex<double> simpleBW::eval(double s) const {
	double num = _parameters[0] * _parameters[1];
	std::complex<double> den = std::complex<double>(_parameters[0]*_parameters[0] - s, -_parameters[0]*_parameters[1]);
	return num/den;
}
//------------------------------------------------------------------------------
stepLike::stepLike(double sMin, double sMax) :
	massShape("stepLike", {sMin, sMax}, {"sMin","sMax"}) {}

std::complex<double> stepLike::eval(double s) const {
	if (s > _parameters[0] and s <= _parameters[1]) {
		return std::complex<double>(1.,0.);
	}
	return std::complex<double>(0.,0.);
}
//------------------------------------------------------------------------------
constant::constant() :
	massShape("constant", {}, {}) {};

std::complex<double> constant::eval(double s) const {
	return std::complex<double>(s/s,0.);
}
//------------------------------------------------------------------------------
zeroMode0pp::zeroMode0pp(double s, double m2) :
	massShape("zeroMode0pp", {s, m2}, {"s", "m2"}) {};

std::complex<double> zeroMode0pp::eval(double s12) const {
	double retVal = (_parameters[0] - 3*s12 + 3*_parameters[1])/8;
//	std::cout << "Called zeroMode0pp with " << s12 << " giving " << retVal << std::endl;
	return std::complex<double>(retVal, 0.);
}
//------------------------------------------------------------------------------
zeroMode1mm::zeroMode1mm(double s, double m2) :
	massShape("zeroMode0pp", {s, m2}, {"s", "m2"}) {};

std::complex<double> zeroMode1mm::eval(double s12) const {
	double retVal = _parameters[0]/(_parameters[0] + s12 - _parameters[1]);
//	std::cout << "Called zeroMode1mm with " << s12 << " giving " << retVal << std::endl;
	return std::complex<double>(retVal, 0.);
}

