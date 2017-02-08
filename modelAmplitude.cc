#include"modelAmplitude.h"
#include<iostream>
modelAmplitude::modelAmplitude(std::vector<std::complex<double> > transitionAmplitudes, std::vector<std::shared_ptr<amplitude> > amplitudes, std::vector<double> normalizations) {
	if (amplitudes.size() == 0) {
		std::cerr << "modelAmplitude::modelAmplitude(...): ERROR: No amplitudes given." << std::endl;
		throw;
	}
	_nAmpl = amplitudes.size();
	_kinSignature = amplitudes[0]->kinSignature();
	for (std::shared_ptr<amplitude> a : amplitudes) {
		if (not ( *(a->kinSignature()) == *_kinSignature)) {
			std::cerr << "modelAmplitude::modelAmplitude(...): ERROR: Different kinematic signatures encountered" << std::endl;
			throw;
		}
	}
	_amplitudes = amplitudes;
	if (amplitudes.size() != transitionAmplitudes.size()) {
		std::cerr << "modelAmplitude::modelAmplitude(...): ERROR: size of trnasitions amplitudes does not match" << std::endl;
		throw;
	}
	_transitionAmplitudes = transitionAmplitudes;
	if (amplitudes.size() != normalizations.size()) {
		std::cerr << "modelAmplitude::modelAmplitude(...): ERROR: size of normalizations does not match" << std::endl;
		throw;		
	}
	_normalizations = normalizations;
}

std::complex<double> modelAmplitude::ampl(const std::vector<double>& kin) const {
	if (kin.size() != _kinSignature->nKin()) {
		std::cerr << "modelAmplitude::ampl(...): ERROR: Number of kinematic variables does not match (" << kin.size() << " != " << _kinSignature->nKin() << "). Returning zero." << std::endl;
		return std::complex<double>(0.,0.);
	}

	std::complex<double> retVal(0.,0.);
	for (size_t a = 0; a < _nAmpl; ++a) {
		retVal += _transitionAmplitudes[a] * _amplitudes[a]->eval(kin) * _normalizations[a];
	}
	return retVal;
}

bool modelAmplitude::setTransitionAmplitude(size_t n, std::complex<double> amp) {
	if (n >= _nAmpl) {
		return false;
	}
	_transitionAmplitudes[n] = amp;
	return true;
}
