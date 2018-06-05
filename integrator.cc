#include"integrator.h"
#include"utils.h"
#include<iostream>
#include<fstream>
integrator::integrator(size_t integralPoints, std::shared_ptr<generator> pointGenerator, const std::vector<std::shared_ptr<amplitude> >& amplitudes, std::shared_ptr<efficiencyFunction>& efficiency):
	_isIntegrated(false), _nAmpl(amplitudes.size()), _nPoints(integralPoints), _amplitudes(amplitudes), _generator(pointGenerator), _efficiency(efficiency), _integralMatrix(), _accCorrIntegralMatrix() {
	if (_amplitudes.size() == 0) {
		std::cerr << "integrator::integrator(...): ERROR: No amplitudes given" << std::endl;
		throw;
	}
	_kinSignature = _amplitudes[0]->kinSignature();
	for (std::shared_ptr<amplitude> a : _amplitudes) {
		if (not (*(a->kinSignature()) == *_kinSignature)) {
			std::cerr << "integrator::integrator(...): ERROR: Amplitudes with different signatures found" << std::endl;
			throw;
		}
	}
	if (not (*(_generator->kinSignature()) == *(_amplitudes[0]->kinSignature()))) {
		std::cerr << "integrator::integrator(...): ERROR: Kinematic signatures in amplitudes and generator differ" << std::endl;
		throw;
	}
	if (not (*(_efficiency->kinSignature()) == *(_amplitudes[0]->kinSignature()))) {
		std::cerr << "integrator::integrator(...): ERROR: Kinematic signatures in amplitudes and efficiency differ" << std::endl;
		throw;
	}
}

bool integrator::integrate() {
	_integralMatrix        = std::vector<std::vector<std::complex<double> > > (nAmpl(), std::vector<std::complex<double> > (nAmpl(), std::complex<double>(0.,0.)));
	_accCorrIntegralMatrix = std::vector<std::vector<std::complex<double> > > (nAmpl(), std::vector<std::complex<double> > (nAmpl(), std::complex<double>(0.,0.)));
	for (size_t point = 0; point < _nPoints; ++point) {
		std::vector<double> kin = _generator->generate();
		std::vector<std::complex<double> > ampl(nAmpl(), std::complex<double>(0.,0.));
		bool accepted = false;
		if (_efficiency->call(kin) > random()) {
			accepted = true;
		}
		size_t countAmp = 0;
		for (std::shared_ptr<amplitude> a : _amplitudes) {
			ampl[countAmp] = a->eval(kin);
			++countAmp;
		}
		for (size_t i = 0; i < nAmpl(); ++i) {
			for (size_t j = 0; j < nAmpl(); ++j) {
				_integralMatrix[i][j] += std::conj(ampl[i])*ampl[j];
				if (accepted) {
					_accCorrIntegralMatrix[i][j] += std::conj(ampl[i])*ampl[j];
				}
			}
		}
	}
	for (size_t i = 0; i < nAmpl(); ++i) {
		for (size_t j = 0; j < nAmpl(); ++j) {
			_integralMatrix[i][j]        /= _nPoints;
			_accCorrIntegralMatrix[i][j] /= _nPoints;
		}
	}
	_isIntegrated = true;
	return true;
}

std::vector<std::vector<std::complex<double> > > integrator::getIntegralMatrix(bool accCorr) const {
	if (not _isIntegrated) {
		std::cerr << "integrator::getIntegralMatrix(): ERROR: Not integrated yet" << std::endl;
	}
	if (accCorr){
		return _accCorrIntegralMatrix;
	}
	return _integralMatrix;
}

std::pair<bool, std::complex<double> > integrator::element(size_t i, size_t j, bool accCorr) const {
	if (not _isIntegrated) {
		std::cerr << "integrator::element(): ERROR: Not integrated yet" << std::endl;
		return std::pair<bool, std::complex<double> >(false, std::complex<double>(0.,0.));
	}
	if (i >= nAmpl()) {
		std::cerr << "integrator::element(): ERROR: First index too big" << std::endl;
		return std::pair<bool, std::complex<double> >(false, std::complex<double>(0.,0.));
	}
	if (j >= nAmpl()) {
		std::cerr << "integrator::element(): ERROR: Second index too big" << std::endl;
		return std::pair<bool, std::complex<double> >(false, std::complex<double>(0.,0.));
	}
	if (accCorr){
		return std::pair<bool, std::complex<double> >(true, _accCorrIntegralMatrix[i][j]);	
	}
	return std::pair<bool, std::complex<double> >(true, _integralMatrix[i][j]);	
}

double integrator::totalIntensity(const std::vector<std::complex<double> >& prodAmpl, bool accCorr) const {
	if (not _isIntegrated) {
		std::cerr << "integrator::totalIntensity(...): ERROR: Not integrated yet. Returning zero." << std::endl;
		return 0.;
	}
	if (not prodAmpl.size() == nAmpl()) {
		std::cerr << "integrator::totalIntensity(...): ERROR: Number of prodiction amplitudes does not match. Returning zero." << std::endl;
		return 0.;
	}
	const std::vector<std::vector<std::complex<double> > > &integralMatrix = accCorr ? _accCorrIntegralMatrix : _integralMatrix;
	std::complex<double> retVal(0.,0.);
	for (size_t i = 0; i < nAmpl(); ++i) {
		for (size_t j = 0; j < nAmpl(); ++j) {
			double norm = pow(_integralMatrix[i][i].real() * _integralMatrix[j][j].real(), .5); // Always normalized to phaseSpace
			if (norm == 0.) {
				continue;
			}
			retVal += std::conj(prodAmpl[i]) * integralMatrix[i][j] * prodAmpl[j]/norm;			
		}
	}
	return retVal.real();	
}

std::vector<double> integrator::DtotalIntensity(const std::vector<std::complex<double> >& prodAmpl, bool accCorr) const {
// Return a vector of length 2*nAmpl() for the derivative w.r.t. re, and im : {dRe0, dIm0, dRe1, ..., dImnAmpl()}
	if (not _isIntegrated) {
		std::cerr << "integrator::DtotalIntensity(...): ERROR: Not integrated yet. Returning empty vector." << std::endl;
		return std::vector<double>();
	}
	if (not prodAmpl.size() == nAmpl()) {
		std::cerr << "integrator::DtotalIntensity(...): ERROR: Number of prodiction amplitudes does not match. Returning zero." << std::endl;
		return std::vector<double>();
	}
	std::vector<double> retVal(2*nAmpl(), 0.);
	const std::vector<std::vector<std::complex<double> > > &integralMatrix = accCorr ? _accCorrIntegralMatrix : _integralMatrix;
	for (size_t ai = 0; ai < nAmpl(); ++ai) {
		for (size_t aj = 0; aj < nAmpl(); ++aj) {
			double norm = pow(_integralMatrix[ai][ai].real() * _integralMatrix[aj][aj].real(), .5); // Always normalized to phaseSpace
			if (norm == 0.) {
				continue;
			}
			std::complex<double> factor = integralMatrix[ai][aj] * prodAmpl[aj]/norm*2.;
			retVal[2*ai  ] += factor.real();
			retVal[2*ai+1] += factor.imag();
		}
	}
	return retVal;
}

std::vector<std::vector<double> > integrator::DDtotalIntensity(const std::vector<std::complex<double> >& prodAmpl, bool accCorr) const {
	if (not _isIntegrated) {
		std::cerr << "integrator::DDtotalIntensity(...): ERROR: Not integrated yet. Returning empty vector." << std::endl;
		return std::vector<std::vector<double> >();
	}
	if (not prodAmpl.size() == nAmpl()) {
		std::cerr << "integrator::DDtotalIntensity(...): ERROR: Number of prodiction amplitudes does not match. Returning zero." << std::endl;
		return std::vector<std::vector<double> >();
	}
	std::vector<std::vector<double> > retVal(2*nAmpl(), std::vector<double>(2*nAmpl(),0.));
	const std::vector<std::vector<std::complex<double> > > &integralMatrix = accCorr ? _accCorrIntegralMatrix : _integralMatrix;
	for (size_t ai = 0; ai < nAmpl(); ++ai) {
		for (size_t aj = 0; aj < nAmpl(); ++aj) {
			double norm = pow(_integralMatrix[ai][ai].real() * _integralMatrix[aj][aj].real(), .5); // Always normalized to phaseSpace
			if (norm == 0.) {
				continue;
			}
			std::complex<double> factor = integralMatrix[ai][aj]/norm*2.;
			retVal[2*ai  ][2*aj  ] += factor.real();
			retVal[2*ai  ][2*aj+1] -= factor.imag();
			retVal[2*ai+1][2*aj  ] += factor.imag();
			retVal[2*ai+1][2*aj+1] += factor.real();
		}
	}
	return retVal;
}

bool integrator::writeToFile(std::string fileName) const {
	if (not _isIntegrated) {
		std::cerr << "integrator::writeToFile(...): ERROR: Not integrated yet. Returning empty vector." << std::endl;
		return false;
	}
	std::ofstream outFile;
	outFile.open (fileName.c_str());
	for (size_t i = 0; i < nAmpl(); ++i) {
		for (size_t j = 0; j < nAmpl(); ++j) {
			outFile << _integralMatrix[i][j] << " ";
		}
		outFile << std::endl;
	}
	outFile.close();
	return true;
}

bool integrator::setNpoints(size_t n) {
	_nPoints = n;
	return true;
}

std::pair<bool, std::string> integrator::getWaveName(size_t i) const {
	if (i >= nAmpl()) {
		return std::pair<bool, std::string>(false, "");
	}
	std::string waveName = _amplitudes[i]->name();
	return std::pair<bool, std::string>(true, waveName);
}

std::pair<bool, size_t> integrator::getNpoints() const {
	return std::pair<bool, size_t>(true, _nPoints);
}
