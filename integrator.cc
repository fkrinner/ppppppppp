#include"integrator.h"
#include"utils.h"
#include<iostream>
#include<fstream>
#include<iomanip>
#include<limits>
integrator::integrator(size_t integralPoints, std::shared_ptr<generator> pointGenerator, const std::vector<std::shared_ptr<amplitude> >& amplitudes, std::shared_ptr<efficiencyFunction>& efficiency):
	_isIntegrated(false), _numLim(1.e-15), _nAmpl(amplitudes.size()), _nPoints(integralPoints), _amplitudes(amplitudes), _generator(pointGenerator), _efficiency(efficiency), _integralMatrix(), _accCorrIntegralMatrix(), _realIntegralMatrix(), _realAccCorrMatrix() {
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
	if (!utils::checkComplexDouble()) {
		std::cerr << "integrator::integrator(...): ERROR: Complex double has wrong structure" << std::endl;
		throw;
	}
	_amplitudeCoherenceBorders = {0, _nAmpl};
}

bool integrator::integrate() {
	std::vector<std::vector<std::complex<double> > > integralMatrix       (nAmpl(), std::vector<std::complex<double> > (nAmpl(), std::complex<double>(0.,0.)));
	std::vector<std::vector<std::complex<double> > > accCorrIntegralMatrix(nAmpl(), std::vector<std::complex<double> > (nAmpl(), std::complex<double>(0.,0.)));
	for (size_t point = 0; point < _nPoints; ++point) {
		std::vector<double> kin = _generator->generate();
		std::vector<std::complex<double> > ampl(nAmpl(), std::complex<double>(0.,0.));
		double eff = _efficiency->eval(kin);
		size_t countAmp = 0;
		for (std::shared_ptr<amplitude> a : _amplitudes) {
			ampl[countAmp] = a->eval(kin);
//			if (std::isnan(ampl[countAmp].real()) || std::isnan(ampl[countAmp].imag())) {
//				std::cout << "integrator::loadIntegrals(...): ERROR: NaN amplitude encountered for " << a->name() << std::endl;
//				return false;
//			}
			++countAmp;
		}
		for (size_t i = 0; i < nAmpl(); ++i) {
			for (size_t j = 0; j < nAmpl(); ++j) {
				integralMatrix[i][j] += std::conj(ampl[i])*ampl[j];
				accCorrIntegralMatrix[i][j] += std::conj(ampl[i])*ampl[j]*eff;
			}
		}
		if (point % 1000 == 0) { 
			std::cout << "#" << point << std::endl;
		}
	}
	for (size_t i = 0; i < nAmpl(); ++i) {
		for (size_t j = 0; j < nAmpl(); ++j) {
			integralMatrix[i][j]        /= _nPoints;
			accCorrIntegralMatrix[i][j] /= _nPoints;
		}
	}
	if (not setIntegrals(integralMatrix, accCorrIntegralMatrix)) {
		std::cout << "integrator::loadIntegrals(...): ERROR: Could not set integrated matrices" << std::endl;
		return false;
	}
	return true;
}

bool integrator::loadIntegrals(const std::string& psFileName, const std::string& accFileName) {
	// Trust the content of the files, as long as the dimensions match
	std::pair<bool, std::vector<std::vector<std::complex<double> > > > ps_load = utils::readComplexMatrixFromTextFile(psFileName, _amplitudes.size());
	if (!ps_load.first) {
		std::cout << "integrator::loadIntegrals(...): ERROR: Could not load ps matrix from '" << psFileName << "'" << std::endl;
		return false;
	}
	std::pair<bool, std::vector<std::vector<std::complex<double> > > > ac_load = utils::readComplexMatrixFromTextFile(accFileName, _amplitudes.size());
	if (!ac_load.first) {
		std::cout << "integrator::loadIntegrals(...): ERROR: Could not load ac matrix from '" << accFileName << "'" <<  std::endl;
		return false;
	}
	if (not setIntegrals(ps_load.second, ac_load.second)) {
		std::cout << "integrator::loadIntegrals(...): ERROR: Could not set loaded matrices" << std::endl;
		return false;
	}
	return true;
}

bool integrator::setIntegrals(const std::vector<std::vector<std::complex<double> > >& ps_integral, const std::vector<std::vector<std::complex<double> > >& ac_integral) {
	if (ps_integral.size() != _nAmpl) {
		std::cout << "integrator::setIntegrals(...): ERROR: ps_integral has the wrong size" << std::endl;
		return false;
	}
	if (ac_integral.size() != _nAmpl) {
		std::cout << "integrator::setIntegrals(...): ERROR: ac_integral has the wrong size" << std::endl;
		return false;
	}
	for (size_t a = 0; a < _nAmpl; ++a) {
		if (ps_integral[a].size() != _nAmpl) {
			std::cout << "integrator::setIntegrals(...): ERROR: ps_integral has the wrong size (line " << a << ")" << std::endl;
			return false;
		}	
		if (ac_integral[a].size() != _nAmpl) {
			std::cout << "integrator::setIntegrals(...): ERROR: ac_integral has the wrong size (line " << a << ")" << std::endl;
			return false;
		}	
	}
	for (size_t a_i = 0; a_i < _nAmpl; ++a_i) {
		if (ps_integral[a_i][a_i].imag() != 0.) {
			std::cout << "integrator::setIntegrals(...): ERROR: ps_integral has non-vanisihng imaginary part on the diagonal: " << ps_integral[a_i][a_i] << "(wave '" << _amplitudes[a_i]->name() << "')" << std::endl;
			return false;
		}
		if (ac_integral[a_i][a_i].imag() != 0.) {
			std::cout << "integrator::setIntegrals(...): ERROR: ac_integral has non-vanisihng imaginary part on the diagonal: " << ac_integral[a_i][a_i]  << "(wave '" << _amplitudes[a_i]->name() << "')" << std::endl;
			return false;
		}
		for (size_t a_j = 0; a_j < a_i; ++a_j) {
			if (norm(ps_integral[a_i][a_j] - std::conj(ps_integral[a_j][a_i])) > _numLim ) {
				std::cout << "integrator::setIntegrals(...): ERROR: ps_integral is not hermitian: " << ps_integral[a_i][a_j] << " vs. " << ps_integral[a_j][a_i] << std::endl;
				return false;
			}
			if (norm(ac_integral[a_i][a_j] - std::conj(ac_integral[a_j][a_i])) > _numLim ) {
				std::cout << "integrator::setIntegrals(...): ERROR: ac_integral is not hermitian: " << ac_integral[a_i][a_j] << " vs. " << ac_integral[a_j][a_i] << std::endl;
				return false;
			}
		}
	}
	bool nanAmpl = false;
	for (size_t a = 0; a < _nAmpl; ++a) {
		if (std::isnan(ac_integral[a][a].real()) || std::isnan(ac_integral[a][a].imag())) {
			nanAmpl = true;
			std::cout << "integrator::setIntegrals(...): ERROR: NaN in ac integral for wave '" << _amplitudes[a]->name() << "': " << ac_integral[a][a]  << std::endl;
		}
	}
       for (size_t a = 0; a < _nAmpl; ++a) {
                if (std::isnan(ps_integral[a][a].real()) || std::isnan(ps_integral[a][a].imag())) {
                        nanAmpl = true;
                        std::cout << "integrator::setIntegrals(...): ERROR: NaN in ps integral for wave '" << _amplitudes[a]->name() << "': " << ps_integral[a][a] << std::endl;
                }
        }

	if (nanAmpl) {
		std::cout <<  "integrator::setIntegrals(...): ERROR: NaN amplitudes encountered" << std::endl;
		return false;
	}
	_integralMatrix = ps_integral;
	_accCorrIntegralMatrix =  ac_integral;
	if (not makeRealMatrices()) {
		std::cout <<  "integrator::setIntegrals(...): ERROR: Could bit create real versions of the matrices" << std::endl;
		return false;
	}
	_isIntegrated = true;
	return true;
}

bool integrator::makeRealMatrices() {
	_realIntegralMatrix = std::vector<std::vector<double> >(2*_nAmpl, std::vector<double>(2*_nAmpl, 0.));
	_realAccCorrMatrix = std::vector<std::vector<double> >(2*_nAmpl, std::vector<double>(2*_nAmpl, 0.));
	for (size_t a_i = 0; a_i < _nAmpl; ++a_i) {
		for (size_t a_j = 0; a_j < _nAmpl; ++a_j) {
			double norm = pow((_integralMatrix[a_i][a_i] * _integralMatrix[a_j][a_j]).real(),.5);
			if (norm == 0.) {
				continue;
			}
			_realIntegralMatrix[2*a_i  ][2*a_j  ] = _integralMatrix[a_i][a_j].real()/norm;
			_realIntegralMatrix[2*a_i  ][2*a_j+1] =-_integralMatrix[a_i][a_j].imag()/norm;
			_realIntegralMatrix[2*a_i+1][2*a_j  ] = _integralMatrix[a_i][a_j].imag()/norm;
			_realIntegralMatrix[2*a_i+1][2*a_j+1] = _integralMatrix[a_i][a_j].real()/norm;

			_realAccCorrMatrix[2*a_i  ][2*a_j  ] = _accCorrIntegralMatrix[a_i][a_j].real()/norm;
			_realAccCorrMatrix[2*a_i  ][2*a_j+1] =-_accCorrIntegralMatrix[a_i][a_j].imag()/norm;
			_realAccCorrMatrix[2*a_i+1][2*a_j  ] = _accCorrIntegralMatrix[a_i][a_j].imag()/norm;
			_realAccCorrMatrix[2*a_i+1][2*a_j+1] = _accCorrIntegralMatrix[a_i][a_j].real()/norm;
		}
	}
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

std::pair<bool, std::vector<double> > integrator::getNormalizations(bool accCorr) const {
	if (!_isIntegrated) {
		std::cout << "integrator::getNormalizations(...): ERROR: Cannot get normalizations before integration" << std::endl;
		return std::pair<bool, std::vector<double> >(false, std::vector<double>());
	}
	std::vector<double> retVal(_nAmpl);
	std::complex<double> val;
	for (size_t a = 0; a < _nAmpl; ++a) {
		if (accCorr) {
			val = _accCorrIntegralMatrix[a][a];
		} else {
			val = _integralMatrix[a][a];
		}
		if (val.imag() != 0.) {
			std::cout << "integrator::getNormalizations(...): ERROR: non-vanishing imaginary part on the diagonal (wave: '" << _amplitudes[a]->name() << "')" << std::endl;
			return std::pair<bool, std::vector<double> >(false, std::vector<double>());
		}
		if (val.real() == 0.) {
			retVal[a] = 0.;	
		} else {
			retVal[a] = 1./pow(val.real(), .5);
		}
	}
	return std::pair<bool, std::vector<double> >(true, retVal);
}

double integrator::totalIntensity(const std::vector<std::complex<double> >& prodAmpl, bool accCorr) const {
	if (not _isIntegrated) {
		std::cerr << "integrator::totalIntensity(...): ERROR: Not integrated yet. Returning zero." << std::endl;
		return 0.;
	}
	if (prodAmpl.size() != nAmpl()) {
		std::cerr << "integrator::totalIntensity(...): ERROR: Number of prodiction amplitudes does not match. Returning zero." << std::endl;
		return 0.;
	}
//	const std::vector<std::vector<std::complex<double> > > &integralMatrix = accCorr ? _accCorrIntegralMatrix : _integralMatrix;
//	std::complex<double> retVal(0.,0.);
//	std::cout << "-------------------------------------" << std::endl;
//	for (size_t i = 0; i < nAmpl(); ++i) {
//		for (size_t j = 0; j < nAmpl(); ++j) {
//			double norm = pow(_integralMatrix[i][i].real() * _integralMatrix[j][j].real(), .5); // Always normalized to phaseSpace
//			if (norm == 0.) {
//				continue;
//			}
//			if (std::conj(prodAmpl[i]) * integralMatrix[i][j] * prodAmpl[j]/norm !=  0.) {
//				std::cout << i << " " << j << "resse "  <<  std::conj(prodAmpl[i]) << " " << integralMatrix[i][j]/norm << " " << prodAmpl[j] << std::endl;
//			}
//			retVal += std::conj(prodAmpl[i]) * integralMatrix[i][j] * prodAmpl[j]/norm;			
//		}
//	}
	double alternative = 0.;
	const double *reIm = (double*)&prodAmpl[0]; // Trust this at the moment
	const std::vector<std::vector<double> > &realMatrix = accCorr ? _realAccCorrMatrix : _realIntegralMatrix;
	for (size_t a_i = 0; a_i < 2*_nAmpl; ++a_i) {
		for (size_t a_j = 0; a_j < 2*_nAmpl; ++a_j) {
//			if ( reIm[a_i]*realMatrix[a_i][a_j]*reIm[a_j] != 0.) {
//				std::cout << a_i << " " << a_j << "innse " << reIm[a_i]<< " " <<realMatrix[a_i][a_j] << " " << reIm[a_j] << std::endl;
//			}
			alternative += reIm[a_i]*realMatrix[a_i][a_j]*reIm[a_j];
		}
	}
	return alternative;	
}

std::vector<double> integrator::DtotalIntensity(const std::vector<std::complex<double> >& prodAmpl, bool accCorr) const {
// Return a vector of length 2*nAmpl() for the derivative w.r.t. re, and im : {dRe0, dIm0, dRe1, ..., dImnAmpl()}
	if (not _isIntegrated) {
		std::cerr << "integrator::DtotalIntensity(...): ERROR: Not integrated yet. Returning empty vector." << std::endl;
		return std::vector<double>();
	}
	if (prodAmpl.size() != nAmpl()) {
		std::cerr << "integrator::DtotalIntensity(...): ERROR: Number of prodiction amplitudes does not match. Returning zero." << std::endl;
		return std::vector<double>();
	}
	std::vector<double> retVal(2*nAmpl(), 0.);
	const std::vector<std::vector<double> > &matrix = accCorr ? _realAccCorrMatrix : _realIntegralMatrix;
	const double* reIm=  (double*)&prodAmpl[0];
	for (size_t a_i = 0; a_i < 2*_nAmpl; ++a_i) {
		for (size_t a_j = 0; a_j < 2*_nAmpl; ++ a_j) {
			retVal[a_i] += (matrix[a_i][a_j] + matrix[a_j][a_i])*reIm[a_j];
		}
	}

//	const std::vector<std::vector<std::complex<double> > > &integralMatrix = accCorr ? _accCorrIntegralMatrix : _integralMatrix;
//	for (size_t j = 0; j < nAmpl(); ++j) {
//		for (size_t i = 0; i < nAmpl(); ++i) {
//			double norm = pow(_integralMatrix[i][i].real() * _integralMatrix[j][j].real(), .5); // Always normalized to phaseSpace
//			if (norm == 0.) {
//				continue;
//			}
//			std::complex<double> factor = integralMatrix[i][j] * prodAmpl[j]/norm*2.;
//			retVal[2*i  ] += factor.real();
//			retVal[2*i+1] += factor.imag();
//
//		}
//	}
	return retVal;
}

std::vector<std::vector<double> > integrator::DDtotalIntensity(const std::vector<std::complex<double> >& prodAmpl, bool accCorr) const {
	if (not _isIntegrated) {
		std::cerr << "integrator::DDtotalIntensity(...): ERROR: Not integrated yet. Returning empty vector." << std::endl;
		return std::vector<std::vector<double> >();
	}
	if (prodAmpl.size() != nAmpl()) {
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

bool integrator::writeToFile(std::string fileName, bool accCorr) const {
	if (not _isIntegrated) {
		std::cerr << "integrator::writeToFile(...): ERROR: Not integrated yet. Returning empty vector." << std::endl;
		return false;
	}
	std::ofstream outFile;
	outFile.open (fileName.c_str());
	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (size_t i = 0; i < nAmpl(); ++i) {
		for (size_t j = 0; j < nAmpl(); ++j) {
			if (accCorr) {
				outFile << _accCorrIntegralMatrix[i][j] << " ";
			} else {
				outFile << _integralMatrix[i][j] << " ";
			}
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

bool integrator::setCoherenceBorders(std::vector<size_t>& borders) {
	if (borders[0] != 0) {
		std::cout << "integrator::setCoherenceBorders(...): ERROR: First border has to be zero" << std::endl;
		return false;
	}
	_amplitudeCoherenceBorders = std::vector<size_t>(borders.size(),0);
	for (size_t i = 1; i < borders.size(); ++i) {
		if (borders[i] <= _amplitudeCoherenceBorders[i-1]) {
			std::cout << "integrator::setCoherenceBorders(...): ERROR: Borders are nor ordered" << std::endl;
			return false;
		}
		 _amplitudeCoherenceBorders[i] = borders[i];
	}
	if (_amplitudeCoherenceBorders[_amplitudeCoherenceBorders.size()] != _nAmpl) {
		std::cout << "integrator::setCoherenceBorders(...): ERROR: Last border has to be _nAmpl" << std::endl;
		return false;
	}
	return true;
}

bool integrator::addIncoherentSector(std::shared_ptr<integrator> sector) {
	if (not isIntegrated()) {
		std::cout << "integrator::addIncoherentSector(...): ERROR: Cannot add second sector before integration" << std::endl;
		return false;
	}
	if (not sector->isIntegrated()) {
		std::cout << "integrator::addIncoherentSector(...): ERROR: Cannot add non-integrated sector" << std::endl;
		return false;
	}
	if (not (*_kinSignature == *(sector->kinSignature()))) {
		std::cout << "integrator::addIncoherentSector(...): ERROR: Kinematic signatures have to match" << std::endl;
		return false;
	}
	if (_nPoints != sector->getNpoints().second) {
		std::cout << "integrator::addIncoherentSector(...): ERROR: Number of points has to match" << std::endl;
		return false;
	}
	
	size_t nNew   = sector->nAmpl();
	size_t nTotal = _nAmpl + nNew;
	std::vector<std::vector<std::complex<double> > > new_ps_integral(nTotal, std::vector<std::complex<double> >(nTotal, std::complex<double>(0.,0.)));
	std::vector<std::vector<std::complex<double> > > new_ac_integral(nTotal, std::vector<std::complex<double> >(nTotal, std::complex<double>(0.,0.)));
	for (size_t i = 0; i < _nAmpl; ++i) {
		for (size_t j = 0; j < _nAmpl; ++j) {
			new_ps_integral[i][j] = _integralMatrix[i][j];
			new_ac_integral[i][j] = _accCorrIntegralMatrix[i][j];
		}	
	}
	for (size_t i = 0; i < nNew; ++i) {
		for (size_t j = 0; j < nNew; ++j) {
			new_ps_integral[_nAmpl+i][_nAmpl+j] = sector->element(i,j,false).second;
			new_ac_integral[_nAmpl+i][_nAmpl+j] = sector->element(i,j,true).second;
		}
	}
	for (size_t i : sector->getCoherenceBorders()) {
		if (i > 0) { // Do not add the leading 0
			_amplitudeCoherenceBorders.push_back(_nAmpl + i);
		}
	}
	_nAmpl = nTotal;
	if (not setIntegrals(new_ps_integral, new_ac_integral)) {
		std::cout << "integrator::addIncoherentSector(...): ERROR: Could not set new integral matrices" << std::endl;
		return false;
	}
	std::vector<std::shared_ptr<amplitude> > newAmplitudes;
	for (std::shared_ptr<amplitude> a : _amplitudes) {
		newAmplitudes.push_back(a);
	}
	for (std::shared_ptr<amplitude> a : (*sector)._amplitudes) {
		newAmplitudes.push_back(a);
	}
	_amplitudes = newAmplitudes;
	return true;	
}
