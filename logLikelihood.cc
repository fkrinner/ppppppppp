#include"logLikelihood.h"
#include<string>
#include<iostream>
#include<nlopt.hpp>

double ff_nlopt(const std::vector<double> &x, std::vector<double> &grad, void* f_data) {
	logLikelihood* inst = reinterpret_cast<logLikelihood*>(f_data);
	double ret   = inst->nloptCall(x,grad);
	return ret;
}

logLikelihood::logLikelihood(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral):
	_fixFirstPhase(false), _nCalls(0), _nCallsPrint(10000), _nAmpl(amplitudes.size()), _nPoints(0), _integral(integral), _amplitudes(amplitudes), _points() {
	if (_amplitudes.size() == 0) {
		std::cerr << "logLikelihood::logLikelihood(...): ERROR: No amplitudes given" << std::endl;
		throw;
	}	
	_kinSignature = _amplitudes[0]->kinSignature();
	if (not (_kinSignature == integral->kinSignature())) {
		std::cerr << "logLikelihood::logLikelihood(...): ERROR: Kinematic signature of integral differs" << std::endl;
		throw;
	}
	for (std::shared_ptr<amplitude> a : _amplitudes) {
		if (not(a->kinSignature() == _kinSignature)) {
			std::cerr << "logLikelihood::logLikelihood(...): ERROR: Amplitudes have different kinematic signatures" << std::endl;
			throw;
		}
	}
	if (_nAmpl != _integral->nAmpl()) {
		std::cerr  << "logLikelihood::logLikelihood(...): ERROR: Number of amplitudes does not match in integral" << std::endl;
	}
	for (size_t a = 0; a < _nAmpl; ++a) {
		std::pair<bool, std::string> waveName =  _integral->getWaveName(a);
		if (not waveName.first) {
			std::cerr  << "logLikelihood::logLikelihood(...): ERROR: Could not get wave name from integral" << std::endl;
			throw;
		}
		if (not (_amplitudes[a]->name() == waveName.second)) {
			std::cerr  << "logLikelihood::logLikelihood(...): ERROR: wave names do not match" << std::endl;
			throw;
		}
	}
	if (not _integral->isIntegrated()) {
		if (not _integral->integrate()) {
			std::cerr << "logLikelihood::logLikelihood(...): ERROR: Could not obtain integral" << std::endl;
			throw;
		}
	}
}

std::pair<double, std::vector<std::complex<double> > > logLikelihood::fit(std::vector<std::complex<double> >& parameters) {
	nlopt::opt opt = nlopt::opt(nlopt::LD_LBFGS, getNpar());
	opt.set_ftol_abs(1.E-2);
	opt.set_min_objective(ff_nlopt, this);
	std::vector<double> startPars(getNpar());
	size_t skip = 0;
	if (_fixFirstPhase) {
		startPars[0] = parameters[0].real();
		skip = 1;
	} else {
		startPars[0] = parameters[0].real();
		startPars[1] = parameters[0].imag();
		skip = 2;
	}
	for (size_t a = 1; a < _nAmpl; ++a) {
		startPars[2*a-2+skip] = parameters[a].real();
		startPars[2*a-1+skip] = parameters[a].imag();
	}
	double bestVal;
	opt.optimize(startPars, bestVal);
	std::vector<std::complex<double> >retVal(_nAmpl);
	if (_fixFirstPhase) { // skip is already set accordingly
		retVal[0] = std::complex<double>(startPars[0], 0.);
	} else {
		retVal[0] = std::complex<double>(startPars[0], startPars[1]);
	}
	for (size_t a = 1; a < _nAmpl; ++a) {
		retVal[a] = std::complex<double>(startPars[2*a-2+skip], startPars[2*a-1+skip]);
	}
//	std::cout << "The fit gave a best value of " << bestVal << std::endl;
	return std::pair<double, std::vector<std::complex<double> > >(bestVal, retVal);
}

double logLikelihood::eval(std::vector<std::complex<double> >& prodAmps) const {
	if (not prodAmps.size() == _nAmpl) {
		std::cerr << "logLikelihood::eval(...): ERROR: Number of production amplitudes does not match" << std::endl;
		throw; // Throw here, since a ll = 0. could confuse the program
	}
	double ll = 0.;
//#pragma omp parallel for reduction(+:ll)
	for (size_t p = 0; p < _nPoints; ++p) {
		std::complex<double> ampl (0.,0.);
		for (size_t a = 0; a < _nAmpl; ++a) {
			ampl += _points[p][a] * prodAmps[a];
		}
		ll += log(norm(ampl));
	}
	ll -= _integral->totalIntensity(prodAmps);
	++_nCalls;
	if (_nCalls%_nCallsPrint == 0) {
		std::cout << "call #" << _nCalls << " like " << -ll << std::endl;
	}	
	return -ll;
}

std::vector<double> logLikelihood::Deval(std::vector<std::complex<double> >& prodAmps) const {
	if (not prodAmps.size() == _nAmpl) {
		std::cerr << "logLikelihood::eval(...): ERROR: Number of production amplitudes does not match" << std::endl;
		throw; // Throw here, since a ll = 0. could confuse the program
	}
	std::vector<double> retVal(2*_nAmpl, 0.);
	for (size_t p = 0; p < _nPoints; ++p) {
		std::complex<double> ampl(0.,0.);
		for (size_t a = 0; a < _nAmpl; ++a) {
			ampl += _points[p][a] * prodAmps[a];
		}
		double intens = norm(ampl);
		ampl = std::conj(ampl); // Will only be used conjugated
		for (size_t a = 0; a < _nAmpl; ++a) {
			std::complex<double> factor = _points[p][a] * ampl/intens * 2.;
			retVal[2*a  ] -= factor.real(); // One factor of -1 since, the NEGATIVE likelihood is used 
			retVal[2*a+1] += factor.imag(); // Two factors of -1: Complex i*i = -1 and since the NEGATIVE likelihood is used
		}
	}
	std::vector<double> Dintegral = _integral->DtotalIntensity(prodAmps);
	for (size_t a = 0; a < 2*_nAmpl; ++a) {
		retVal[a] += Dintegral[a];
	}
	return retVal;
}

std::vector<std::vector<double> > logLikelihood::DDeval(std::vector<std::complex<double> >& prodAmps) const {
	if (not prodAmps.size() == _nAmpl) {
		std::cerr << "logLikelihood::eval(...): ERROR: Number of production amplitudes does not match" << std::endl;
		throw; // Throw here, since a ll = 0. could confuse the program
	}
	std::vector<std::vector<double> > retVal(2*_nAmpl, std::vector<double> (2*_nAmpl, 0.));
	for (size_t p = 0; p < _nPoints; ++p) {
		std::complex<double> ampl(0.,0.);
		for (size_t a = 0; a < _nAmpl; ++a) {
			ampl += _points[p][a] * prodAmps[a];
		}
		double intens = norm(ampl);
		ampl = std::conj(ampl); // Will only be used conjugated
		for(size_t ai = 0; ai < _nAmpl; ++ai) {
			std::complex<double> factori = _points[p][ai] * ampl/intens * 2.;
			for (size_t aj = 0; aj < _nAmpl; ++aj) {
				std::complex<double> factor = _points[p][ai] * std::conj(_points[p][aj])/intens*2.;
				retVal[2*ai  ][2*aj  ] -= factor.real();
				retVal[2*ai  ][2*aj+1] -= factor.imag();
				retVal[2*ai+1][2*aj  ] += factor.imag();
				retVal[2*ai+1][2*aj+1] -= factor.real();

				std::complex<double> factorj = -_points[p][aj] * ampl/intens *2.; // -1/intens and then the same factor;

				retVal[2*ai  ][2*aj  ] -= factori.real() * factorj.real();
				retVal[2*ai  ][2*aj+1] += factori.real() * factorj.imag();
				retVal[2*ai+1][2*aj  ] += factori.imag() * factorj.real();
				retVal[2*ai+1][2*aj+1] -= factori.imag() * factorj.imag();
			}
		}
	}
	std::vector<std::vector<double> > DDintegral = _integral->DDtotalIntensity(prodAmps);
	for (size_t i = 0; i < 2*_nAmpl; ++i) {
		for (size_t j = 0; j < 2*_nAmpl; ++j) {
			retVal[i][j] += DDintegral[i][j];
		}
	}
	return retVal;
}

double logLikelihood::nloptCall(const std::vector<double> &x, std::vector<double> &grad) const {
	if (not x.size() == getNpar()) {
		std::cerr << "logLikelihood::nloptCall(...): ERROR: Nomber of amplitudes does not match" << std::endl;
		throw;
	}
	std::vector<std::complex<double> > prodAmps(_nAmpl);
	size_t skip = 0;
	if (_fixFirstPhase) {
		prodAmps[0] = std::complex<double>(x[0], 0.);
		skip = 1;
	} else {
		prodAmps[0] = std::complex<double>(x[0], x[1]);
		skip = 2;
	}
	for (size_t a = 1; a < _nAmpl; ++a) {
		prodAmps[a] = std::complex<double>(x[2*a-2+skip],x[2*a-1+skip]);
	}
	if (grad.size() > 0) {
		std::vector<double> diff = Deval(prodAmps);
		grad[0] = diff[0];
		if (not _fixFirstPhase) { // Skip is already set accordingly
			grad[1] = diff[1];
		}
		for (size_t a = 1; a < _nAmpl; ++a) {
			grad[2*a-2+skip] = diff[2*a  ];
			grad[2*a-1+skip] = diff[2*a+1];
		}
	}
	return eval(prodAmps);
}

size_t logLikelihood::getNpar() const {
	size_t retVal = 2*_nAmpl;
	if (_fixFirstPhase) {
		retVal -= 1;
	}
	return retVal;
}

bool logLikelihood::loadDataPoints(const std::vector<std::vector<double> >& dataPoints) {
	if (dataPoints.size() == 0) {
		std::cerr << "logLikelihood::addDataPoints(...): ERROR: Not data points given" << std::endl;
		return false;
	}
	_nPoints = dataPoints.size();
	_points = std::vector<std::vector<std::complex<double> > > (_nPoints, std::vector<std::complex<double> > (_nAmpl));
	for(size_t a = 0; a < _nAmpl; ++a) {
		std::pair<bool, std::complex<double> > diag = _integral->element(a,a);
		if (not diag.first) {
			std::cerr << "logLikelihood::addDataPoints(...): ERROR: Could not get diagonal integral" << std::endl;
			return false;
		}
		double norm = 0.;
		if (not diag.second.real() == 0.) {
			norm = 1./pow(diag.second.real(), .5);
		}
		for (size_t p = 0; p < _nPoints; ++p) {
			_points[p][a] = _amplitudes[a]->eval(dataPoints[p])*norm;
		}
	}
	return true;
}

