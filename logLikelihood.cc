#include"logLikelihood.h"
#include<string>
#include<nlopt.hpp>


double ff_nlopt(const std::vector<double> &x, std::vector<double> &grad, void* f_data) {
	logLikelihood* inst = reinterpret_cast<logLikelihood*>(f_data);
	double ret   = inst->nloptCall(x,grad);
	return ret;
}

logLikelihoodBase::logLikelihoodBase (size_t nAmpl, std::shared_ptr<kinematicSignature> kinSignature, std::shared_ptr<integrator> integral) : 
	_fixFirstPhase(false), _extended(true), _kinSignature(kinSignature), _nCalls(0), _nCallsPrint(10000), _nAmpl(nAmpl), _nPoints(0), _integral(integral) {
	if (_nAmpl == 0) {
		std::cout << "logLikelihoodBase::logLikelihoodBase(...): ERROR: No amplitudes given" << std::endl;
		throw;
	}
	if (_nAmpl != _integral->nAmpl()) {
		std::cout << "logLikelihoodBase::logLikelihoodBase(...): ERROR: Number of amplitudes does not match" << std::endl;	
		throw;
	}
	if (not (*(_integral->kinSignature()) == *_kinSignature)) {
		std::cout << "logLikelihoodBase::logLikelihoodBase(...): ERROR: Kinematic signatures doe not match" << std::endl;
		throw;	
	}
	if (not _integral->isIntegrated()) {
		if (not _integral->integrate()) {
			std::cerr << "logLikelihoodBase::logLikelihoodBase(...): ERROR: Could not obtain integral" << std::endl;
			throw;
		}
	}
}

std::pair<double, std::vector<std::complex<double> > > logLikelihoodBase::fitNlopt(std::vector<std::complex<double> >& parameters) {
	nlopt::opt opt = nlopt::opt(nlopt::LD_LBFGS, getNpar());
	opt.set_ftol_abs(1.E-14);
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

double logLikelihoodBase::nloptCall(const std::vector<double> &x, std::vector<double> &grad) const {
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

size_t logLikelihoodBase::getNpar() const {
	size_t retVal = 2*_nAmpl;
	if (_fixFirstPhase) {
		retVal -= 1;
	}
	return retVal;
}

/*
std::pair<double, std::vector<std::complex<double> > > logLikelihood::fitROOT(std::vector<std::complex<double> >& parameters) {
	std::unique_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
	min->SetMaxFunctionCalls(100000);
	min->SetMaxIterations(100000);
	min->SetTolerance(1.);
	ROOT::Math::Functor functor(*this,getNpar());
	min->SetFunction(functor);

	size_t skip = 1;
	std::complex<double> avg(0.,0.);
	for (std::complex<double> a : parameters) {
		avg += a;
	}
	avg /= nAmpl();
	double step = pow(std::norm(avg), .5)/10.;
	min->SetVariable(0, "real0", parameters[0].real(), step);
	if (not _fixFirstPhase) {
		++skip;
		min->SetVariable(1, "imag0", parameters[1].imag(), step);
	}
	for (size_t a = 1; a < nAmpl(); ++a) {
		std::string realName = std::string("real") + std::to_string(a);
		std::string imagName = std::string("imag") + std::to_string(a);
		min->SetVariable(2*a-2+skip, realName, parameters[a].real(), step);
		min->SetVariable(2*a-1+skip, imagName, parameters[a].imag(), step);
	}
	min->Minimize();
	const double *xs = min->X();
	std::vector<std::complex<double> > bestPar(nAmpl());
	if (_fixFirstPhase) {
		bestPar[0] = std::complex<double>(xs[0], 0.);
	} else {
		bestPar[0] = std::complex<double>(xs[0], xs[1]);
	}
	for (size_t a = 1; a < nAmpl(); ++a) {
		bestPar[a] = std::complex<double>(xs[2*a-2+skip], xs[2*a-1+skip]);
	}
	double retVal = eval(bestPar);

//	std::cout << "Calculating hessian" << std::endl;
//	min->Hesse();
//	std::cout << "Hessian calculated" << std::endl;
//	double* hessianValues = 0;
//	if (not min->GetHessianMatrix(hessianValues)) {
//		std::cout << "logLikelihood::fitROOT(...): ERROR: Could not get hessian values" << std::endl;
//		throw;
//	}
//	_hessian = std::vector<std::vector<double> >(2*nAmpl(), std::vector<double>(2*nAmpl(), 0.));
//	size_t nLineHessian = 2*_nAmpl;
//	if (_fixFirstPhase) {
//		nLineHessian -= 1;
//	}
//	_hessian[0][0] = hessianValues[0];
//	if (not _fixFirstPhase) {
//		_hessian[0][1] = hessianValues[1];
//		_hessian[1][0] = hessianValues[nLineHessian];
//		_hessian[1][1] = hessianValues[nLineHessian+1];
//	}
//	for (size_t a1 = 1; a1 < nAmpl(); ++a1) {
//		_hessian[0][2*a1  ] = hessianValues[2*a1  ];
//		_hessian[0][2*a1+1] = hessianValues[2*a1+1];
//
//		_hessian[2*a1  ][0] = hessianValues[(2*a1  )*nLineHessian];
//		_hessian[2*a1+1][0] = hessianValues[(2*a1+1)*nLineHessian];
//		if (not _fixFirstPhase) {
//			_hessian[1][2*a1  ] = hessianValues[2*a1   + nLineHessian];
//			_hessian[1][2*a1+1] = hessianValues[2*a1+1 + nLineHessian];
//
//			_hessian[2*a1  ][1] = hessianValues[1+(2*a1  )*nLineHessian];
//			_hessian[2*a1+1][1] = hessianValues[1+(2*a1+1)*nLineHessian];
//		}
//		for (size_t a2 = 1; a2 < nAmpl(); ++a2) {
//			std::cout << a1 << " " << a2 << " " << nAmpl() << std::endl;
//			_hessian[2*a1  ][2*a2  ] = hessianValues[2*a1  +(2*a2  )*nLineHessian];
//			_hessian[2*a1  ][2*a2+1] = hessianValues[2*a1  +(2*a2+1)*nLineHessian];
//			_hessian[2*a1+1][2*a2  ] = hessianValues[2*a1+1+(2*a2  )*nLineHessian];
//			_hessian[2*a1+1][2*a2+1] = hessianValues[2*a1+1+(2*a2+1)*nLineHessian];
//		}
//	}

	return std::pair<double, std::vector<std::complex<double> > >(retVal, bestPar);
}

double logLikelihood::operator()(const double* parameters) const {
	std::vector<std::complex<double> > prodAmps(nAmpl());
	size_t skip = 0;
	if (_fixFirstPhase) {
		skip = 1;
		prodAmps[0] = std::complex<double>(parameters[0], 0.);

	} else {
		skip = 2;
		prodAmps[0] = std::complex<double>(parameters[0], parameters[1]);
	}
	for (size_t a = 1; a < nAmpl(); ++a){
		prodAmps[a] = std::complex<double>(parameters[2*a-2+skip],parameters[2*a-1+skip]);
	}
	return eval(prodAmps);
}
*/

logLikelihood::logLikelihood(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral) :
	logLikelihoodBase(amplitudes.size(), integral->kinSignature(), integral), _nSect(1), _amplitudes(amplitudes), _points(), _amplitudeCoherenceBorders(0), _contributingWaves(0) {
	if (_amplitudes.size() == 0) {
		std::cerr << "logLikelihood::logLikelihood(...): ERROR: No amplitudes given" << std::endl;
		throw;
	}
	_kinSignature = _amplitudes[0]->kinSignature();
	if (not (*_kinSignature == *(integral->kinSignature()))) {
		std::cerr << "logLikelihood::logLikelihood(...): ERROR: Kinematic signature of integral differs" << std::endl;
		throw;
	}
	for (std::shared_ptr<amplitude> a : _amplitudes) {
		if (not(*(a->kinSignature()) == *_kinSignature)) {
			std::cerr << "logLikelihood::logLikelihood(...): ERROR: Amplitudes have different kinematic signatures" << std::endl;
			throw;
		}
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
	if (not setCoherenceBorders(_integral->getCoherenceBorders())) {
		std::cerr << "logLikelihood::logLikelihood(...): ERROR: Could not set coherence borders" << std::endl;
		throw;
	}
}

double logLikelihood::eval(std::vector<std::complex<double> >& prodAmps) const {
	if (prodAmps.size() != _nAmpl) {
		std::cerr << "logLikelihood::eval(...): ERROR: Number of production amplitudes does not match" << std::endl;
		throw; // Throw here, since a ll = 0. could confuse the program
	}
	double ll = 0.;
//#pragma omp parallel for reduction(+:ll)
	for (size_t p = 0; p < _nPoints; ++p) {
		size_t upperSectorBorder = _amplitudeCoherenceBorders[1]; // {0, b1, b2, ..., _nAmpl}.size() = nSector +1
		double intens = 0.;
		std::complex<double> ampl (0.,0.);
		for (size_t a : _contributingWaves[p]) {
			if (a >= upperSectorBorder) {
				intens += norm(ampl);
				ampl = std::complex<double>(0.,0.);
				upperSectorBorder = _amplitudeCoherenceBorders[getSector(a)+1];
			}
			ampl += _points[p][a] * prodAmps[a];
		}
		intens += norm(ampl); // Do also for the last sector
		ll += log(intens);
	}
	if (_extended) {
		ll -= _integral->totalIntensity(prodAmps, true);
	} else {
		ll -= log(_integral->totalIntensity(prodAmps, true))*(double)_nPoints;
	}
	++_nCalls;
	if (_nCalls%_nCallsPrint == 0) {
		std::cout << "call #" << _nCalls << " like " << -ll << std::endl;
	}
	return -ll;
}

std::vector<double> logLikelihood::Deval(std::vector<std::complex<double> >& prodAmps) const {
	if (prodAmps.size() != _nAmpl) {
		std::cerr << "logLikelihood::Deval(...): ERROR: Number of production amplitudes does not match" << std::endl;
		throw; // Throw here, since a ll = 0. could confuse the program
	}
	std::vector<double> retVal(2*_nAmpl, 0.);
	for (size_t p = 0; p < _nPoints; ++p) {
		size_t sector = 0;
		size_t upperSectorBorder = _amplitudeCoherenceBorders[1];
		double intens = 0.;
		std::vector<std::complex<double> > sectorAmpls(_nSect, std::complex<double>(0.,0.));
		for (size_t a : _contributingWaves[p]) {
			if (a >= upperSectorBorder) {
				intens += norm(sectorAmpls[sector]);
				sectorAmpls[sector]= std::conj(sectorAmpls[sector]);
				sector = getSector(a);
				upperSectorBorder = _amplitudeCoherenceBorders[sector+1];
			}
			sectorAmpls[sector] += _points[p][a] * prodAmps[a];
		}
		intens += norm(sectorAmpls[sector]); // Do also for the last sector
		sectorAmpls[sector]= std::conj(sectorAmpls[sector]);
//		ampl = std::conj(ampl); // Will only be used conjugated
		sector = 0;
		upperSectorBorder = _amplitudeCoherenceBorders[1];
		for (size_t a : _contributingWaves[p]) {
			if (a >= upperSectorBorder) {
				sector = getSector(a);
				upperSectorBorder = _amplitudeCoherenceBorders[sector+1];
			}
			std::complex<double> factor = _points[p][a] * sectorAmpls[sector]/intens * 2.;
			retVal[2*a  ] -= factor.real(); // One factor of -1 since, the NEGATIVE likelihood is used 
			retVal[2*a+1] += factor.imag(); // Two factors of -1: Complex i*i = -1 and since the NEGATIVE likelihood is used
		}
	}
	std::vector<double> Dintegral = _integral->DtotalIntensity(prodAmps, true);
	if (_extended) {
		for (size_t a = 0; a < 2*_nAmpl; ++a) {
			retVal[a] += Dintegral[a];
		}
	} else {
		double totalIntens = _integral->totalIntensity(prodAmps, true);
		for (size_t a = 0; a< 2*_nAmpl; ++a) {
			retVal[a] += Dintegral[a]/totalIntens*(double)_nPoints;
		}
	}
	return retVal;
}

std::vector<std::vector<double> > logLikelihood::DDeval(std::vector<std::complex<double> >& prodAmps) const {
	if (prodAmps.size() != _nAmpl) {
		std::cerr << "logLikelihood::DDeval(...): ERROR: Number of production amplitudes does not match" << std::endl;
		throw; // Throw here, since a ll = 0. could confuse the program
	}
	std::vector<std::vector<double> > retVal(2*_nAmpl, std::vector<double> (2*_nAmpl, 0.));
	for (size_t p = 0; p < _nPoints; ++p) {
		size_t sector = 0;
		size_t upperSectorBorder = _amplitudeCoherenceBorders[1];
		double intens = 0.;
		std::vector<std::complex<double> > sectorAmpls(_nSect, std::complex<double>(0.,0.));
		for (size_t a : _contributingWaves[p]) {
			if (a >= upperSectorBorder) {
				intens += norm(sectorAmpls[sector]);
				sectorAmpls[sector]= std::conj(sectorAmpls[sector]);
				sector = getSector(a);
				upperSectorBorder = _amplitudeCoherenceBorders[sector+1];
			}
			sectorAmpls[sector] += _points[p][a] * prodAmps[a];
		}
		intens += norm(sectorAmpls[sector]); // Do also for the last sector
		sectorAmpls[sector]= std::conj(sectorAmpls[sector]);

		size_t sector_i = 0;
		size_t upperSectorBorder_i = _amplitudeCoherenceBorders[1];
		for(size_t ai : _contributingWaves[p]) {
			if (ai >= upperSectorBorder_i) {
				sector_i = getSector(ai);
				upperSectorBorder_i = _amplitudeCoherenceBorders[sector_i+1];
			}
			std::complex<double> factori = _points[p][ai] * sectorAmpls[sector_i]/intens * 2.;

			size_t sector_j = 0;
			size_t upperSectorBorder_j = _amplitudeCoherenceBorders[1];
			for (size_t aj : _contributingWaves[p]) {
				if (aj >= upperSectorBorder_j) {
					sector_j = getSector(aj);
					upperSectorBorder_j = _amplitudeCoherenceBorders[sector_j+1];
				}
				if (sector_i == sector_j) {
					std::complex<double> factor = _points[p][ai] * std::conj(_points[p][aj])/intens*2.;
					retVal[2*ai  ][2*aj  ] -= factor.real();
					retVal[2*ai  ][2*aj+1] -= factor.imag();
					retVal[2*ai+1][2*aj  ] += factor.imag();
					retVal[2*ai+1][2*aj+1] -= factor.real();
				}
				std::complex<double> factorj = -_points[p][aj] * sectorAmpls[sector_j]/intens *2.; // -1/intens and then the same factor;

				retVal[2*ai  ][2*aj  ] -= factori.real() * factorj.real();
				retVal[2*ai  ][2*aj+1] += factori.real() * factorj.imag();
				retVal[2*ai+1][2*aj  ] += factori.imag() * factorj.real();
				retVal[2*ai+1][2*aj+1] -= factori.imag() * factorj.imag();
			}
		}
	}
	std::vector<std::vector<double> > DDintegral = _integral->DDtotalIntensity(prodAmps, true);
	if (_extended) {
		for (size_t i = 0; i < 2*_nAmpl; ++i) {
			for (size_t j = 0; j < 2*_nAmpl; ++j) {
				retVal[i][j] += DDintegral[i][j];
			}
		}
	} else {
		double totalIntens = _integral->totalIntensity(prodAmps, true);
		std::vector<double> Dintegral = _integral->DtotalIntensity(prodAmps, true);
		for (size_t i = 0; i < 2*_nAmpl; ++i) {
			for (size_t j = 0; j < 2*_nAmpl; ++j) {
				retVal[i][j] += (DDintegral[i][j]/totalIntens - Dintegral[i]*Dintegral[j]/totalIntens/totalIntens)*(double)_nPoints;
			}
		}
	}
	return retVal;
}

bool logLikelihood::loadDataPoints(const std::vector<std::vector<double> >& dataPoints) {
	if (dataPoints.size() == 0) {
		std::cerr << "logLikelihood::loadDataPoints(...): ERROR: Not data points given" << std::endl;
		return false;
	}
	_nPoints = dataPoints.size();
	_points  = std::vector<std::vector<std::complex<double> > > (_nPoints, std::vector<std::complex<double> > (_nAmpl));
	_contributingWaves= std::vector<std::vector<size_t> >(_nPoints, std::vector<size_t>(_nAmpl));
	std::vector<size_t> counts(_nPoints, 0);
	for(size_t a = 0; a < _nAmpl; ++a) {
		std::pair<bool, std::complex<double> > diag = _integral->element(a,a);
		if (not diag.first) {
			std::cerr << "logLikelihood::loadDataPoints(...): ERROR: Could not get diagonal integral" << std::endl;
			return false;
		}
		double norm = 0.;
		if (not diag.second.real() == 0.) {
			norm = 1./pow(diag.second.real(), .5);
		}
		for (size_t p = 0; p < _nPoints; ++p) {
			_points[p][a] = _amplitudes[a]->eval(dataPoints[p])*norm;
			if (_points[p][a] != std::complex<double>(0.,0.)) {
				_contributingWaves[p][counts[p]] = a;
				++counts[p];
			}
		}
	}
	for (size_t p = 0; p < _nPoints; ++p) {
		_contributingWaves[p].resize(counts[p]);
	}
	return true;
}

bool logLikelihood::setCoherenceBorders(std::vector<size_t> borders) {
	if (borders[0] != 0) {
		std::cout << "logLikelihood::setCoherenceBorders(...): ERROR: First border has to be zero" << std::endl;
		return false;
	}
	_amplitudeCoherenceBorders = std::vector<size_t>(borders.size(),0);
	for (size_t i = 1; i < borders.size(); ++i) {
		if (borders[i] <= _amplitudeCoherenceBorders[i-1]) {
			std::cout << "logLikelihood::setCoherenceBorders(...): ERROR: Borders are nor ordered" << std::endl;
			return false;
		}
		 _amplitudeCoherenceBorders[i] = borders[i];
	}
	if (_amplitudeCoherenceBorders[_amplitudeCoherenceBorders.size()-1] != _nAmpl) {
		std::cout << "logLikelihood::setCoherenceBorders(...): ERROR: Last border has to be _nAmpl ( = " << _nAmpl << " != "  << _amplitudeCoherenceBorders[_amplitudeCoherenceBorders.size()-1] << " )" << std::endl;
		return false;
	}
	_nSect = _amplitudeCoherenceBorders.size()-1;
	return true;
}

size_t logLikelihood::getSector(size_t a) const {
	for (size_t s = 0; s < _nSect; ++s) {
		if (_amplitudeCoherenceBorders[s] <= a && _amplitudeCoherenceBorders[s+1] > a) {
			return s;
		}
	}
	std::cout << "logLikelihood::getSector(...): ERROR: No sector found for amplitude index " << a << " (_nAmpl = " << _nAmpl << ")" << std::endl;
	throw;
	return 0;
}

logLikelihoodAllFree::logLikelihoodAllFree (std::vector<double> binning, std::vector<std::shared_ptr<angularDependence> > freedAmplitudes, std::shared_ptr<integrator> integral) :
	logLikelihoodBase ((binning.size()-1)*freedAmplitudes.size(), integral->kinSignature(), integral), 
	_nBins(binning.size()-1), 
	_binning(binning), 
	_eventsPerBin(std::vector<std::vector<size_t> >(binning.size()-1, std::vector<size_t>(binning.size() - 1, 0))),
	_amplitudesInBin(std::vector<std::vector<std::vector<std::complex<double> > > >(binning.size()-1, 
	                 std::vector<std::vector<std::complex<double> > >(binning.size()-1, 
	                 std::vector<std::complex<double> >(freedAmplitudes.size(), 
	                 std::complex<double>(0.,0.))))) {
	for (size_t b = 0; b < _nBins; ++b) {
		if (_binning[b+1] <= _binning[b]) {
			std::cout << "logLikelihoodAllFree::logLikelihoodAllFree(...): ERROR: Binning not ordered" << std::endl;
			throw;
		}
	}
	if (not _kinSignature->nKin() == 3) {
		std::cout << "Number of kinemtaic variables is not three... unknown prcess... Abortring..." << std::endl;
		throw;
	}
}

double logLikelihoodAllFree::eval(std::vector<std::complex<double> >& prodAmps) const {
	std::cout << "logLikelihoodAllFree::eval(...): NOT IMPLEMENTED " << prodAmps.size() << std::endl;
	throw;
	return 0;
}

std::vector<double> logLikelihoodAllFree::Deval(std::vector<std::complex<double> >& prodAmps) const {
	std::cout << "logLikelihoodAllFree::Deval(...): NOT IMPLEMENTED " << prodAmps.size() << std::endl;
	throw;
	std::vector<double>();
}

std::vector<std::vector<double> >  logLikelihoodAllFree::DDeval(std::vector<std::complex<double> >& prodAmps) const {
	std::cout << "logLikelihoodAllFree::DDeval(...): NOT IMPLEMENTED " << prodAmps.size() << std::endl;
	throw;
	std::vector<std::vector<double> >();
}

bool logLikelihoodAllFree::loadDataPoints(const std::vector<std::vector<double> >& dataPoints) {
	if (dataPoints.size() == 0) {
		std::cerr << "logLikelihood::loadDataPoints(...): ERROR: Not data points given" << std::endl;
		return false;
	}
	const double s = dataPoints[0][0];
	_nPoints = dataPoints.size();
	std::vector<std::vector<std::pair<double, double> > > averageMasses(_nBins, std::vector<std::pair<double, double> >(_nBins, std::pair<double, double>(0.,0.)));
	for (const std::vector<double>& point : dataPoints) {
		std::pair<bool, std::pair<size_t, size_t> > bin = findBin(point);
		if (not bin.first) {
			std::cout << "logLikelihood::loadDataPoints(...): ERROR: Bin not found" << std::endl;
			return false;
		}
		_eventsPerBin[bin.second.first][bin.second.second] += 1;
		averageMasses[bin.second.first][bin.second.second].first  += point[1];
		averageMasses[bin.second.first][bin.second.second].second += point[2];
	}
	for (size_t b1 = 0; b1 < _nBins; ++b1) {
		for (size_t b2 = 0; b2 < _nBins; ++b2) {
			size_t nBin = _eventsPerBin[b1][b2];
			if (nBin == 0) {
				continue;
			}
			std::vector<double> binPoint = {s, averageMasses[b1][b2].first/nBin, averageMasses[b1][b2].second/nBin};
			std::cout << "logLikelihood::loadDataPoints(...): NOT IMPLEMENTED: Not fully implemented" << std::endl;
			throw;
		}
	}
	return true;
}

std::pair<bool, std::pair<size_t, size_t> > logLikelihoodAllFree::findBin(const std::vector<double>& point) const {
	bool found1 = false;
	size_t bin1 = 0;
	bool found2 = false;
	size_t bin2 = 0;
	if (not point.size() == _kinSignature->nKin()) {
		std::cout << "logLikelihoodAllFree::findBin(...): Number of kinematic variables does not match" << std::endl;
		return std::pair<bool, std::pair<size_t, size_t> >(false, std::pair<size_t, size_t>(0,0));
	}
	const double s12 = point[1];
	const double s13 = point[2];
	for (size_t b = 0; b < _nBins; ++b) {
		const double sMin = _binning[b];
		const double sMax = _binning[b+1];
		if (sMin < s12 and s12 <= sMax) {
			if (found1) {
				std::cout << "Two bins found for s_{12}... Returning false" << std::endl;
				return  std::pair<bool, std::pair<size_t, size_t> >(false, std::pair<size_t, size_t>(0,0));
			}
			found1 = true;
			bin1   = b;
		}
		if (sMin < s13 and s13 <= sMax) {
			if (found2) {
				std::cout << "Two bins found for s_{13}... Returning false" << std::endl;
				return  std::pair<bool, std::pair<size_t, size_t> >(false, std::pair<size_t, size_t>(0,0));
			}
			found2 = true;
			bin2   = b;
		}
		if (found1 and found2) { // This breaks the total check of the binnig, i.e. if a secod bin would be found, but ordering is checked in the constuctor anyway...
			break;
		}
	}
	return std::pair<bool, std::pair<size_t, size_t> >(found1 and found2, std::pair<size_t, size_t>(bin1, bin2));
}
