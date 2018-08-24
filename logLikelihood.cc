#include"logLikelihood.h"

#include<string>
#include<limits>
#include<cmath>        // std::abs

#include<nlopt.hpp>
//#include<fstream>

double ff_nlopt(const std::vector<double> &x, std::vector<double> &grad, void* f_data) {
	logLikelihood* inst = reinterpret_cast<logLikelihood*>(f_data);
	double ret   = inst->nloptCall(x,grad);
	return ret;
}

logLikelihoodBase::logLikelihoodBase (size_t nAmpl, std::shared_ptr<kinematicSignature> kinSignature, std::shared_ptr<integrator> integral) : 
	_extended(true), _kinSignature(kinSignature), _nCalls(0), _nCallsPrint(10000), _nAmpl(nAmpl), _nPoints(0), _nScale(0), _integral(integral), _fixedParameters(), _copyParameters() {
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
//	opt.set_ftol_abs(1.E-14);
	opt.set_min_objective(ff_nlopt, this);
	std::vector<double> startPars = getFinalParams(prodAmpsToFullParams(parameters));
	double bestVal;
	int optimization_result  = opt.optimize(startPars, bestVal);
	std::cout << "optimization status: " <<  optimization_result << std::endl;
	if (optimization_result == 1) {
		std::cout << " - nlopt status: NLOPT_SUCCESS" << std::endl;
	} else if ( optimization_result == 2 ) {
		std::cout << " - nlopt status: NLOPT_STOPVAL_REACHED" << std::endl;
	} else if ( optimization_result == 3 ) {
		std::cout << " - nlopt status: NLOPT_FTOL_REACHED" << std::endl;
	} else if ( optimization_result == 4 ) {
		std::cout << " - nlopt status: NLOPT_XTOL_REACHED" << std::endl;
	} else if ( optimization_result == 5 ) {
		std::cout << " - nlopt status: NLOPT_MAXEVAL_REACHED" << std::endl;
	} else if ( optimization_result == 6 ) {
		std::cout << " - nlopt status: NLOPT_MAXTIME_REACHED" << std::endl;
	} else if ( optimization_result == -1 ) {
		std::cout << " - nlopt status: NLOPT_FAILURE" << std::endl;
	} else if ( optimization_result == -2 ) {
		std::cout << " - nlopt status: NLOPT_INVALID_ARGS" << std::endl;
	} else if ( optimization_result == -3 ) {
		std::cout << " - nlopt status: NLOPT_OUT_OF_MEMORY" << std::endl;
	} else if ( optimization_result == -4 ) {
		std::cout << " - nlopt status: NLOPT_ROUNDOFF_LIMITED" << std::endl;
	} else if ( optimization_result == -5 ) {
		std::cout << " - nlopt status: NLOPT_FORCED_STOP" << std::endl;
	} else {
		std::cout << " - Unknown nlopt status: " << optimization_result << std::endl;
	}

	std::vector<std::complex<double> > retVal = fullParamsToProdAmps(getFullParameters(startPars));
//	std::cout << "The fit gave a best value of " << bestVal << std::endl;
	return std::pair<double, std::vector<std::complex<double> > >(bestVal, retVal);
}

double logLikelihoodBase::nloptCall(const std::vector<double> &x, std::vector<double> &grad) const {
	if (not x.size() == getNpar()) {
		std::cerr << "logLikelihood::nloptCall(...): ERROR: Nomber of amplitudes does not match" << std::endl;
		throw;
	}
	std::vector<std::complex<double> > prodAmps = fullParamsToProdAmps(getFullParameters(x));
	if (grad.size() > 0) {
		grad = makeFinalGradient(x,Deval(prodAmps));
	}
	return eval(prodAmps);
}

std::vector<double> logLikelihoodBase::makeFinalGradient(const std::vector<double>& params, const std::vector<double>& fullGradient) const {

	const size_t nPar         = getNpar();
	const size_t nParTot      = getNparTot();
	const size_t nParNonScale = nParTot - _fixedParameters.size();

	if (params.size() != nPar) {
		std::cout << "logLikelihoodBase::makeFinalGradient(...): ERROR: Wrong parameter size." << std::endl;
		throw;
	}

	std::vector<double> fullParams;
	if (_nScale > 0) {
		fullParams = getFullParameters(params);
	}

	std::vector<double> scaleGradients(_nScale, 0.);
	std::vector<double> modifiedGradient = fullGradient;
	for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) { // First: modify the entries of the gradient corresponding to the actual parameters & collect the gradient for the scale parameters
		int scl = std::get<2>(copy);
		double scaler = std::get<3>(copy);
		if (scl > -1) {
			scaleGradients[scl] += fullGradient[std::get<0>(copy)] * fullParams[std::get<1>(copy)]*scaler;
			scaler *= params[nParNonScale + scl];
		}
		modifiedGradient[std::get<1>(copy)] += fullGradient[std::get<0>(copy)] * scaler;
	}

	modifiedGradient = getFinalParams(modifiedGradient); // Second: Cut away the fixed parameters from the grandient (this leaves the gradient entries for the scale parameters with odd values)
	for (size_t s = 0; s < _nScale; ++s) {
		modifiedGradient[nParNonScale+s] = scaleGradients[s]; // Third: Fix the grandient entries for the scale parameters 
	}
	return modifiedGradient;

}

std::vector<std::vector<double> > logLikelihoodBase::makeFinalHessian(const std::vector<double>& params, const std::vector<std::vector<double> >& fullHessian) const {

	const size_t nPar = getNpar();
	const size_t nParTot = getNparTot();
	const size_t nParNonScale = nParTot - _fixedParameters.size();

	std::vector<double> fullParams;
	if (_nScale > 0) {
		fullParams = getFullParameters(params);
	}
	std::vector<std::vector<double> > modifiedHessian(fullHessian.size());
// Do eveything for lines first (interprete the hessian as a list of gradients)
	for (size_t p = 0; p < fullHessian.size(); ++p) {
		std::vector<double> scaleGradients(_nScale, 0.);
		std::vector<double> modifiedGradient = fullHessian[p];
		for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) { // First: modify the entries of the gradient corresponding to the actual parameters & collect the gradient for the scale parameters
			int scl = std::get<2>(copy);
			double scaler = std::get<3>(copy);
			if (scl > -1) {
				scaleGradients[scl] += fullHessian[p][std::get<0>(copy)] * fullParams[std::get<1>(copy)]*scaler;
				scaler *= params[nParNonScale + scl];
			}
			modifiedGradient[std::get<1>(copy)] += fullHessian[p][std::get<0>(copy)] * scaler;
		}
		modifiedGradient = getFinalParams(modifiedGradient); // Second: Cut away the fixed parameters from the grandient (this leaves the gradient entries for the scale parameters with odd values)
		for (size_t s = 0; s < _nScale; ++s) {
			modifiedGradient[nParNonScale+s] = scaleGradients[s]; // Third: Fix the grandient entries for the scale parameters 
		}
		modifiedHessian[p] = modifiedGradient;
	}
// The do the same for columns (then turn around this interpretation)
	std::vector<std::vector<double> > retVal(nPar);
	for (size_t p = 0; p < nPar; ++p) {
		std::vector<double> scaleGradients(_nScale, 0.);
		std::vector<double> modifiedGradient(nParTot);
		for (size_t i = 0; i < nParTot; ++i) {
		 	modifiedGradient[i] = modifiedHessian[i][p];
		}
		for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) { // First: modify the entries of the gradient corresponding to the actual parameters & collect the gradient for the scale parameters
			int scl = std::get<2>(copy);
			double scaler = std::get<3>(copy);
			if (scl > -1) {
				scaleGradients[scl] += modifiedHessian[std::get<0>(copy)][p] * fullParams[std::get<1>(copy)]*scaler;
				scaler *= params[nParNonScale + scl];
			}
			modifiedGradient[std::get<1>(copy)] += modifiedHessian[std::get<0>(copy)][p] * scaler;
		}
		modifiedGradient = getFinalParams(modifiedGradient); // Second: Cut away the fixed parameters from the grandient (this leaves the gradient entries for the scale parameters with odd values)
		for (size_t s = 0; s < _nScale; ++s) {
			modifiedGradient[nParNonScale+s] = scaleGradients[s]; // Third: Fix the grandient entries for the scale parameters 
		}
		retVal[p] = modifiedGradient;
	}
	return retVal;
}


size_t logLikelihoodBase::getNparTot() const {
	size_t retVal = 2*_nAmpl;
	return retVal;
}

size_t logLikelihoodBase::getNpar() const {
	size_t retVal = getNparTot() - _fixedParameters.size() + _nScale;
	return retVal;
}

std::vector<double> logLikelihoodBase::getFullParameters(const std::vector<double>& params) const {
	const size_t nPar = getNpar();
	const size_t nParTot = getNparTot();
	if (params.size() != nPar) {
		std::cout << "logLikelihoodBase::getFullParameters(...): ERROR: Parameter size mismatch: " << params.size() << " should be " << nPar << std::endl;
		throw;
	}
	std::vector<double> retVal(nParTot);
	size_t nextFixed = nParTot;
	if (_fixedParameters.size() > 0) {
		nextFixed = _fixedParameters[0].first;
	}
	size_t fixCount  = 0;
	size_t parCount  = 0;
	for (size_t p = 0; p < nParTot; ++p) {
		if (p == nextFixed) {
			retVal[p] = _fixedParameters[fixCount].second;
			++fixCount;
			if (fixCount == _fixedParameters.size()) {
				nextFixed = nParTot;
			} else {
				nextFixed = _fixedParameters[fixCount].first;
			}
		} else {
			retVal[p] = params[parCount];
			++parCount;
		}
	}
	for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) {
		retVal[std::get<0>(copy)] = retVal[std::get<1>(copy)] * std::get<3>(copy);
		if (std::get<2>(copy) > -1) {
			retVal[std::get<0>(copy)] *= params[parCount + std::get<2>(copy)];
		}	
	}
//	parCount += _nScale; // do this, if a third type of parameter should appear
	return retVal;
}

std::vector<double> logLikelihoodBase::getFinalParams(const std::vector<double>& params) const {
	const size_t nPar = getNpar();
	const size_t nParTot = getNparTot();
	if (params.size() != nParTot) {
		std::cout << "logLikelihoodBase::getFinalParams(...): ERROR: Parameter size mismatch: " << params.size() << " should be " << nParTot << std::endl;
		throw;
	}
	std::vector<double> retVal(nPar);
	size_t nextFixed = nParTot;
	if (_fixedParameters.size() > 0) {
		nextFixed = _fixedParameters[0].first;
	}
	size_t fixCount  = 0;
	size_t parCount  = 0;
	for (size_t p = 0; p < nParTot; ++p) {
		if (p == nextFixed) {
			++fixCount;
			if (fixCount == _fixedParameters.size()) {
				nextFixed = nParTot;
			} else {
				nextFixed = _fixedParameters[fixCount].first;
			}
		} else {
			retVal[parCount] = params[p];
			++parCount;
		}
	}
	std::vector<bool> scalesSet(_nScale, false);
	size_t nScalesSet = 0;
	for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) {
		int scl = std::get<2>(copy);
		if (scl > -1) {
			if (!scalesSet[scl]) {
				double from = params[std::get<1>(copy)];
				if (from == 0.) {
					continue;
				}
				double to   = params[std::get<0>(copy)];
				retVal[parCount+scl] = to/from/std::get<3>(copy);
				scalesSet[scl] = true;
				++nScalesSet;
			}
		}
		if (nScalesSet == _nScale) {
			break;
		}
	}
//	parCount += _nScale; // do this, if a third type of parameter should appear
	return retVal;
}

bool logLikelihoodBase::fixParameter(size_t nPar, double value) {
	if (nPar >= getNparTot()) {
		std::cout << "logLikelihoodBase::fixParameter(...): ERROR: Parameter number too large." << std::endl;
		return false;
	}
	std::vector<std::pair<size_t, double> > newFixedPars;
	bool appended = false;
	for (std::pair<size_t, double> fixed : _fixedParameters) {
		if (fixed.first == nPar) {
			std::cout << "logLikelihoodBase::fixParameter(...): ERROR: Parameter " << nPar << " already fixed" << std::endl;
			return false;
		}
		if (fixed.first > nPar and !appended) {
			newFixedPars.push_back(std::pair<size_t, double>(nPar, value));
			appended = true;
		}
		newFixedPars.push_back(fixed);
	}
	if (!appended)  {
		newFixedPars.push_back(std::pair<size_t, double>(nPar, value));
	}
	_fixedParameters = newFixedPars;
	return true;
}

bool logLikelihoodBase::addCopyParameter(size_t nPar, size_t copyFrom, int nScale, double scaler) {
	if (nScale < -1) {
		std::cout << "logLikelihoodBase::addCopyParameter(...): ERROR: Invalid value for nScale" << std::endl;
	}
	if (scaler == 0.) {
		std::cout << "logLikelihoodBase::addCopyParameter(...): ERROR: Encountered a scaler of zero (fix parameter to zero instead)" << std::endl;
		return false;
	}
	const size_t nParTot = getNparTot();
	if (nPar >= nParTot){
		std::cout << "logLikelihoodBase::addCopyParameter(...): ERROR: nPar > nParTot" << std::endl;
		return false;
	}
	if (copyFrom >= nParTot) {
		std::cout << "logLikelihoodBase::addCopyParameter(...): ERROR: nPar > nParTot" << std::endl;
		return false;		
	}
	for (std::tuple<size_t, size_t, int, double> copy : _copyParameters) {
		if (std::get<0>(copy) == copyFrom) {
			std::cout << "logLikelihoodBase::addCopyParameter(...): ERROR: Cannot copy from a copied parameter" << std::endl;
			return false;
		}
	}
	if (!fixParameter(nPar, std::numeric_limits<double>::quiet_NaN())) { // Set to NaN, since it will have to be overwritten
		std::cout << "logLikelihoodBase::addCopyParameter(...): ERROR: Could not fix the parameter to NaN" << std::endl;
		return false;
	}
	if (nScale > -1) {
		if (nScale > (int)_nScale) {
			std::cout << "logLikelihoodBase::addCopyParameter(...): ERROR: nScale skipped (Would introduce an unused scale parameter) (nScale = " << nScale << ", _nScale = "<< _nScale << ")" << std::endl;
			return false;
		}
		if ((size_t)nScale == _nScale) {
			++_nScale; // Add additional scale parameter
		}
		
	}
	_copyParameters.push_back(std::tuple<size_t, size_t, int, double>(nPar, copyFrom, nScale, scaler));
	return true;
}

std::vector<std::complex<double> > logLikelihoodBase::fullParamsToProdAmps(const std::vector<double>& params) const {
	if (params.size() != getNparTot()) {
		std::cout << "logLikelihoodBase::fullParamsToProdAmps(...): ERROR: Parameter size mismatch: " << params.size() << " should be " << getNparTot() << std::endl;
		throw;
	}

	std::vector<std::complex<double> > prodAmps(_nAmpl);
	for(size_t a = 0; a < _nAmpl; ++a) {
		prodAmps[a] = std::complex<double>(params[2*a], params[2*a+1]);
	}
	return prodAmps;
}


std::vector<double> logLikelihoodBase::prodAmpsToFullParams(const std::vector<std::complex<double> >& prodAmps) const {
	if (prodAmps.size() != _nAmpl) {
		std::cout << "logLikelihoodBase::prodAmpsToFullParams(...): ERROR: Parameter size mismatch: " << prodAmps.size() << " should be " << _nAmpl << std::endl;
		throw;
	}
	std::vector<double> params(2*_nAmpl);
	for (size_t a = 0; a < _nAmpl; ++a) {
		params[2*a  ] = prodAmps[a].real();
		params[2*a+1] = prodAmps[a].imag();
	}
	return params;
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
	logLikelihoodBase(amplitudes.size(), integral->kinSignature(), integral), _nSect(1), _maximumIncoherentIntensity(0.), _amplitudes(amplitudes), _points(), _amplitudeCoherenceBorders(0), _contributingWaves(0),
        _nStore(0), _storageLocation(0), _lastPar(), _storeEvals(), _storeSteps() {
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
//	std::ofstream outFile;
//	std::string outFileName = "was_geschicht_"+std::to_string(_nAmpl)+".deleteMe";
//	outFile.open(outFileName.c_str());

	double apa = 0.;
	for (std::complex<double> PA : prodAmps) {
		apa += norm(PA);
	}
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
//double imll = ll;
//std::cout << " i n t e r m e d i a r y " << imll << std::endl;

	if (_extended) {
		ll -= _integral->totalIntensity(prodAmps, true);
	} else {
		ll -= log(_integral->totalIntensity(prodAmps, true))*(double)_nPoints;
	}
//std::cout << " d e l t a l l " << ll - imll << std::endl;
//std::cout << " - - - - - - - - - - - - - -" << std::endl;
	if (_nStore > 0) {
		if (_lastPar.size() > 0) {
			_storeEvals[_storageLocation] = -ll;
			double step = 0.;
			for (size_t a = 0; a < prodAmps.size(); ++a) {
				step += norm(prodAmps[a] - _lastPar[a]);
			}
			_storeSteps[_storageLocation] = pow(step,.5);
			++_storageLocation;
			if (_storageLocation == _nStore) {
				_storageLocation = 0;
			}
		}
		_lastPar = prodAmps;
	}
	if (_maximumIncoherentIntensity >0. ){
		double incoherentIntensity = 0.;
		for (std::complex<double> pa : prodAmps) {
			incoherentIntensity += norm(pa);
		}
		if (incoherentIntensity > _maximumIncoherentIntensity) {
//			std::cout << " es ist soweit " << incoherentIntensity << std::endl;
			ll -= pow(incoherentIntensity - _maximumIncoherentIntensity,2);
		}
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
//	_integral->setNonDiagToZero(); 
	std::vector<double> Dintegral = _integral->DtotalIntensity(prodAmps, true);
	
	double totInt = _integral->totalIntensity(prodAmps, true);
	double delta = 1.e-6;
	for (size_t a = 0; a < _nAmpl; ++a) {
		prodAmps[a] += std::complex<double>(delta,0.);
//		std::cout << Dintegral[2*a] << " | | | " << (_integral->totalIntensity(prodAmps, true)-totInt)/delta  << std::endl;
//		Dintegral[2*a  ] = (_integral->totalIntensity(prodAmps, true)-totInt)/delta;
		prodAmps[a] += std::complex<double>(-delta,delta);
//		std::cout << Dintegral[2*a+1] << " | | | " << (_integral->totalIntensity(prodAmps, true)-totInt)/delta << std::endl;
//		Dintegral[2*a+1] = (_integral->totalIntensity(prodAmps, true)-totInt)/delta;
		prodAmps[a] += std::complex<double>(0.,-delta);
	}
	totInt += totInt;
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
	if (_maximumIncoherentIntensity > 0. ){
		double incoherentIntensity = 0.;
		for (std::complex<double> pa : prodAmps) {
			incoherentIntensity += norm(pa);
		}
		if (incoherentIntensity > _maximumIncoherentIntensity) {
			double fakk = 4*(incoherentIntensity - _maximumIncoherentIntensity);
			for (size_t a = 0; a < _nAmpl; ++a) {
				retVal[2*a  ] += fakk*prodAmps[a].real();
				retVal[2*a+1] += fakk*prodAmps[a].imag();
			}
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
	if (_maximumIncoherentIntensity > 0. ){
		double incoherentIntensity = 0.;
		for (std::complex<double> pa : prodAmps) {
			incoherentIntensity += norm(pa);
		}
		if (incoherentIntensity > _maximumIncoherentIntensity) {
			double fakk = 4*(incoherentIntensity - _maximumIncoherentIntensity);
			for (size_t a = 0; a < _nAmpl; ++a) {
				retVal[2*a  ][2*a  ] += fakk;
				retVal[2*a+1][2*a+1] += fakk;
				for (size_t b = 0; b < _nAmpl; ++b) {
					retVal[2*a  ][2*b  ] += 8*prodAmps[a].real()*prodAmps[b].real();
					retVal[2*a+1][2*b  ] += 8*prodAmps[a].imag()*prodAmps[b].real();
					retVal[2*a  ][2*b+1] += 8*prodAmps[a].real()*prodAmps[b].imag();
					retVal[2*a+1][2*b+1] += 8*prodAmps[a].imag()*prodAmps[b].imag();
				}
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
		std::pair<bool, std::complex<double> > diag = _integral->element(a,a, false);
		if (not diag.first) {
			std::cerr << "logLikelihood::loadDataPoints(...): ERROR: Could not get diagonal integral" << std::endl;
			return false;
		}
		double norm = 0.;
		if (not diag.second.real() == 0.) {
			norm = 1./pow(diag.second.real(), .5);
		}
		for (size_t p = 0; p < _nPoints; ++p) {
			_points[p][a] = _amplitudes[a]->eval(dataPoints[p]);
			if (norm == 0. && _points[p][a] != 0.) {
				std::cout << "logLikelihood::loadDataPoints(...): ERROR: Zero-norm wave has non-vanishing amplitude" << std::endl;
				return false;
			}
			_points[p][a] *= norm;
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

bool logLikelihood::setNstore(size_t nStore) {
	_nStore = nStore;
	_storageLocation = 0;
	_storeEvals = std::vector<double>(_nStore,0.);
	_storeSteps = std::vector<double>(_nStore,0.);
	return true;
}

std::pair<std::vector<double>, std::vector<double> > logLikelihood::getStoredPoints() const {
	std::pair<std::vector<double>, std::vector<double> > retVal(std::vector<double>(_nStore,0.), std::vector<double>(_nStore,0.));
	size_t  location = _storageLocation;
	for (size_t s = 0; s < _nStore; ++s) {
		retVal.first[s]  = _storeEvals[location];
		retVal.second[s] = _storeSteps [location];
		if (location > 0) {
			--location;
		} else {
			location = _nStore - 1;
		}
	}
	return retVal;
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
		_eventsPerBin[bin.second.first][bin.second.second]        += 1;
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
		if (found1 and found2) { // This breaks the total check of the binnig, i.e. if a second bin would be found, but ordering is checked in the constuctor anyway...
			break;
		}
	}
	return std::pair<bool, std::pair<size_t, size_t> >(found1 and found2, std::pair<size_t, size_t>(bin1, bin2));
}
