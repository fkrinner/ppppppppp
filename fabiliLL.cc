#include"fabiliLL.h"

fabiliLL::fabiliLL(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral) :
	logLikelihood(amplitudes,integral), 
	fabiliFunction() {}

double fabiliLL::scalarEval(const std::vector<double>& parameters) const {
        std::vector<std::complex<double> > prodAmps(_nAmpl);
        size_t skip       = 1;
        if (_fixFirstPhase) {
                prodAmps[0] = std::complex<double>(parameters[0], 0.);
        } else {
                ++skip;
                prodAmps[0] = std::complex<double>(parameters[0], parameters[1]);
        }
        for (size_t a = 1; a < _nAmpl; ++a) {
                prodAmps[a] = std::complex<double>(parameters[2*a-2+skip], parameters[2*a-1+skip]);
        }

	return logLikelihood::eval(prodAmps);
}

evalType fabiliLL::eval(const std::vector<double>& parameters) const {
	const size_t dim = fabiliLL::dim();
	if (parameters.size() != dim) {
		std::cout << "Wrong size of parameters" << std::endl;
		throw;
	}
	std::vector<std::complex<double> > prodAmps(_nAmpl);
	size_t skip       = 1;
	if (_fixFirstPhase) {
		prodAmps[0] = std::complex<double>(parameters[0], 0.);
	} else {
		++skip;
		prodAmps[0] = std::complex<double>(parameters[0], parameters[1]);
	}
	for (size_t a = 1; a < _nAmpl; ++a) {
		prodAmps[a] = std::complex<double>(parameters[2*a-2+skip], parameters[2*a-1+skip]);
	}
	double ll = 0.;
	std::vector<double> Dll(2*_nAmpl, 0.);
	std::vector<std::vector<double> > DDll(2*_nAmpl, std::vector<double>(2*_nAmpl, 0.));

	for (size_t p = 0; p < _nPoints; ++p) {
		std::complex<double> ampl (0.,0.);
		for (size_t a = 0; a < _nAmpl; ++a) {
			ampl += _points[p][a] * prodAmps[a];
		}
		double intens = norm(ampl);
		ll -= log(intens);
		ampl = std::conj(ampl); // Will only be used conjugated
		for (size_t ai = 0; ai < _nAmpl; ++ai) {
			std::complex<double> factori = _points[p][ai] * ampl/intens * 2.;
			Dll[2*ai  ] -= factori.real(); // One factor of -1 since, the NEGATIVE likelihood is used 
			Dll[2*ai+1] += factori.imag(); // Two factors of -1: Complex i*i = -1 and since the NEGATIVE likelihood is used
			for (size_t aj = 0; aj < _nAmpl; ++aj) {
				std::complex<double> factor = _points[p][ai] * std::conj(_points[p][aj])/intens*2.;
				DDll[2*ai  ][2*aj  ] -= factor.real();
				DDll[2*ai  ][2*aj+1] -= factor.imag();
				DDll[2*ai+1][2*aj  ] += factor.imag();
				DDll[2*ai+1][2*aj+1] -= factor.real();

				std::complex<double> factorj = -_points[p][aj] * ampl/intens *2.; // -1/intens and then the same factor;

				DDll[2*ai  ][2*aj  ] -= factori.real() * factorj.real();
				DDll[2*ai  ][2*aj+1] += factori.real() * factorj.imag();
				DDll[2*ai+1][2*aj  ] += factori.imag() * factorj.real();
				DDll[2*ai+1][2*aj+1] -= factori.imag() * factorj.imag();
			}
		}
	}
	double integralIntensity      = _integral->totalIntensity(prodAmps, true);
	std::vector<double> Dintegral = _integral->DtotalIntensity(prodAmps, true);
	std::vector<std::vector<double> > DDintegral = _integral->DDtotalIntensity(prodAmps, true);
	if (_extended) {
		ll += integralIntensity;
		for (size_t i = 0; i < 2*_nAmpl; ++i) {
			Dll[i] += Dintegral[i];
			for (size_t j = 0; j < 2*_nAmpl; ++j) {
				DDll[i][j] += DDintegral[i][j];
			}
		}
	} else {
		ll += log(integralIntensity)*(double)_nPoints;
		for (size_t i = 0; i< 2*_nAmpl; ++i) {
			Dll[i] += Dintegral[i]/integralIntensity*(double)_nPoints;
			for (size_t j = 0; j < 2*_nAmpl; ++j) {
				DDll[i][j] += (DDintegral[i][j]/integralIntensity - Dintegral[i]*Dintegral[j]/integralIntensity/integralIntensity)*(double)_nPoints;
			}
		}
	}
	++_nCalls;
	evalType retVal;
	retVal.value        = ll;
	retVal.gradient     = Eigen::VectorXd::Zero(dim);
	retVal.hessian      = Eigen::MatrixXd::Zero(dim, dim);
	for (size_t ai = 0; ai < _nAmpl; ++ai) {
		size_t iR = 2*ai;
		size_t iI = 2*ai+1;
		if (ai > 0. && _fixFirstPhase) {
			--iR;
			--iI;
		}
		retVal.gradient(iR) = Dll[2*ai  ];
		retVal.gradient(iI) = Dll[2*ai+1];
		for (size_t aj = 0; aj < _nAmpl; ++aj) {
			size_t jR = 2*aj;
			size_t jI = 2*aj+1;
			if (aj > 0 && _fixFirstPhase) {
				--jR;
				--jI;
			}
			retVal.hessian(iR,jR) = DDll[2*ai  ][2*aj  ];
			retVal.hessian(iR,jI) = DDll[2*ai  ][2*aj+1];
			retVal.hessian(iI,jR) = DDll[2*ai+1][2*aj  ];
			retVal.hessian(iI,jI) = DDll[2*ai+1][2*aj+1];
		}
	}
	return retVal;
}
