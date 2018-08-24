#include"fabiliLL_openMP.h"
#include"omp.h"

fabiliLL_openMP::fabiliLL_openMP(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral) :
	logLikelihood(amplitudes,integral), 
	fabiliFunction() {}

double fabiliLL_openMP::scalarEval(const std::vector<double>& parameters) const {
        std::vector<std::complex<double> > prodAmps = fullParamsToProdAmps(getFullParameters(parameters));
	return logLikelihood::eval(prodAmps);
}

evalType fabiliLL_openMP::eval(const std::vector<double>& parameters) const {
	const size_t nPar = getNpar();
	if (parameters.size() != nPar) {
		std::cout << "Wrong size of parameters" << std::endl;
		throw;
	}
        std::vector<std::complex<double> > prodAmps = fullParamsToProdAmps(getFullParameters(parameters));
	const size_t maxNthread = omp_get_max_threads();
//	std::cout << "maximum number of threads is " << maxNthread << std::endl;

	std::vector<double> ll_thread = std::vector<double>(maxNthread, 0.);
	std::vector<std::vector<double> > Dll_thread(maxNthread, std::vector<double>(2*_nAmpl, 0.));
	std::vector<std::vector<std::vector<double> > > DDll_thread(maxNthread, std::vector<std::vector<double> >(2*_nAmpl, std::vector<double>(2*_nAmpl, 0.)));
#pragma omp parallel for
	for (size_t p = 0; p < _nPoints; ++p) {
		const size_t ID = omp_get_thread_num();
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
		ll_thread[ID] -= log(intens);
//		#pragma omp parallel for

		size_t sector_i = 0;
		size_t upperSectorBorder_i = _amplitudeCoherenceBorders[1];		
		for (size_t ai : _contributingWaves[p]) {
			if (ai >= upperSectorBorder_i) {
				sector_i = getSector(ai);
				upperSectorBorder_i = _amplitudeCoherenceBorders[sector_i+1];
			}
			std::complex<double> factori = _points[p][ai] * sectorAmpls[sector_i]/intens * 2.;
			Dll_thread[ID][2*ai  ] -= factori.real(); // One factor of -1 since, the NEGATIVE likelihood is used 
			Dll_thread[ID][2*ai+1] += factori.imag(); // Two factors of -1: Complex i*i = -1 and since the NEGATIVE likelihood is used

			size_t sector_j = 0;
			size_t upperSectorBorder_j = _amplitudeCoherenceBorders[1];
			for (size_t aj : _contributingWaves[p]) {
				if (aj >= upperSectorBorder_j) {
					sector_j = getSector(aj);
					upperSectorBorder_j = _amplitudeCoherenceBorders[sector_j+1];
				}
				if (sector_i == sector_j) {
					std::complex<double> factor = _points[p][ai] * std::conj(_points[p][aj])/intens*2.;
					DDll_thread[ID][2*ai  ][2*aj  ] -= factor.real();
					DDll_thread[ID][2*ai  ][2*aj+1] -= factor.imag();
					DDll_thread[ID][2*ai+1][2*aj  ] += factor.imag();
					DDll_thread[ID][2*ai+1][2*aj+1] -= factor.real();
				}
				std::complex<double> factorj = -_points[p][aj] * sectorAmpls[sector_j]/intens *2.; // -1/intens and then the same factor;

				DDll_thread[ID][2*ai  ][2*aj  ] -= factori.real() * factorj.real();
				DDll_thread[ID][2*ai  ][2*aj+1] += factori.real() * factorj.imag();
				DDll_thread[ID][2*ai+1][2*aj  ] += factori.imag() * factorj.real();
				DDll_thread[ID][2*ai+1][2*aj+1] -= factori.imag() * factorj.imag();
			}
		}
	}
	double ll = 0.;
	std::vector<double> Dll(2*_nAmpl, 0.);
	std::vector<std::vector<double> >DDll(2*_nAmpl, std::vector<double>(2*_nAmpl, 0.));
	for (size_t ID = 0; ID < maxNthread; ++ID) {
		ll += ll_thread[ID];
		for (size_t a1 = 0; a1 < 2*_nAmpl; ++a1) {
			Dll[a1] += Dll_thread[ID][a1];
			for (size_t a2 = 0; a2 < 2*_nAmpl; ++a2) {
				DDll[a1][a2] += DDll_thread[ID][a1][a2];
			}
		}
	}


	double integralIntensity      = _integral->totalIntensity(prodAmps, true);
//	std::cout << "integral inensity " << integralIntensity << std::endl;
	std::vector<double> Dintegral = _integral->DtotalIntensity(prodAmps, true);
	std::vector<std::vector<double> > DDintegral = _integral->DDtotalIntensity(prodAmps,true);
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

	DDll                = makeFinalHessian(parameters, DDll); // To this first, since the full  gradient is needed
	Dll                 = makeFinalGradient(parameters, Dll);

	retVal.gradient     = Eigen::VectorXd::Zero(nPar);
	retVal.hessian      = Eigen::MatrixXd::Zero(nPar, nPar);
	for (size_t i = 0; i < nPar; ++i) {
		retVal.gradient(i) = Dll[i];
		retVal.gradient(i) = Dll[i];
		for (size_t j = 0; j < nPar; ++j) {
			retVal.hessian(i,j) = DDll[i][j];
		}
	}
	return retVal;
}
