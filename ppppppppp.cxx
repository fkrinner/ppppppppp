#include<iostream>
#include<memory>
#include<vector>
#include <chrono>
#include <ctime>
#include"massShape.h" 
#include"angularDependence.h"
#include"amplitude.h"
#include"integrator.h"
#include"generator.h"
#include"modelAmplitude.h"
#include"modelGenerator.h"
#include"utils.h"
#include"logLikelihood.h"

int main() {
	size_t seed = size_t( time(NULL) );
	srand( seed );
//	utils::opening();

	const double mPi =  .13957018;
	const double mD0 = 1.86484;

	size_t nPoints = 10000;
	size_t integralPoints = 1000000;

	bool freedIsobar = true;

	std::vector<double> fsMasses = {mPi, mPi, mPi};
	std::shared_ptr<threeParticleMassGenerator> gen = std::make_shared<threeParticleMassGenerator>(mD0, fsMasses);

	std::shared_ptr<simpleBW> f0  = std::make_shared<simpleBW>(.98, .1);
	std::shared_ptr<simpleBW> rho = std::make_shared<simpleBW>(.77, .16);
	std::shared_ptr<simpleBW> f2  = std::make_shared<simpleBW>(1.27, .2);

	std::shared_ptr<sameMassZeroS> Sangle = std::make_shared<sameMassZeroS>(mPi);
	std::shared_ptr<sameMassOneP>  Pangle = std::make_shared<sameMassOneP>(mPi);
	std::shared_ptr<sameMassTwoD>  Dangle = std::make_shared<sameMassTwoD>(mPi);

	std::shared_ptr<threeParticleIsobaricAmplitude> Swave = std::make_shared<threeParticleIsobaricAmplitude>(true, "0mp0ppPiS", f0,  Sangle);
	std::shared_ptr<threeParticleIsobaricAmplitude> Pwave = std::make_shared<threeParticleIsobaricAmplitude>(true, "0mp1mmPiP", rho, Pangle);
	std::shared_ptr<threeParticleIsobaricAmplitude> Dwave = std::make_shared<threeParticleIsobaricAmplitude>(true, "0mp2ppPiP", f2,  Dangle);

	std::vector<std::shared_ptr<amplitude> > amplitudes = {Dwave, Pwave, Swave};
	size_t nAmpl = amplitudes.size();
	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integralPoints, gen, amplitudes);
	std::cout << "Starting integration" << std::endl;
	integral->integrate();
	std::cout << "Finished integration" << std::endl;

	std::vector<double> normalizations;
	std::vector<std::complex<double> > transitionAmplitudes;
	for (size_t a = 0; a < nAmpl; ++a) {
		std::pair<bool, std::complex<double> > diag = integral->element(a,a);
		if (not diag.first) {
			std::cerr << "Could not get diagonal element" << std::endl;
			return 1;
		}
		normalizations.push_back(1./pow(diag.second.real(), .5));
		transitionAmplitudes.push_back(std::complex<double>(utils::random2(), utils::random2()));
	}


	std::shared_ptr<modelAmplitude> model = std::make_shared<modelAmplitude>(transitionAmplitudes, amplitudes, normalizations);
	modelGenerator generator(model, gen);
	std::cout << "Starting generation" << std::endl;
	std::vector<std::vector<double> > generatedPoints = generator.generateDataPoints(nPoints, nPoints); // Won't loose anything, if the burn-in is as long as the sample
	std::cout << "Finished generation" << std::endl;

/*	double s23base = mD0*mD0 + 3*mPi*mPi;
	for (std::vector<double>& pp : generatedPoints) {
		double s12 = pp[1];
		double s13 = pp[2];
		double s23 = s23base - s12 - s13;
		std::cout << s12 << " " << s13 << " " << s23 << std::endl;
	}
	return 0;*/

	double mass = 2*mPi;
	std::vector<double> binning(1,mass*mass);
	double binWidth = .04;
	while (mass < mD0 - mPi) {
		mass += binWidth;
		binning.push_back(mass*mass);
	}

	std::vector<std::shared_ptr<amplitude> > amplitudesFit = amplitudes;
	std::shared_ptr<integrator> integralFit = integral;

	if (freedIsobar) {
		amplitudesFit = {Dwave};
		for (size_t b = 0; b < binning.size() - 1; ++b) {
			std::shared_ptr<stepLike> step  = std::make_shared<stepLike>(binning[b],binning[b+1]);
			std::shared_ptr<threeParticleIsobaricAmplitude> stepWave = std::make_shared<threeParticleIsobaricAmplitude>(true, std::string("0mp0pp[") + std::to_string(b) + std::string("]PiS"), step, Sangle);
			amplitudesFit.push_back(stepWave);
		}
		for (size_t b = 0; b < binning.size() - 1; ++b) {
			std::shared_ptr<stepLike> step  = std::make_shared<stepLike>(binning[b],binning[b+1]);
			std::shared_ptr<threeParticleIsobaricAmplitude> stepWave = std::make_shared<threeParticleIsobaricAmplitude>(true, std::string("0mp1mm[") + std::to_string(b) + std::string("]PiP"), step, Pangle);
			amplitudesFit.push_back(stepWave);
		}
		integralFit = std::make_shared<integrator>(integralPoints, gen, amplitudesFit);
		std::cout << "Starting integration 2" << std::endl;
		integralFit->integrate();
		std::cout << "Finished integration 2" << std::endl;
		logLikelihood ll(amplitudesFit, integralFit);
	}

	logLikelihood ll(amplitudesFit, integralFit);

	integralFit->writeToFile("./integral.dat");

	ll.setFixFirstPhase(true);
	std::cout << "Start loading data points" << std::endl;
	ll.loadDataPoints(generatedPoints);
	std::cout << "Finished loading data points" << std::endl;

	double avgAmpl = pow(nPoints, .5);
	size_t nTries = 1;
	double bestLike = 1./0.;
	std::vector<std::complex<double> > bestVals;
	std::cout << "Start fitting" << std::endl;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	std::vector<std::complex<double> >pa(ll.nAmpl());
	for (size_t _ = 0; _ < nTries; ++_) {
		for (size_t a = 0; a < ll.nAmpl(); ++a) {
			pa[a] = avgAmpl * std::complex<double>(utils::random2(), utils::random2());
		}	
		std::pair<double, std::vector<std::complex<double> > > retVal = ll.fit(pa);
		if (retVal.first < bestLike) {
			std::cout << "try " << _ << " bestLike improved from " << bestLike << " to " << retVal.first << std::endl;
			bestLike = retVal.first;
			bestVals = retVal.second;
		}
	}
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout << "Finished fitting " << std::ctime(&end_time) << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
	std::vector<std::complex<double> > doppel(1,transitionAmplitudes[0] * bestVals[1] / transitionAmplitudes[1] / bestVals[0]);
	if (not freedIsobar) {
		doppel.push_back(transitionAmplitudes[1] * bestVals[2] / transitionAmplitudes[2] / bestVals[1]);
	}
	std::cout << "Doppel:";
	for (std::complex<double> d : doppel){
		std::cout << " " << d;
	}
	std::cout << std::endl;

	double intens = integralFit->totalIntensity(bestVals);
	std::cout << intens << " events described, " << nPoints << " points generated" << std::endl;
	std::cout << "The best likelihood is " << bestLike << std::endl;
	if (freedIsobar) {
		for (size_t a = 0; a < ll.nAmpl(); ++a) {
			std::complex<double> amp = bestVals[a];
			std::cout << a << " " << norm(amp) << " " << amp.real() << " " << amp.imag() << std::endl;
		}
	}
	return 0;
}
/* Testcode for hessians...
	double delta = 0.00001;
	std::complex<double> deltaR(delta,0.);
	std::complex<double> deltaI(0.,delta);
	std::vector<double> val = ll.Deval(transitionAmplitudes);
	std::vector<std::vector<double> > hess = ll.DDeval(transitionAmplitudes);
	for(size_t i = 0; i < nAmpl; ++i) {
		transitionAmplitudes[i] += deltaR;
		std::vector<double> DR = ll.Deval(transitionAmplitudes);
		transitionAmplitudes[i] += deltaI - deltaR;
		std::vector<double> DI = ll.Deval(transitionAmplitudes);
		transitionAmplitudes[i] -= deltaI;
		for (size_t j = 0; j < nAmpl; ++j) {
			double D = (DR[2*j] - val[2*j])/delta;
			std::cout << "DRDR " << D << " " << hess[2*i][2*j] << std::endl;
			D = (DR[2*j+1] - val[2*j+1])/delta;
			std::cout << "DIDR " << D << " " << hess[2*i][2*j+1] << std::endl;
			D = (DI[2*j] - val[2*j])/delta;
			std::cout << "DRDI " << D << " " << hess[2*i+1][2*j] << std::endl;
			D = (DI[2*j+1] - val[2*j+1])/delta;
			std::cout << "DIDI " << D << " " << hess[2*i+1][2*j+1] << std::endl;
			std::cout << "-------------------------------------" << std::endl;
		}
		std::cout << "====================================" << std::endl;
	}
	return 0;
*/

