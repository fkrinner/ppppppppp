#include<iostream>
#include<memory>
#include<vector>
#include<chrono>
#include<ctime>
#include<fstream>
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
	const double mPi =  .13957018;
	const double mDp = 1.8696200;

	const double width = .04;

	const size_t integralPoints = 1000000*(0.04/width)*(0.04/width);

	double m = 2*mPi;

	std::vector<double> binning = {m*m};
	while (m < mDp - mPi) {
		m += width;
		binning.push_back(m*m);
	}
	const size_t nBins = binning.size() - 1;

	std::vector<double> fsMasses = {mPi, mPi, mPi};

	std::shared_ptr<simpleBW> BW  = std::make_shared<simpleBW>( 1., .2 );

	std::shared_ptr<threeParticleMassGenerator> generator = std::make_shared<threeParticleMassGenerator>(mDp, fsMasses);
	std::shared_ptr<sameMassZeroS> Sangle = std::make_shared<sameMassZeroS>(mPi);
	std::shared_ptr<sameMassOneP>  Pangle = std::make_shared<sameMassOneP>(mPi);


	std::shared_ptr<threeParticleIsobaricAmplitude> Swave = std::make_shared<threeParticleIsobaricAmplitude>(true, std::string("SBW"), BW, Sangle);
	std::shared_ptr<threeParticleIsobaricAmplitude> Pwave = std::make_shared<threeParticleIsobaricAmplitude>(true, std::string("PBW"), BW, Pangle);
	std::vector<std::shared_ptr<amplitude> > amplitudes;

	std::shared_ptr<efficiencyFunction> efficiency = std::make_shared<threeParticlPerfectEfficiency>(Swave->kinSignature());

	for (size_t b = 0; b < nBins; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binning[b], binning[b+1]);
		std::shared_ptr<threeParticleIsobaricAmplitude> stepWave = std::make_shared<threeParticleIsobaricAmplitude>(true, std::string("Swave[") + std::to_string(b) + std::string("]"), step, Sangle);
		amplitudes.push_back(stepWave);
	}
	for (size_t b = 0; b < nBins; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binning[b], binning[b+1]);
		std::shared_ptr<threeParticleIsobaricAmplitude> stepWave = std::make_shared<threeParticleIsobaricAmplitude>(true, std::string("Pwave[") + std::to_string(b) + std::string("]"), step, Pangle);
		amplitudes.push_back(stepWave);
	}
//	amplitudes.push_back(Swave);
//	amplitudes.push_back(Pwave);
	
	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integralPoints, generator, amplitudes, efficiency);
	integral->integrate();
	integral->writeToFile("./integrals0mp.dat");
	return 0;
}
