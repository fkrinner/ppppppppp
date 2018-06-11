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
#include"efficiencyFunction.h"
#include"modelAmplitude.h"
#include"modelGenerator.h"
#include"utils.h"
#include"logLikelihood.h"
#include"constants.h"
#include"getBELLEdata.h"

int main() {
	size_t seed = size_t( time(NULL) );
	std::cout << "Seed: " << seed << std::endl;
	srand(seed);

	const int softpionSign = 0;

	const size_t nTries         = 10;
	const size_t integralPoints = 300000*10*10; // Ten times data set (ecicciency is 0.1)

	std::cout << nTries << " tries; " << integralPoints << " integral points" << std::endl;
	std::vector<double> fsMasses = {mPi, mKs, mPi};

	std::shared_ptr<BELLE_S> S12 = std::make_shared<BELLE_S>(12);	
	std::shared_ptr<BELLE_S> S13 = std::make_shared<BELLE_S>(13);	
	std::shared_ptr<BELLE_S> S23 = std::make_shared<BELLE_S>(23);	

	std::shared_ptr<BELLE_P> P12 = std::make_shared<BELLE_P>(12);	
	std::shared_ptr<BELLE_P> P13 = std::make_shared<BELLE_P>(13);	
	std::shared_ptr<BELLE_P> P23 = std::make_shared<BELLE_P>(23);	

	std::shared_ptr<BELLE_D> D12 = std::make_shared<BELLE_D>(12);	
	std::shared_ptr<BELLE_D> D13 = std::make_shared<BELLE_D>(13);	
	std::shared_ptr<BELLE_D> D23 = std::make_shared<BELLE_D>(23);	
	
	const double binWidth = 0.04;

	double m = mPi+ mKs;	
	std::vector<double> binningKpi = {m*m};
	while (m < mD0 - mPi) {
		m += binWidth;
		binningKpi.push_back(m*m);
	};

	m = 2*mPi;	
	std::vector<double> binningPiPi = {m*m};
	while (m < mD0 - mKs) {
		m += binWidth;
		binningPiPi.push_back(m*m);
	};

	std::vector<std::shared_ptr<amplitude> > amplitudes = {};


// - - - - - - - K piRight S-wave
	for (size_t b = 0; b < binningKpi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpi[b], binningKpi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(12, std::string("KpiRight[") + std::to_string(b) + std::string("]PiS"), step, S12, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - K piRight P-wave
	for (size_t b = 0; b < binningKpi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpi[b], binningKpi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(12, std::string("KpiRight[") + std::to_string(b) + std::string("]PiP"), step, P12, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - K piRight D-wave
	for (size_t b = 0; b < binningKpi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpi[b], binningKpi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(12, std::string("KpiRight[") + std::to_string(b) + std::string("]PiD"), step, D12, fsMasses);
		amplitudes.push_back(stepWave);
	}
// ===========================================
// - - - - - - - K piWrong S-wave
	for (size_t b = 0; b < binningKpi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpi[b], binningKpi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(23, std::string("KpiWrong[") + std::to_string(b) + std::string("]PiS"), step, S23, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - K piWrong P-wave
	for (size_t b = 0; b < binningKpi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpi[b], binningKpi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(23, std::string("KpiWrong[") + std::to_string(b) + std::string("]PiP"), step, P23, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - K piRight D-wave
	for (size_t b = 0; b < binningKpi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpi[b], binningKpi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(23, std::string("KpiWrong[") + std::to_string(b) + std::string("]PiD"), step, D23, fsMasses);
		amplitudes.push_back(stepWave);
	}
// ============================================
// - - - - - - - pi pi S-wave
	for (size_t b = 0; b < binningPiPi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningPiPi[b], binningPiPi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(13, std::string("piPi[") + std::to_string(b) + std::string("]KS"), step, S13, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - pi pi P-wave
	for (size_t b = 0; b < binningPiPi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningPiPi[b], binningPiPi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(13, std::string("piPi[") + std::to_string(b) + std::string("]KP"), step, P13, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - pi pi D-wave
	for (size_t b = 0; b < binningPiPi.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningPiPi[b], binningPiPi[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(13, std::string("piPi[") + std::to_string(b) + std::string("]KD"), step, D13, fsMasses);
		amplitudes.push_back(stepWave);
	}

	std::cout << amplitudes.size() << " waves in the model" << std::endl;

	std::shared_ptr<threeParticleMassGenerator> gen = std::make_shared<threeParticleMassGenerator>(mD0, fsMasses, S12->kinSignature());
	std::shared_ptr<efficiencyFunction> efficiency  = std::make_shared<BELLE_DtoKpipi_efficiency>();
	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integralPoints, gen, amplitudes, efficiency);

//	std::cout << "Starting integration" << std::endl;
//	integral->integrate();
//	std::cout << "Finished integration" << std::endl;

	if (!integral->loadIntegrals("ps_integral.dat", "ac_integral.dat")) {
		std::cout << "Could not load integrals" << std::endl;
		return 1;
	}
	std::vector<std::vector<double> > dataPoints = getBELLEevents("./BELLE_data.root", softpionSign);
	logLikelihood ll(amplitudes, integral);
	ll.loadDataPoints(dataPoints);
	ll.setExtended(true);
	const size_t nData = dataPoints.size();
	const size_t nAmpl = amplitudes.size();
	std::vector<std::complex<double> > startValues(nAmpl);
	for (size_t a = 0; a < nAmpl; ++a) {
		startValues[a] = std::complex<double>(pow(random()*nData,.5), pow(random()*nData,.5));
	}
	std::cout << "Start fitting" << std::endl;
	std::pair<double, std::vector<std::complex<double> > > retVal = ll.fitNlopt(startValues);
	std::cout << "Finished" << std::endl;
	std::ofstream outFile;
	std::string outFileName = std::string("./BELLEfit_pi") + std::to_string(softpionSign) + std::string("_seed")+ std::to_string(seed) + std::string("_ll") + std::to_string(retVal.first) + std::string(".dat");
	outFile.open(outFileName.c_str());
	for (std::complex<double>& amp : retVal.second) {
		outFile << amp << std::endl;
	}
	outFile.close();

	std::vector<std::vector<double> > hessian = ll.DDeval(retVal.second);
	outFileName = std::string("./BELLEfit_hessian_") + std::to_string(seed) + std::string(".dat");
	outFile.open(outFileName.c_str());
	for (size_t i = 0; i < 2*ll.nAmpl(); ++i) {
	for (size_t j = 0; j < 2*ll.nAmpl(); ++j) {
			outFile << hessian[i][j] << " ";
		}
		outFile << std::endl;
	}
	outFile.close();

	return 0;
}
