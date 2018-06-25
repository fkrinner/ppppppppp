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
#include"fabiliLL.h"
#include"fabiliLL_openMP.h"
#include"fabili.h"

std::complex<double> BW(const double& m2, const double&  m0, const double& G0) {
	return std::complex<double>(m0*G0,0.)/std::complex<double>(m0*m0 - m2, -m0*G0);
}

int main(int argc, char *argv[]) {
	(void) argv;

	size_t seed = size_t( time(NULL) );
	std::cout << "Seed: " << seed << std::endl;
	srand(seed);

	const int softpionSign = 0;

	const size_t integralPoints = 300000*10*10; // Ten times data set (ecicciency is 0.1)

	std::cout << integralPoints << " integral points" << std::endl;
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
	size_t nAmpl = amplitudes.size();
	std::shared_ptr<threeParticleMassGenerator> gen = std::make_shared<threeParticleMassGenerator>(mD0, fsMasses, S12->kinSignature());
	std::shared_ptr<efficiencyFunction> efficiency  = std::make_shared<BELLE_DtoKpipi_efficiency>();
	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integralPoints, gen, amplitudes, efficiency);

//	std::cout << "Starting integration" << std::endl;
//	integral->integrate();
//	std::cout << "Finished integration" << std::endl;

	if (!integral->loadIntegrals("ps_integral_corr.dat", "ac_integral_corr.dat")) {
		std::cout << "Could not load integrals" << std::endl;
		return 1;
	}
	std::vector<std::vector<double> > dataPoints = getBELLEevents("./BELLE_data.root", softpionSign);
	size_t nData = pow(dataPoints.size(),.5);
/*	const size_t nData = dataPoints.size();
	const size_t nAmpl = amplitudes.size();
*/	std::vector<std::complex<double> > startValues(nAmpl);
	
	// Only crude estimates for the resonance parameters to have nice start vlaues with continuous phase motion

	const double mK_0 = 0.824;
	const double GK_0 = 0.478;
	const double mK_1 = 0.895;
	const double GK_1 = 0.047;
	const double mK_2 = 1.430;
	const double GK_2 = 0.109;

	const double mf_0 = 0.500;
	const double Gf_0 = 0.500;
	const double mRho = 0.77526;
	const double Grho = 0.1491;
	const double mf_2 = 1.2755;
	const double Gf_2 = 0.1867;

	size_t count = 0;
	double m2;
	std::complex<double> coeff(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpi.size()-1;++b) {
		m2 = (binningKpi[b] + binningKpi[b+1])/2;
		startValues[count] = coeff*BW(m2, mK_0, GK_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpi.size()-1;++b) {
		m2 = (binningKpi[b] + binningKpi[b+1])/2;
		startValues[count] = coeff*BW(m2, mK_1, GK_1);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpi.size()-1;++b) {
		m2 = (binningKpi[b] + binningKpi[b+1])/2;
		startValues[count] = coeff*BW(m2, mK_2, GK_2);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpi.size()-1;++b) {
		m2 = (binningKpi[b] + binningKpi[b+1])/2;
		startValues[count] = coeff*BW(m2, mK_0, GK_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpi.size()-1;++b) {
		m2 = (binningKpi[b] + binningKpi[b+1])/2;
		startValues[count] = coeff*BW(m2, mK_1, GK_1);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpi.size()-1;++b) {
		m2 = (binningKpi[b] + binningKpi[b+1])/2;
		startValues[count] = coeff*BW(m2, mK_2, GK_2);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningPiPi.size()-1;++b) {
		m2 = (binningPiPi[b] + binningPiPi[b+1])/2;
		startValues[count] = coeff*BW(m2, mf_0, Gf_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningPiPi.size()-1;++b) {
		m2 = (binningPiPi[b] + binningPiPi[b+1])/2;
		startValues[count] = coeff*BW(m2, mRho, Grho);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningPiPi.size()-1;++b) {
		m2 = (binningPiPi[b] + binningPiPi[b+1])/2;
		startValues[count] = coeff*BW(m2, mf_2, Gf_2);
		++count;
	}
	if (count != nAmpl) {
		std::cout << "Number of initialized start values does not match " << count << " != " << nAmpl;
		return 1;
	}
	std::cout << "Start fitting " << argc << std::endl;

	std::ofstream outFile;

//	outFile.open("startValueTest.deleteMe");
//	for (size_t a = 0; a < startValues.size(); ++a) {
//		outFile << a << " " << startValues[a].real() << " " << startValues[a].imag() << std::endl;
//	}
//	outFile.close();

	bool parallel = false;
	if (argc > 1) {
		parallel = true;
	}

	if (parallel) {
		std::cout << "Running in parallel mode (results not written so far)" << std::endl;
//		std::shared_ptr<fabiliLL> ll = std::make_shared<fabiliLL>(amplitudes, integral);
		std::shared_ptr<fabiliLL_openMP> ll = std::make_shared<fabiliLL_openMP>(amplitudes, integral);
		ll->loadDataPoints(dataPoints);
		ll->setExtended(true);

		fabili ff(ll);
		ff.setStepSize(.1);
		std::vector<double> ppp;
		for (std::complex<double>& a : startValues) {
			ppp.push_back(a.real());
			ppp.push_back(a.imag());
		}

		ff.minimize(ppp);

		std::cout << "Finished" << std::endl;
	} else {
		std::cout << "Running in normal mode" << std::endl;

		logLikelihood ll(amplitudes, integral);
		ll.loadDataPoints(dataPoints);
		ll.setExtended(true);

		std::pair<double, std::vector<std::complex<double> > > retVal = ll.fitNlopt(startValues);


		std::string outFileName = std::string("./BELLEfit_pi") + std::to_string(softpionSign) + std::string("_nsv_seed")+ std::to_string(seed) + std::string("_ll") + std::to_string(retVal.first) + std::string(".dat");
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
	}
	return 0;
}
