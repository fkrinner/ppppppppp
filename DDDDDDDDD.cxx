#include<iostream>
#include<memory>
#include<vector>
#include<chrono>
#include<ctime>
#include<fstream>
#include<iomanip>
#include<limits>
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
	std::cout.precision(17);
	const size_t integralPoints = 300000*10*10; // Ten times data set (efficiency is 0.1)

	std::ofstream outFile;
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
	std::string integralFileName = "integral_finer_binning";

	double m = mPi+ mKs;
	std::vector<double> binningKpiS = {m*m};
	while (m < mD0 - mPi) {
		if (false) {
			m += binWidth/2;
		} else {
			m += binWidth;
		}
		binningKpiS.push_back(m*m);
	}

	m = mPi+ mKs;
	std::vector<double> binningKpiP = {m*m};
	while (m < mD0 - mPi) {
		if (m > .789 and m < .999) {
			m += binWidth/4;
		} else {
			m += binWidth;
		}
		binningKpiP.push_back(m*m);
	}

	m = mPi+ mKs;
	std::vector<double> binningKpiD = {m*m};
	while (m < mD0 - mPi) {
		if (m > 1.329 and m < 1.529) {
			m += binWidth/2;
		} else {
			m += binWidth;
		}
		binningKpiD.push_back(m*m);
	}

	m = 2*mPi;
	std::vector<double> binningPiPiS = {m*m};
	while (m < mD0 - mKs) {
		if (m > .919 && m < 1.079) {
			m += binWidth/4;
		} else {
			m += binWidth;
		}
		binningPiPiS.push_back(m*m);
	}

	m = 2*mPi;
	std::vector<double> binningPiPiP = {m*m};
	while (m < mD0 - mKs) {
		if (m > .679 && m < .919) {
			m += binWidth/2;
		} else {
			m += binWidth;
		}
		binningPiPiP.push_back(m*m);
	}

	m = 2*mPi;
	std::vector<double> binningPiPiD = {m*m};
	while (m < mD0 - mKs) {
		if (m > 1.159 && m < 1.39) {
			m += binWidth/2;
		} else {
			m += binWidth;
		}
		binningPiPiD.push_back(m*m);
	}

	std::vector<std::shared_ptr<amplitude> > amplitudes = {};
// - - - - - - - K piRight S-wave
	for (size_t b = 0; b < binningKpiS.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpiS[b], binningKpiS[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(12, std::string("KpiRight[") + std::to_string(b) + std::string("]PiS"), step, S12, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - K piRight P-wave
	for (size_t b = 0; b < binningKpiP.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpiP[b], binningKpiP[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(12, std::string("KpiRight[") + std::to_string(b) + std::string("]PiP"), step, P12, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - K piRight D-wave
	for (size_t b = 0; b < binningKpiD.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpiD[b], binningKpiD[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(12, std::string("KpiRight[") + std::to_string(b) + std::string("]PiD"), step, D12, fsMasses);
		amplitudes.push_back(stepWave);
	}
// ===========================================
// - - - - - - - K piWrong S-wave
	for (size_t b = 0; b < binningKpiS.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpiS[b], binningKpiS[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(23, std::string("KpiWrong[") + std::to_string(b) + std::string("]PiS"), step, S23, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - K piWrong P-wave
	for (size_t b = 0; b < binningKpiP.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpiP[b], binningKpiP[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(23, std::string("KpiWrong[") + std::to_string(b) + std::string("]PiP"), step, P23, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - K piRight D-wave
	for (size_t b = 0; b < binningKpiD.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpiD[b], binningKpiD[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(23, std::string("KpiWrong[") + std::to_string(b) + std::string("]PiD"), step, D23, fsMasses);
		amplitudes.push_back(stepWave);
	}
// ============================================
// - - - - - - - pi pi S-wave
	for (size_t b = 0; b < binningPiPiS.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningPiPiS[b], binningPiPiS[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(13, std::string("piPi[") + std::to_string(b) + std::string("]KS"), step, S13, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - pi pi P-wave
	for (size_t b = 0; b < binningPiPiP.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningPiPiP[b], binningPiPiP[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(13, std::string("piPi[") + std::to_string(b) + std::string("]KP"), step, P13, fsMasses);
		amplitudes.push_back(stepWave);
	}
// - - - - - - - pi pi D-wave
	for (size_t b = 0; b < binningPiPiD.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningPiPiD[b], binningPiPiD[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(13, std::string("piPi[") + std::to_string(b) + std::string("]KD"), step, D13, fsMasses);
		amplitudes.push_back(stepWave);
	}

	std::cout << amplitudes.size() << " waves in the model" << std::endl;
	size_t nAmpl                                    = amplitudes.size();
	std::shared_ptr<threeParticleMassGenerator> gen = std::make_shared<threeParticleMassGenerator>(mD0, fsMasses, S12->kinSignature());
	std::shared_ptr<efficiencyFunction> efficiency  = std::make_shared<BELLE_DtoKpipi_efficiency>();
	std::vector<std::complex<double> > oldMin = utils::readComplexValuesFromTextFile("BELLE_sidebands_pi0_nsv_seed1533124364_ll-3894592.723834_N288410.143028.dat");

	std::shared_ptr<integrator> integral            = std::make_shared<integrator>(integralPoints, gen, amplitudes, efficiency);
	if (!integral->loadIntegrals("ps_" + integralFileName+"_regular.dat", "ac_" + integralFileName+"_regular.dat")) {
		std::cout << "Could not load integrals" << std::endl;
		return 1;
	}

	size_t rsz = 30;

	integral->resize(rsz);
	integral->setNonDiagToZero();
	integral->writeToFile("ps_hae.dat", false);
	integral->writeToFile("ac_hae.dat", true);

	std::vector<std::vector<std::complex<double> > > ps(3, std::vector<std::complex<double> >(3, std::complex<double>(0.,0.)));
	std::vector<std::vector<std::complex<double> > > ac(3, std::vector<std::complex<double> >(3, std::complex<double>(0.,0.)));
	ps[0][0] = 3.;
	ps[1][1] = 4.;
	ps[2][2] = 5.;

	ac[0][0] = 3.;
	ac[1][1] = 5.;
	ac[2][2] = 7.;

	ac[1][2] = std::complex<double>(3., 111.);
	ac[2][1] = std::complex<double>(3.,-111.);

	ac[0][1] = std::complex<double>(123., 37.);
	ac[1][0] = std::complex<double>(123.,-37.);;

	ac[2][1] = std::complex<double> (42., 41.);
	ac[1][2] = std::complex<double> (42.,-41.);


	integral->setIntegrals(ps,ac);
	std::vector<std::complex<double> > pa = oldMin;//  {std::complex<double>(11.,23.), std::complex<double>(13., 29.), std::complex<double> (17., 31.)};
//	pa[0] = std::complex<double>(0.,0.); pa[2] = std::complex<double>(0.,0.);
	pa.resize(rsz);

	double tI = integral->totalIntensity(pa, true);
	std::vector<std::complex<double> > paL = pa;
	std::vector<std::complex<double> > paR = pa;

	double LR =  integral->multiplyLR(paL, paR, true);

	std::vector<double> DL = integral->DLmultiplyLR(paL, paR, true);
	std::vector<double> DR = integral->DRmultiplyLR(paL, paR, true);

	std::vector<double> nDL = integral->DLmultiplyLR(paL, paR, true);
	std::vector<double> nDR = integral->DRmultiplyLR(paL, paR, true);

	double delta = 1.;
	for (size_t a = 0; a< paL.size(); ++a) {
		paL[a] += std::complex<double>(delta, 0.);
		nDL[2*a  ] = (integral->multiplyLR(paL, paR, true) - LR)/delta;
		paL[a] += std::complex<double>(-delta, delta);
		nDL[2*a+1] = (integral->multiplyLR(paL, paR, true) - LR)/delta;
		paL[a] += std::complex<double>(0., -delta);
		std::cout << nDL[2*a]/DL[2*a] << " LLL " << nDL[2*a+1]/DL[2*a+1] << std::endl;
		paR[a] += std::complex<double>(delta, 0.);
		nDR[2*a  ] = (integral->multiplyLR(paL, paR, true) - LR)/delta;
		paR[a] += std::complex<double>(-delta, delta);
		std::cout << "(" << integral->multiplyLR(paL, paR, true) << " - " << LR <<")/" << delta << " = ";
		nDR[2*a+1] = (integral->multiplyLR(paL, paR, true) - LR)/delta;
		std::cout << nDR[2*a+1] << std::endl;
		paR[a] += std::complex<double>(0., -delta);
		std::cout << nDR[2*a]/DR[2*a] << " RRR " << nDR[2*a+1]/DR[2*a+1] << std::endl;
	}

	return 0;

	std::vector<double> D(2*pa.size());
	std::vector<double> DD = integral->DtotalIntensity(pa, true);

//	double deltaSmall = 1.e-5;
	for (size_t i = 0; i < pa.size(); ++i) {
		pa[i] += std::complex<double>(delta, 0.);		
		D[2*i  ] = (integral->totalIntensity(pa, true)-tI)/delta;
		pa[i] += std::complex<double>(-delta, delta);
		D[2*i+1] = (integral->totalIntensity(pa, true)-tI)/delta;
		pa[i] += std::complex<double>(0., -delta);

//		pa[i] += std::complex<double>(deltaSmall, 0.);		
//		D[2*i  ] = (integral->totalIntensity(pa, true)-tI)/deltaSmall;
//		pa[i] += std::complex<double>(-deltaSmall, deltaSmall);
//		D[2*i+1] = (integral->totalIntensity(pa, true)-tI)/deltaSmall;
//		pa[i] += std::complex<double>(0., -deltaSmall);

		std::cout << D[2*i  ] / DD[2*i ] << std::endl <<  D[2*i+1] / DD[2*i+1] << std::endl;
	}
	return 0;


	bool doIntegration = false;

	if (doIntegration) {
		std::cout << "Starting integration" << std::endl;
		integral->integrate();
		std::cout << "Finished integration" << std::endl;

		integral->writeToFile("ps_"+integralFileName+".dat", false);
		integral->writeToFile("ac_"+integralFileName+".dat", true);

		std::cout << "Integral files written... finished" << std::endl;
		return 0;
	}
	if (!integral->loadIntegrals("ps_" + integralFileName+"_regular.dat", "ac_" + integralFileName+"_regular.dat")) {
		std::cout << "Could not load integrals" << std::endl;
		return 1;
	}

//	std::string BGfileName = "./BELLE_sidebands_pi0_nsv_seed1532622088_ll-3893334.857917.dat";
//	std::string BGfileName = "./BELLE_sidebands_pi0_nsv_seed1532955121_ll-3894864.545268_N288407.215450.dat";
//	std::vector<std::complex<double> > backgroundCoefficients = utils::readComplexValuesFromTextFile(BGfileName);

	std::vector<double> norms;
	{
		std::pair<bool, std::vector<double> > retNorm = integral->getNormalizations();
		if (!retNorm.first) {
			std::cout << "ERROR: Could not get normalizations" << std::endl;
			return 0;
		}
		norms = retNorm.second;
	}
	

//	std::shared_ptr<amplitude> bg_amplitude = std::make_shared<modelAmplitude> (backgroundCoefficients, amplitudes, norms, "backgroundParamneterization");
//	std::vector<std::shared_ptr<amplitude> > bg_amplitudes (1,bg_amplitude);
//	std::shared_ptr<integrator> integral_bg = std::make_shared<integrator>(integralPoints, gen, bg_amplitudes, efficiency);

//	integral_bg->integrate();
//	integral_bg->writeToFile("ps_bg_integral.dat", false);
//	integral_bg->writeToFile("ac_bg_integral.dat", true);
//	return 0;

//	if (!integral_bg->loadIntegrals("ps_bg_integral.dat","ac_bg_integral.dat")) {
//		std::cout << "Could not load background integral" << std::endl;
//		return 1;
//	}

//	std::cout << "Starting integration" << std::endl;
//	integral_bg->integrate();
//	std::cout << "Finished integration" << std::endl;

//	integral_bg->writeToFile("ps_bg_integral.dat", false);
//	integral_bg->writeToFile("ac_bg_integral.dat", true);

//	std::cout << "Integral files written... finished" << std::endl;

	bool signalEvents = false;

	std::string dataFileName;
	if (signalEvents) {
		dataFileName = "./BELLE_data.root";
	} else {
		dataFileName = "./BELLE_bothSidebands.root";
	}
	std::vector<std::vector<double> > dataPoints = getBELLEevents(dataFileName, softpionSign);

	const size_t nData = dataPoints.size();
	std::vector<std::complex<double> > startValues(nAmpl);


	// Only crude estimates for the resonance parameters to have nice start values with continuous phase motion

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
//	double width;

//	nData = 250.;

	std::complex<double> coeff(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpiS.size()-1;++b) {
		m2    = (binningKpiS[b] + binningKpiS[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_0, GK_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpiP.size()-1;++b) {
		m2 = (binningKpiP[b] + binningKpiP[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_1, GK_1);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpiD.size()-1;++b) {
		m2 = (binningKpiD[b] + binningKpiD[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_2, GK_2);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpiS.size()-1;++b) {
		m2 = (binningKpiS[b] + binningKpiS[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_0, GK_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpiP.size()-1;++b) {
		m2 = (binningKpiP[b] + binningKpiP[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_1, GK_1);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningKpiD.size()-1;++b) {
		m2 = (binningKpiD[b] + binningKpiD[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_2, GK_2);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningPiPiS.size()-1;++b) {
		m2 = (binningPiPiS[b] + binningPiPiS[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mf_0, Gf_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningPiPiP.size()-1;++b) {
		m2 = (binningPiPiP[b] + binningPiPiP[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mRho, Grho);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2()*nData);
	for (size_t b = 0; b < binningPiPiD.size()-1;++b) {
		m2 = (binningPiPiD[b] + binningPiPiD[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mf_2, Gf_2);
		++count;
	}
	if (count != nAmpl) {
		std::cout << "Number of initialized start values does not match " << count << " != " << nAmpl;
		return 1;
	}
//	outFile.open("startValueTest.deleteMe");
//	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
//	for (size_t a = 0; a < startValues.size(); ++a) {
//		outFile << a << " " << startValues[a].real() << " " << startValues[a].imag() << std::endl;
//	}
//	outFile.close();

	bool parallel = false;
	if (argc > 1) {
		parallel = true;
	}

	std::string outFileNameBase;
	if (signalEvents) {
		outFileNameBase = "./BELLEfit";
	} else {
		outFileNameBase = "./BELLE_sidebands";
	}

	if (parallel) {
		std::cout << "Running in fabili mode" << std::endl;
		std::shared_ptr<fabiliLL> ll = std::make_shared<fabiliLL>(amplitudes, integral);
//		std::shared_ptr<fabiliLL_openMP> ll = std::make_shared<fabiliLL_openMP>(amplitudes, integral);
		ll->loadDataPoints(dataPoints);
		ll->setExtended(true);

		fabili ff(ll);
		ff.setStepSize(.1);
		std::vector<double> ppp;
		for (std::complex<double>& a : startValues) {
			ppp.push_back(a.real());
			ppp.push_back(a.imag());
		}

		std::pair<bool, std::vector<double> >retVal = ff.minimize(ppp, 500);
		if (!retVal.first) {
			std::cout << "Minimization failed" << std::endl;
			return 1;
		}
		std::vector<std::complex<double> > fitResultProdAmps;
		for (size_t i = 0; i < retVal.second.size()/2; ++i) {
			fitResultProdAmps.push_back(std::complex<double>(retVal.second[2*i], retVal.second[2*i+1]));
		}
		std::cout << "Finished" << std::endl;
		evalType finalEval = ll->eval(retVal.second);
		std::cout << "Events in model: " << integral->totalIntensity(fitResultProdAmps,false) << " ::: " << integral->totalIntensity(fitResultProdAmps,true) << std::endl;

		std::string outFileName = outFileNameBase + std::string("_pi") + std::to_string(softpionSign) + std::string("_nsv_fabili_seed")+ std::to_string(seed) + std::string("_ll") + std::to_string(finalEval.value) + std::string(".dat");
		outFile.open(outFileName.c_str());
		outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
		for (std::complex<double>& amp : fitResultProdAmps) {
			outFile << amp << std::endl;
		}
		outFile.close();

		outFileName =  outFileNameBase + std::string("_hessian_") + std::to_string(seed) + std::string(".dat");
		outFile.open(outFileName.c_str());
		outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
		for (size_t i = 0; i < 2*ll->nAmpl(); ++i) {
			for (size_t j = 0; j < 2*ll->nAmpl(); ++j) {
				outFile << finalEval.hessian(i,j) << " ";
			}
			outFile << std::endl;
		}
		outFile.close();
	} else {
		std::cout << "Running in normal mode" << std::endl;

		logLikelihood ll(amplitudes, integral);
		ll.loadDataPoints(dataPoints);
		ll.setExtended(true);
//		ll.setNstore(10000);

		double minEval = ll.eval(oldMin);
		double numGradElement = 0.;
		double delta   = 1.e-6;
		std::vector<double> minGrad = ll.Deval(oldMin);
		std::vector<double> numGrad(minGrad.size());
		for (size_t a = 0; a < oldMin.size(); ++a) {
			oldMin[a] += std::complex<double>(delta,0.);
			numGradElement = (ll.eval(oldMin) - minEval)/delta;
			std::cout << numGradElement << " | | | " << minGrad[2*a  ] << std::endl;
			oldMin[a] += std::complex<double>(-delta,delta);
			numGradElement = (ll.eval(oldMin) - minEval)/delta;
			std::cout << numGradElement << " | | | " << minGrad[2*a+1] << std::endl;
			oldMin[a] += std::complex<double>(0.,-delta);
			numGradElement += numGradElement;
		}
		return 0;

		std::pair<double, std::vector<std::complex<double> > > retVal = ll.fitNlopt(startValues);


		std::string outFileName = outFileNameBase + std::string("_pi") + std::to_string(softpionSign) + std::string("_nsv_seed")+ std::to_string(seed) + std::string("_ll") + std::to_string(retVal.first) + std::string("_N")+std::to_string(integral->totalIntensity(retVal.second, true)) + std::string(".dat");
		outFile.open(outFileName.c_str());
		outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
		for (std::complex<double>& amp : retVal.second) {
			outFile << amp << std::endl;
		}
		outFile.close();
		std::vector<std::vector<double> > hessian = ll.DDeval(retVal.second);
		outFileName = outFileNameBase + std::string("_hessian_") + std::to_string(seed) + std::string(".dat");
		outFile.open(outFileName.c_str());
		outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
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
