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

#include"TH2D.h"
#include"TFile.h"

std::complex<double> BW(const double& m2, const double&  m0, const double& G0) {
	return std::complex<double>(m0*G0,0.)/std::complex<double>(m0*m0 - m2, -m0*G0);
}

int main(int argc, char *argv[]) {
	(void) argv;
	TH2D ful = TH2D("full_dalitz", "fullDalitz", 100, 0.,2.2, 100 ,0.3 , 3.3);
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

	double m = mPi + mKs;
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
	size_t fixNbinSwaves = 10;
	std::vector<size_t> fixAmplitudes;
	std::vector<std::shared_ptr<amplitude> > amplitudes = {};
// - - - - - - - K piRight S-wave
	for (size_t b = 0; b < binningKpiS.size() - 1; ++b) {
		std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningKpiS[b], binningKpiS[b+1]);
		std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave  = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(12, std::string("KpiRight[") + std::to_string(b) + std::string("]PiS"), step, S12, fsMasses);
		if (b == fixNbinSwaves) {
			fixAmplitudes.push_back(amplitudes.size()); // do this befre push back, so the index is right (Otherwirse .size() - 1)
		}
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
		if (b == fixNbinSwaves) {
			fixAmplitudes.push_back(amplitudes.size()); // do this befre push back, so the index is right (Otherwirse .size() - 1)
		}
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
		if (b == fixNbinSwaves) {
			fixAmplitudes.push_back(amplitudes.size()); // do this befre push back, so the index is right (Otherwirse .size() - 1)
		}
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

	std::shared_ptr<integrator> integral            = std::make_shared<integrator>(integralPoints, gen, amplitudes, efficiency);
	if (!integral->loadIntegrals("ps_" + integralFileName+"_regular.dat", "ac_" + integralFileName+"_regular.dat")) {
		std::cout << "Could not load integrals" << std::endl;
		return 1;
	}

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

//	std::shared_ptr<mCosTintensPolynomial> bg_amplitude = std::make_shared<mCosTintensPolynomial>(std::make_shared<kinematicSignature>(2), "dalitzBackground_version1.pol", 1.905, fsMasses, 13);
	std::shared_ptr<dalitzPolynomialAmplitude> bg_amplitude = std::make_shared<dalitzPolynomialAmplitude>(std::make_shared<kinematicSignature>(2), "dalitzBackground_version2.pol", .405, 1.845, .407, 2.057);
//	bg_amplitude->setMlimits(std::pair<double,double>(0.26,1.44913767462));
	std::vector<std::shared_ptr<amplitude> > bg_amplitudes = {bg_amplitude};
//	for (int i = 0; i < ful.GetNbinsX(); ++i) {
//		double x = ful.GetXaxis()->GetBinCenter(i+1);
//		for (int j = 0; j < ful.GetNbinsX(); ++j) {
//			double y = ful.GetYaxis()->GetBinCenter(j+1);
//			std::vector<double> kkiinn = {mD0*mD0, y,x};
//			if (gen->isValidPoint(kkiinn)) {
//				std::complex<double> val = bg_amplitude->eval(kkiinn);
//				ful.SetBinContent(i+1, j+1, norm(val));
//			}
//		}
//	}
//	TFile* ouFile = new TFile("lloll.root","RECREATE");
//	ful.Write();
//	ouFile->Close();
//	return 0;

	std::shared_ptr<integrator> integral_bg = std::make_shared<integrator>(integralPoints, gen, bg_amplitudes, efficiency);

//	integral_bg->integrate();
//	integral_bg->writeToFile("ps_bg_integral.dat", false);
//	integral_bg->writeToFile("ac_bg_integral.dat", true);
//	return 0;

	if (!integral_bg->loadIntegrals("ps_integral_dalitzBackground_version1.dat","ac_integral_dalitzBackground_version1.dat")) {
		std::cout << "Could not load background integral" << std::endl;
		return 1;
	}

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
		dataFileName = "./BELLE_bothSidebandsHigherMD.root";
	}
	std::vector<std::vector<double> > dataPoints = getBELLEevents(dataFileName, softpionSign);
//	dataPoints.resize(1000);

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

	std::complex<double> coeff(utils::random2()*nData,utils::random2());
	for (size_t b = 0; b < binningKpiS.size()-1;++b) {
		m2    = (binningKpiS[b] + binningKpiS[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_0, GK_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2());
	for (size_t b = 0; b < binningKpiP.size()-1;++b) {
		m2 = (binningKpiP[b] + binningKpiP[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_1, GK_1);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2());
	for (size_t b = 0; b < binningKpiD.size()-1;++b) {
		m2 = (binningKpiD[b] + binningKpiD[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_2, GK_2);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2());
	for (size_t b = 0; b < binningKpiS.size()-1;++b) {
		m2 = (binningKpiS[b] + binningKpiS[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_0, GK_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2());
	for (size_t b = 0; b < binningKpiP.size()-1;++b) {
		m2 = (binningKpiP[b] + binningKpiP[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_1, GK_1);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2());
	for (size_t b = 0; b < binningKpiD.size()-1;++b) {
		m2 = (binningKpiD[b] + binningKpiD[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mK_2, GK_2);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2());
	for (size_t b = 0; b < binningPiPiS.size()-1;++b) {
		m2 = (binningPiPiS[b] + binningPiPiS[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mf_0, Gf_0);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2());
	for (size_t b = 0; b < binningPiPiP.size()-1;++b) {
		m2 = (binningPiPiP[b] + binningPiPiP[b+1])/2;
//		width = binningKpiS[b+1] - binningKpiS[b];
		startValues[count] = coeff*BW(m2, mRho, Grho);
		++count;
	}
	coeff = std::complex<double>(utils::random2()*nData,utils::random2());
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
	double scaller = integral->totalIntensity(startValues)/nData/20;
	scaller = pow(1./scaller,.5);

	scaller = 0.00001;

	for (size_t a = 0; a < startValues.size(); ++a) {
		startValues[a] *= scaller;
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

//	startValues = utils::readComplexValuesFromTextFile("BELLE_sidebands_pi0_nsv_seed1533231138_ll-3891379.772211_N288558.816340.dat");

	amplitudes.push_back(bg_amplitude);
	integral->addIncoherentSector(integral_bg);
	double regionAreaFactor = .5;//0.15625; // Signal region is smaller than sidemabnd region by this cator (2 x 0.03) vs. (9.6 x 2*(0.02))
	startValues.push_back(std::complex<double>(pow(288661.997726*regionAreaFactor,.5),0.));

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
		ll.setExtended(true  );
//		ll.setNstore(10000   );
//		ll.setMaxII(10.*nData);

		for (size_t a : fixAmplitudes) {
			ll.fixParameter(2*a  ,0.);
			ll.fixParameter(2*a+1,0.);
			startValues[a] = std::complex<double>(0.,0.);
		}
		ll.addCopyParameter(22, 3);
		ll.addCopyParameter(23, 4, 0, 10000.);
		ll.addCopyParameter(25, 6, 0, 10000.);
		ll.addCopyParameter(100, 80, 1, 500.);

		std::cout << "Start fitting" << std::endl;
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
