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

#include"TH1D.h"
#include"TH2D.h"
#include"TFile.h"

std::vector<double> get_binning(const double m_min, const double m_max, const double bin_width, const double finer_min, const double finer_max, const int finer_factor) {
	double m = m_min;
	std::vector<double> binning = {m*m};
	while (m < m_max) {
		if (m > finer_min && m < finer_max) {
			m += bin_width/finer_factor;
		} else {
			m += bin_width;
		}
		binning.push_back(m*m);
	}
	return binning;
}

std::vector<std::vector<double > > get_all_binnings(const double mD0, const double mPi, const double mKs) {
	const double bin_width = 0.04;
	std::vector<std::vector<double> > retVal(9);
	retVal[0] = get_binning(  mPi + mKs, mD0 - mPi, bin_width, 20.   , 20.,    1); // No finer binning here (K pi S)
	retVal[1] = get_binning(  mPi + mKs, mD0 - mPi, bin_width,   .789,   .999, 4); // K pi P
	retVal[2] = get_binning(  mPi + mKs, mD0 - mPi, bin_width,  1.329,  1.529, 2); // K pi D
	retVal[3] = get_binning(  mPi + mKs, mD0 - mPi, bin_width, 20.   , 20.,    1); // No finer binning here (K pi S)
	retVal[4] = get_binning(  mPi + mKs, mD0 - mPi, bin_width,   .789,   .999, 4); // K pi P
	retVal[5] = get_binning(  mPi + mKs, mD0 - mPi, bin_width,  1.329,  1.529, 2); // K pi D
	retVal[6] = get_binning(2*mPi      , mD0 - mKs, bin_width,   .919,  1.079, 4); // pi pi S
	retVal[7] = get_binning(2*mPi      , mD0 - mKs, bin_width,   .679,   .919, 2); // pi pi P
	retVal[8] = get_binning(2*mPi      , mD0 - mKs, bin_width,  1.159,  1.39,  2); // pi pi D
	return retVal;
}

std::vector<std::shared_ptr<massShape> > get_isobars(const size_t waveIndex, const double mD0, const double mPi, const double mKs, bool writeHists = false) {
	std::vector<std::shared_ptr<massShape> > isobars;
	if (waveIndex > 8) {
		std::cout << "get_isobars(...): ERROR: Invalid wave index (>8) " << waveIndex << std::endl;
		return isobars;
	}
	double mRho      =  .77526; double Grho      = .1478;
//	double mRhoPrime = 1.465;   double GrhoPrime = .4; // This is never on-shell?? How can we treat this?
	double mf2       = 1.2751;  double Gf2       = .1851;
	double mOmega    =  .78265; double Gomega    = .00849;
//	double mK892     =  .89166; double GK892     = .0508;
	double mK892     =  .8937;  double GK892     = .0472;
//	double mK01430   = 1.425;   double GK01430   = .270;
	double mK21430   = 1.4256;  double GK21430   = .0985;
	double mK1410    = 1.414;   double GK1410    = .232;
	double mK1680    = 1.718;   double GK1680    = .322;
//	Ordering: "a", "r", "M0", "G0", "phiF", "phiR", "phiRsin", "F", "R", "MMax" // phiRsin = phiR ("Typo in some paper?? " - S. Wallner); chose MMax
//	std::vector<double> LASS_parameters    = {0.113, -33.8, 1.441, 0.193, utils::degToRad(0.1), utils::degToRad(-109.7), 0.96, 1.};
	std::vector<double> LASS_parameters    = {0.113, -33.797, 1.4405, 0.1926, utils::degToRad(0.099), utils::degToRad(-109.695), 0.955, 1.};
	std::vector<double> KMatrix_parameters = {  8.5212, utils::degToRad(  68.505), // As in the BELLE note
	                                           12.1895, utils::degToRad(  23.950), 
	                                           29.1463, utils::degToRad(  -0.106), 
	                                           10.7457, utils::degToRad( -51.893), 
	                                            0.    , utils::degToRad(   0.   ), 
	                                            8.0441, utils::degToRad(-125.963), 
	                                           26.2981, utils::degToRad(-152.324), 
	                                           33.0346, utils::degToRad( -93.228), 
	                                           26.1735, utils::degToRad(-121.405), 
	                                            0.    , utils::degToRad(   0.   ), 
	                                          -0.07};
//	std::vector<double> KMatrix_parameters = {  9.3 , utils::degToRad( -78.7), // As in the BarBar paper [https://arxiv.org/pdf/0804.2089.pdf]
//	                                           10.89, utils::degToRad(-159.1), 
//	                                           24.2 , utils::degToRad( 168.0), 
//	                                            9.16, utils::degToRad(  90.5), 
//	                                            0.  , utils::degToRad(   0. ), 
//	                                            7.94, utils::degToRad(  73.9), 
//	                                            2.0 , utils::degToRad( -18. ), 
//	                                            5.1 , utils::degToRad(  33. ), 
//	                                            3.23, utils::degToRad(   4.8), 
//	                                            0.  , utils::degToRad(   0. ), 
//	                                          -0.07};

	if (waveIndex == 0 || waveIndex == 3) { // Kpi S
		isobars.push_back(std::make_shared<BELLE_LASS>(LASS_parameters, mPi, mKs));
	} else if (waveIndex == 1 || waveIndex == 4) { // K pi P
		isobars.push_back(std::make_shared<BELLEbreitWigner>("K*(892)", mK892, GK892, 1, mD0, mPi, mKs, mPi));
		isobars.push_back(std::make_shared<BELLEbreitWigner>("K*(1410)", mK1410, GK1410, 1, mD0, mPi, mKs, mPi));
		if (waveIndex == 4) { // This is only used in one combination
			isobars.push_back(std::make_shared<BELLEbreitWigner>("K*(1680)", mK1680, GK1680, 1, mD0, mPi, mKs, mPi));
		}
	} else if (waveIndex == 2 || waveIndex == 5) { // K pi D
		isobars.push_back(std::make_shared<BELLEbreitWigner>("K2*(1430)", mK21430, GK21430, 2, mD0, mPi, mKs, mPi));
	} else if (waveIndex == 6) { // pi pi S
		isobars.push_back(std::make_shared<BELLE_KMatrix>(KMatrix_parameters)); // // // // // Replace me
	} else if (waveIndex == 7) { // pi pi P
		isobars.push_back(std::make_shared<BELLEbreitWigner>("omega", mOmega, Gomega, 1, mD0, mKs, mPi, mPi));
		isobars.push_back(std::make_shared<BELLEbreitWigner>("rho(770)", mRho, Grho, 1, mD0, mKs, mPi, mPi));
//		isobars.push_back(std::make_shared<BELLEbreitWigner>("rho'(1450)", mRhoPrime, GrhoPrime, 1, mD0, mKs, mPi, mPi));
	} else if (waveIndex == 8) { // pi pi D
		isobars.push_back(std::make_shared<BELLEbreitWigner>("f2(1270)", mf2, Gf2, 2, mD0, mKs, mPi, mPi));
	} else {
		std::cout << "get_isobars(...): ERROR: Invalid wave index " << waveIndex << std::endl;
		return isobars;
	}
	if (writeHists) {
		TFile* outFile = new TFile("isobars.root", "UPDATE");
		double mThresh = mPi + mKs;
		if (waveIndex > 5) {
			mThresh = 2*mPi;
		}
		for (std::shared_ptr<massShape> ms : isobars) {
			const size_t nBins = 1000;
			TH1D hist(ms->name().c_str(), ms->name().c_str(), nBins, 0., mD0);
			for (size_t b = 0; b < nBins; ++b) {
				double m = hist.GetXaxis()->GetBinCenter(b+1);
				if (m < mThresh) {
					continue;
				}
				hist.SetBinContent(b+1, std::norm(ms->eval(m*m)));
			}
			hist.Write();
		}
		outFile->Close();
	}
	return isobars;
}

std::vector<std::shared_ptr<amplitude> > get_model(std::vector<bool> free_map, const double mD0, const double mPi, const double mKs, bool cp_conj = false) {

	const std::vector<double> fs_masses = {mPi, mKs, mPi};

	std::vector<std::shared_ptr<amplitude> > model;

	const size_t piKrightIndex = cp_conj ? 23 : 12;
	const size_t piKwrongIndex = cp_conj ? 12 : 23;
	const size_t piPiIndex     = 13;

	std::shared_ptr<BELLE_S> S_right = std::make_shared<BELLE_S>(piKrightIndex);
	std::shared_ptr<BELLE_S> S_pipi  = std::make_shared<BELLE_S>(piPiIndex);
	std::shared_ptr<BELLE_S> S_wrong = std::make_shared<BELLE_S>(piKwrongIndex);

	std::shared_ptr<BELLE_P> P_right = std::make_shared<BELLE_P>(piKrightIndex);
	std::shared_ptr<BELLE_P> P_pipi  = std::make_shared<BELLE_P>(piPiIndex);
	std::shared_ptr<BELLE_P> P_wrong = std::make_shared<BELLE_P>(piKwrongIndex);

	std::shared_ptr<BELLE_D> D_right = std::make_shared<BELLE_D>(piKrightIndex);
	std::shared_ptr<BELLE_D> D_pipi  = std::make_shared<BELLE_D>(piPiIndex);
	std::shared_ptr<BELLE_D> D_wrong = std::make_shared<BELLE_D>(piKwrongIndex);

	std::vector<std::vector<double> > binnings = get_all_binnings(mD0, mPi, mKs);
	std::vector<double> binning;

	std::vector<size_t> isob_combinations = {piKrightIndex, piKrightIndex, piKrightIndex, piKwrongIndex, piKwrongIndex, piKwrongIndex, piPiIndex, piPiIndex, piPiIndex};
	std::vector<std::string> prefixes = {"KpiRight[","KpiRight[","KpiRight[","KpiWrong[","KpiWrong[","KpiWrong[", "piPi[", "piPi[", "piPi["};
	std::vector<std::string> suffixes = {"]PiS", "]PiP", "]PiD","]PiS", "]PiP", "]PiD", "]KS", "]KP", "]KD"};
	std::vector<std::shared_ptr<angularDependence> > angular_dependences = {S_right, P_right, D_right, S_wrong, P_wrong, D_wrong, S_pipi, P_pipi, D_pipi};
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for (size_t d = 0; d < 9; ++d) {
		if (free_map[d]) {
			binning = binnings[d];
			for (size_t b = 0; b < binning.size() - 1; ++b) {
				std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binning[b], binning[b+1]);
				std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> stepWave = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(isob_combinations[d], prefixes[d] + std::to_string(b) + suffixes[d], step, angular_dependences[d], fs_masses);
				model.push_back(stepWave);
			}
		} else {
			std::vector<std::shared_ptr<massShape> > isobars = get_isobars(d, mD0, mPi, mKs, true);
			for (std::shared_ptr<massShape>& isobar : isobars) {
				std::shared_ptr<threeParticlaIsobaricAmplitudeNoBose> wave = std::make_shared<threeParticlaIsobaricAmplitudeNoBose>(isob_combinations[d], prefixes[d] + isobar->name() + suffixes[d], isobar, angular_dependences[d], fs_masses);
				model.push_back(wave);
			}
		}
	}
	return model;
}

std::shared_ptr<amplitude> get_bg_amplitude() {
	std::shared_ptr<dalitzPolynomialAmplitude> bg_amplitude = std::make_shared<dalitzPolynomialAmplitude>(std::make_shared<kinematicSignature>(2), "dalitzBackground_version2.pol", .405, 1.845, .407, 2.057);
	return bg_amplitude;
}

int main(int argc, char* argv[]) {
	TH2D ful = TH2D("full_dalitz", "fullDalitz", 100, 0.,2.2, 100 ,0.3 , 3.3);
	size_t seed = size_t( time(NULL) );
	std::cout << "Dnine::main(...): INFO: Seed: " << seed << std::endl;
	srand(seed);

	const int softpionSign = 0;

	const double f_sig  = 6.75192e-01;// +- 9.11090e-04
	const double f_rand = 8.95779e-01;// +- 2.14394e-03
	const double f_CP   = 0.492;

	const double total_signal_coefficient = f_sig + (1.-f_sig)*f_rand*(1.-f_CP);
	const double total_CP_coefficient    = (1.-f_sig)*f_rand*f_CP;
	const double total_bg_coefficient    = (1.-f_sig) * (1.-f_rand);

//	std::vector<bool> free_map = {true, true, true, true, true, true, true, true, true};
	std::vector<bool> free_map = {false, false, false, false, false, false, false, false , false};
	std::vector<std::string> wave_names = {"KpiRightS", "KpiRightP", "KpiRightD", "KpiWrongS", "KpiWrongP", "KpiWrongD", "piPiS", "piPiP", "piPiD"};
	for (int i = 1; i < argc; ++i) {
		int waveIndex = atoi(argv[i]);
		if (waveIndex < 0 || waveIndex >= (int)free_map.size()) {
			std::cout << "Dnine::main(...): ERROR: Invalid wave index to free: " << waveIndex << std::endl;
			return 1;
		}
		std::cout << "Dnine::main(...): INFO: Freeing wave '" << wave_names[i] << "' (index " << waveIndex << ")" << std::endl;
		free_map[i] = true;
	}

	const size_t integral_points = 6000*10*10*10; // Merged two 30000000 integrals for the model // Comment for ease bugfix
	const std::vector<double> fs_masses = {mPi, mKs, mPi};
	std::vector<size_t> wave_numbers;
	std::vector<size_t> fix_amplitudes;
	size_t n_fix_S = 10;
	bool first_freed_S_wave = true;
	{
		std::vector<std::vector<double> > binnings = get_all_binnings(mD0, mPi, mKs);
		for (size_t d = 0; d < 9; ++d) {
			if ((d == 0 || d == 3 || d == 6) && free_map[d]) { // If it's an S-wave
				if (first_freed_S_wave) {
					first_freed_S_wave = false;
				} else {
					size_t cumul_n_waves = 0;
					for (size_t a : wave_numbers) {
						cumul_n_waves += a;
					}
					fix_amplitudes.push_back(cumul_n_waves + n_fix_S);
				}
			}
			if (free_map[d]) {
				wave_numbers.push_back(binnings[d].size() - 1);
			} else {
				wave_numbers.push_back(get_isobars(d, mD0, mPi, mKs).size());
			}
		}
	}

	std::shared_ptr<threeParticleMassGenerator> generator = std::make_shared<threeParticleMassGenerator>(mD0, fs_masses, std::make_shared<kinematicSignature>(2));
	std::shared_ptr<efficiencyFunction> efficiency        = std::make_shared<BELLE_DtoKpipi_efficiency>();

	std::vector<std::shared_ptr<amplitude> > model = get_model(free_map, mD0, mPi, mKs);
	size_t n_model = model.size();
	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integral_points, generator, model, efficiency);
	std::string integral_file_name = "integral_model_"; 
	for (bool free : free_map) {
		integral_file_name += std::to_string(free);
	}

	if (!integral->loadIntegrals("ps_"+integral_file_name+".dat","ac_"+integral_file_name+".dat")) {
		integral->integrate();
		integral->writeToFile("ps_"+integral_file_name+".dat", false);
                integral->writeToFile("ac_"+integral_file_name+".dat", true);
	}

	std::vector<std::shared_ptr<amplitude> > model_cp = get_model(free_map, mD0, mPi, mKs, true);
	std::shared_ptr<integrator> integral_cp = std::make_shared<integrator>(integral_points, generator, model, efficiency);
	std::string integral_cp_file_name =  "integral_cp_model_"; 
	for (bool free : free_map) {
		integral_cp_file_name += std::to_string(free);
	}

	if (!integral_cp->loadIntegrals("ps_"+integral_cp_file_name+".dat","ac_"+integral_file_name+".dat")) { // The CP conjugated integral is the same, since only 2 Kpi combinations are exchanged.
		integral_cp->integrate();
		integral_cp->writeToFile("ps_"+integral_cp_file_name+".dat", false);
		integral_cp->writeToFile("ac_"+integral_cp_file_name+".dat", true);
	}

	std::shared_ptr<amplitude> bg_amplitude = get_bg_amplitude();
	std::shared_ptr<integrator> integral_bg = std::make_shared<integrator>(integral_points, generator, std::vector<std::shared_ptr<amplitude> >(1,bg_amplitude), efficiency);

	std::string integral_bg_file_name = "integral_bg[" + bg_amplitude->name() + "]";

//	integral_bg->integrate();
//	integral_bg->writeToFile("ps_"+integral_bg_file_name+".dat", false);
//	integral_bg->writeToFile("ac_"+integral_bg_file_name+".dat", true);

	if (!integral_bg->loadIntegrals("ps_"+integral_bg_file_name+".dat","ac_"+integral_bg_file_name+".dat")) {
		std::cout << "Dnine::main(...): ERROR: Could not load bg integral" << std::endl;
		return 1;	
	}
	for (std::shared_ptr<amplitude> a : model_cp) {
		model.push_back(a);
	}
	model.push_back(bg_amplitude);
	if (!integral->addIncoherentSector(integral_cp)) {
		std::cout << "Dnine::main(...): ERROR: Could not add CP sector to integral" << std::endl;
		return 1;
	}
	if (!integral->addIncoherentSector(integral_bg)) {
		std::cout << "Dnine::main(...): ERROR: Could not add bg sector to integral" << std::endl;
		return 1;
	}

	bool signalEvents = true;
	std::string dataFileName;
	if (signalEvents) {
		dataFileName = "./BELLE_data.root";
	} else {
		dataFileName = "./BELLE_bothSidebandsHigherMD.root";
	}
	std::vector<std::vector<double> > dataPoints = getBELLEevents(dataFileName, softpionSign);
	dataPoints = utils::sanitizeDataPoints(dataPoints, fs_masses);
	const size_t nData = dataPoints.size();

	logLikelihood ll(model, integral);
	if (!ll.loadDataPoints(dataPoints)) {
		std::cout << "Dnine::main(...): ERROR: Could not load data points" << std::endl;
		return 1;
	}
	if (!ll.setExtended(true)) {
		std::cout << "Dnine::main(...): ERROR: Could not set extended" << std::endl;
		return 1;
	}
	if (!ll.fixParameter(1, 0.)) {
		std::cout << "Dnine::main(...): ERROR: Could not fix first phase" << std::endl;
		return 1;
	}
	
	for (size_t a : fix_amplitudes) {
		std::cout << "Dnine::main(...): INFO: fix parameter " << a << std::endl;
		if (!ll.fixParameter(2*a  ,0.)) {
			std::cout << "Dnine::main(...): ERROR: Could not fix parameter " << 2*a << " to 0." << std::endl;
			return 1;
		}
		if (!ll.fixParameter(2*a+1,0.)) {
			std::cout << "Dnine::main(...): ERROR: Could not fix parameter " << 2*a+1 << " to 0." << std::endl;
			return 1;
		}
	}
	double CP_scale_factor = pow(total_CP_coefficient/total_signal_coefficient, .5);
	for (size_t a = 0; a < n_model; ++a) {
		std::cout << "Copy " << n_model+a << " from " << a << std::endl;
		if (!ll.addCopyParameter(2*(n_model+a)  , 2*a  , -1, CP_scale_factor)) {
			std::cout << "Dnine::main(...): ERROR: Could not copy parameter " << 2*(n_model+a) << " from " << 2*a << std::endl;
			return 1;
		}
		if (!ll.addCopyParameter(2*(n_model+a)+1, 2*a+1, -1, CP_scale_factor)) {
			std::cout << "Dnine::main(...): ERROR: Could not copy parameter " << 2*(n_model+a)+1 << " from " << 2*a+1 << std::endl;
			return 1;
		}
	}
	std::cout << "Dnine::main(...): INFO: Fix BG parameter " << model.size()-1 << std::endl;
	if (!ll.fixParameter(2*model.size()-2, pow(nData * total_bg_coefficient, .5))) { // Fix imag part of the BG amplitude to deternmined BG ratio
		std::cout << "Dnine::main(...): ERROR: Could not fix BG amplitude's real part (" << 2*model.size()-2 << ") to determiend value" << std::endl;
		return 1;
	}
	if (!ll.fixParameter(2*model.size()-1, 0.)) { // Fix imag part of the BG amplitude to zero
		std::cout << "Dnine::main(...): ERROR: Could not fix BG amplitude's imag part (" << 2*model.size()-1 << ") to zero" << std::endl;
		return 1;
	}
	std::vector<std::complex<double> > startValues;
	for (size_t  a = 0; a < model.size(); ++a) {
		startValues.push_back(std::complex<double>(utils::random2(), utils::random2()));
	}
	double factor = pow(integral->totalIntensity(startValues)/dataPoints.size(), .5);
	for (size_t a = 0; a < model.size(); ++a) {
		startValues[a] /= factor;
	}
	std::pair<double, std::vector<std::complex<double> > > retVal = ll.fitNlopt(startValues);
	
	std::string outFileName = "BELLE_fit_";
	for (bool free : free_map) {
		outFileName += std::to_string(free);
	}
	outFileName += "_"+std::to_string(seed)+".dat";
	std::ofstream outFile;
	outFile.open(outFileName.c_str());
	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (std::complex<double>& amp : retVal.second) {
		outFile << amp << std::endl;
	}
	outFile.close();
	return 0;
}
