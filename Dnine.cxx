#include"Dnine.h"

int main(int argc, char* argv[]) {
	TH2D ful = TH2D("full_dalitz", "fullDalitz", 100, 0.,2.2, 100 ,0.3 , 3.3);
	const size_t seed = size_t( time(NULL) );
	std::cout << "Dnine::main(...): INFO: Seed: " << seed << std::endl;
	srand(seed);

	const int softpionSign = 0;

	const double f_sig  = 6.75192e-01;// +- 9.11090e-04
	const double f_rand = 8.95779e-01;// +- 2.14394e-03
	const double f_CP   = 0.492;

	const double total_signal_coefficient = f_sig + (1.-f_sig)*f_rand*(1.-f_CP);
	const double total_CP_coefficient    = (1.-f_sig)*f_rand*f_CP;
	const double total_bg_coefficient    = (1.-f_sig)*(1.-f_rand);

//	std::cout << total_signal_coefficient << " " << total_CP_coefficient << " " << total_bg_coefficient << std::endl;
//	return 0;

//	std::vector<bool> free_map = {true, true, true, true, true, true, true, true, true};
	std::vector<bool> free_map = {false, false, false, false, false, false, false, false , false};
	std::vector<std::string> wave_names = {"KpiRightS", "KpiRightP", "KpiRightD", "KpiWrongS", "KpiWrongP", "KpiWrongD", "piPiS", "piPiP", "piPiD"};
	for (int i = 1; i < argc; ++i) {
		int waveIndex = atoi(argv[i]);
		if (waveIndex < 0 || waveIndex >= (int)free_map.size()) {
			std::cout << "Dnine::main(...): ERROR: Invalid wave index to free: " << waveIndex << std::endl;
			return 1;
		}
		std::cout << "Dnine::main(...): INFO: Freeing wave '" << wave_names[waveIndex] << "' (index " << waveIndex << ")" << std::endl;
		free_map[waveIndex] = true;
	}
	std::string freeString = "";
	for (bool free : free_map) {
		freeString += std::to_string(free);
	}

	const size_t integral_points = 6000*10*10*10*10; // Merged two 30000000 integrals for the model // Comment for ease bugfix
	const std::vector<double> fs_masses = {mPi, mKs, mPi};
	std::vector<size_t> wave_numbers;
	std::vector<size_t> fix_amplitudes;
	const size_t n_fix_S = 10;
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
	std::vector<std::shared_ptr<amplitude> > fixed_model = get_model({false, false, false, false, false, false, false, false, false}, mD0, mPi, mKs);
	const size_t n_model = model.size();
	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integral_points, generator, model, efficiency);
	std::string integral_file_name = "integral_model_" + freeString; 

	std::shared_ptr<integrator> fixed_integral = std::make_shared<integrator>(integral_points, generator, fixed_model, efficiency);
	if (!fixed_integral->loadIntegrals("./integralFiles/ps_integral_model_000000000.dat","./integralFiles/ac_integral_model_000000000.dat")) {
		std::cout << "Dnine::main(...): ERROR: Could not load fixed model integrals" << std::endl;
		return 1;
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: './integralFiles/ps_integral_model_000000000.dat' and './integralFiles/ac_integral_model_000000000.dat'" <<std::endl;
	}

	if (!integral->loadIntegrals("./integralFiles/ps_"+integral_file_name+".dat","./integralFiles/ac_"+integral_file_name+".dat")) {
		integral->integrate();
		integral->writeToFile("./integralFiles/ps_"+integral_file_name+".dat", false);
                integral->writeToFile("./integralFiles/ac_"+integral_file_name+".dat", true);
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: " << "'./integralFiles/ps_"+integral_file_name+".dat' and './integralFiles/ac_"+integral_file_name+".dat'" <<std::endl;
	}

	std::vector<std::shared_ptr<amplitude> > model_cp = get_model(free_map, mD0, mPi, mKs, true);
	std::shared_ptr<integrator> integral_cp = std::make_shared<integrator>(integral_points, generator, model, efficiency);
	std::string integral_cp_file_name =  "integral_cp_model_" + freeString; 

	if (!integral_cp->loadIntegrals("./integralFiles/ps_"+integral_cp_file_name+".dat","./integralFiles/ac_"+integral_cp_file_name+".dat")) {
		integral_cp->integrate();
		integral_cp->writeToFile("./integralFiles/ps_"+integral_cp_file_name+".dat", false);
		integral_cp->writeToFile("./integralFiles/ac_"+integral_cp_file_name+".dat", true);
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: " << "'./integralFiles/ps_"+integral_cp_file_name+".dat' and './integralFiles/ac_"+integral_cp_file_name+".dat'" <<std::endl;
	}

	std::shared_ptr<amplitude> bg_amplitude = get_bg_amplitude();
	std::shared_ptr<integrator> integral_bg = std::make_shared<integrator>(integral_points, generator, std::vector<std::shared_ptr<amplitude> >(1,bg_amplitude), efficiency);

	std::string integral_bg_file_name = "integral_bg[" + bg_amplitude->name() + "]";

	if (!integral_bg->loadIntegrals("./integralFiles/ps_"+integral_bg_file_name+".dat","./integralFiles/ps_"+integral_bg_file_name+".dat")) { // Use two times ps, since the BG parameterization is already acceptance correctezd...
		std::cout << "Dnine::main(...): WARNING: Could not load bg integral. Integrate" << std::endl;
		integral_bg->integrate();
		integral_bg->writeToFile("./integralFiles/ps_"+integral_bg_file_name+".dat",false);
		integral_bg->writeToFile("./integralFiles/ac_"+integral_bg_file_name+".dat",true);
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: " << "'./integralFiles/ps_"+integral_bg_file_name+".dat' and './integralFiles/ac_"+integral_bg_file_name+".dat'" <<std::endl;
	}

	const bool signalEvents = true;
	std::string dataFileName;
	if (signalEvents) {
		dataFileName = "./BELLE_data.root";
	} else {
		dataFileName = "./BELLE_bothSidebandsHigherMD.root";
	}
	std::vector<std::vector<double> > dataPoints = getBELLEevents(dataFileName, softpionSign);
	dataPoints = utils::sanitizeBELLEdataPoints(dataPoints, fs_masses);
	const size_t nData = dataPoints.size();

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
	
//	for (size_t a = 0; a < 13; ++a) {
//		fix_amplitudes.push_back(28+a);
//	}

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
	std::vector<size_t> fixToZero;
//	fixToZero = {14, 15}; // Remove the K*(1680)

	for (size_t z : fixToZero) {
		if (!ll.fixParameter(z, 0.)) {
			std::cout << "Dnine::main(...): ERROR: Could not fix parameter " << z << " to 0." << std::endl;
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
	const double factor = pow(integral->totalIntensity(startValues)/dataPoints.size(), .5);
	for (size_t a = 0; a < model.size(); ++a) {
		startValues[a] /= factor;
	}
//	std::pair<double, std::vector<std::complex<double> > > retVal = ll.fitNlopt(startValues);
	std::vector<std::complex<double> > BELLE_result_values = getStartValuesForBelleParameters(free_map, integral->getNormalizations(false).second, fixed_integral->getNormalizations(false).second);
	for (size_t a = 0; a < n_model; ++a) {
			BELLE_result_values.push_back(CP_scale_factor * BELLE_result_values[a]);
	}
	BELLE_result_values.push_back(std::complex<double>(pow(nData * total_bg_coefficient, .5),0.));
	const double BELLE_LL = ll.eval(BELLE_result_values);
	std::pair<double, std::vector<std::complex<double> > > retVal(BELLE_LL, BELLE_result_values);
	
	std::string outFileName = "./BELLE_fit_results/BELLE_fit_" + freeString;
	outFileName += "_"+std::to_string(retVal.first)+"_"+std::to_string(seed)+".dat";
	std::ofstream outFile;
	outFile.open(outFileName.c_str());
	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (std::complex<double>& amp : retVal.second) {
		outFile << amp << std::endl;
	}
	outFile.close();

	bool makeAmplitudeDalitzPlots = utils::updateBestLLfile("bestLL_"+freeString, retVal.first, std::to_string(seed));
	if (makeAmplitudeDalitzPlots) {
		std::string outFileNameDalitz = "dalitzPlotResults_"+freeString	+".root";
		std::pair<bool, std::vector<double> > norms_ac = integral->getNormalizations(true);
		std::pair<bool, std::vector<double> > norms    = integral->getNormalizations(false);
		if (!norms.first or !norms_ac.first) {
			std::cout << "Dnine::main(...): ERROR: Could not get normalizations" << std::endl;
			return 1;
		}
		for (size_t a = 0; a < norms.second.size(); ++a) {
//			norms.second[a] = 1.;//norms_ac.second[a];///norms.second[a];
		}
		TFile* outFile = new TFile(outFileNameDalitz.c_str(), "RECREATE");
		TH2D dataHist  = makeDataDalitz(dataPoints, fs_masses);
		dataHist.Write();
		TH2D dataHistOrtho  = makeDataDalitz(dataPoints, fs_masses, 12,23);
		dataHistOrtho.Write();
		const bool singleWavePlots = false;
		if (singleWavePlots) {
			for (size_t a = 0; a < 2*n_model; ++a) {
				TH2D hist = makeDalizFromModel({retVal.second[a]}, {model[a]}, {norms.second[a]}, fs_masses, efficiency);
	//			TH2D hist = makeDalizFromModel({retVal.second[a]}, {model[a]}, fs_masses, efficiency);
				hist.Write();
				TH2D histOrtho = makeDalizFromModel({retVal.second[a]}, {model[a]}, {norms.second[a]}, fs_masses, efficiency, 12,23);
	//			TH2D histOrtho = makeDalizFromModel({retVal.second[a]}, {model[a]}, fs_masses, efficiency, 12,23);
				histOrtho.Write();
			}
		}
		TH2D histAll = makeDalizFromModel(std::vector<std::complex<double> >(&retVal.second[0], &retVal.second[n_model]), 
		                                  std::vector<std::shared_ptr<amplitude> >(&model[0], &model[n_model]), 
		                                  std::vector<double>(&norms.second[0], &norms.second[n_model]), 
		fs_masses, efficiency);
		histAll.Write();
		TH2D histAllOrtho = makeDalizFromModel(std::vector<std::complex<double> >(&retVal.second[0], &retVal.second[n_model]), 
		                                       std::vector<std::shared_ptr<amplitude> >(&model[0], &model[n_model]), 
		                                       std::vector<double>(&norms.second[0], &norms.second[n_model]), 
		fs_masses, efficiency, 12, 23);
		histAllOrtho.Write();
		TH2D histAllCP = makeDalizFromModel(std::vector<std::complex<double> >(&retVal.second[n_model], &retVal.second[2*n_model]), 
		                                    std::vector<std::shared_ptr<amplitude> >(&model[n_model], &model[2*n_model]), 
		                                    std::vector<double>(&norms.second[n_model], &norms.second[2*n_model]),
		fs_masses, efficiency, 12, 13, "_CP");
		histAllCP.Write();
		TH2D histAllOrthoCP = makeDalizFromModel(std::vector<std::complex<double> >(&retVal.second[n_model], &retVal.second[2*n_model]), 
		                                         std::vector<std::shared_ptr<amplitude> >(&model[n_model], &model[2*n_model]), 
		                                         std::vector<double>(&norms.second[n_model], &norms.second[2*n_model]),
		fs_masses, efficiency, 12, 23, "_CP");
		histAllOrthoCP.Write();

		std::shared_ptr<efficiencyFunction> unitEfficiency = std::make_shared<threeParticlPerfectEfficiency>(std::make_shared<kinematicSignature>(2));
		TH2D histBg = makeDalizFromModel({retVal.second[2*n_model]}, {model[2*n_model]}, {norms.second[2*n_model]}, fs_masses, unitEfficiency);
//		TH2D histBg = makeDalizFromModel({retVal.second[2*n_model]}, {model[2*n_model]}, fs_masses, unitEfficiency);
		histBg.Write();
		TH2D histBgOrtho = makeDalizFromModel({retVal.second[2*n_model]}, {model[2*n_model]}, {norms.second[2*n_model]}, fs_masses, unitEfficiency, 12, 23);
//		TH2D histBgOrtho = makeDalizFromModel({retVal.second[2*n_model]}, {model[2*n_model]}, fs_masses, unitEfficiency, 12, 23);
		histBgOrtho.Write();
		outFile->Close();
		std::cout << "Dnine::main(...): INFO: Finished creating Dalitz plots... exit" << std::endl;
		return 0;
	}
	return 0;
}
