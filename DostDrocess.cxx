#include "Dnine.h"
#include <sstream>

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cout << "DostDrocess.cxx::main(): ERROR: No input file given" << std::endl;
		return 1;
	}

	const int softpionSign = 0;

	std::string resultFileName = argv[1];

	const size_t      integral_points = 6000*10*10*10*10; // Merged two 30000000 integrals for the model // Comment for ease bugfix
	double            inputLL         = std::numeric_limits<double>::infinity();
	bool              prior           = false;
	bool              copy            = false;
	std::vector<bool> freeMap         = {false, false, false, false, false, false, false, false, false};
	std::vector<bool> fixToZeroMap    = {false, false, false, false, false, false, false, false, false};

	const double f_sig                    = 6.75192e-01;// +- 9.11090e-04
	const double f_rand                   = 8.95779e-01;// +- 2.14394e-03
	const double f_CP                     = 0.492;
	const double total_signal_coefficient = f_sig + (1.-f_sig)*f_rand*(1.-f_CP);
	const double total_CP_coefficient     = (1.-f_sig)*(1.-f_rand)*f_CP;
//	const double total_bg_coefficient     = (1.-f_sig)*f_rand;
	const double CP_scale_factor          = pow(total_CP_coefficient/total_signal_coefficient, .5);

	std::string freeString;
	std::vector<std::complex<double> > prodAmps = utils::readComplexValuesFromTextFile(resultFileName, false);

	const std::vector<double>                   fs_masses  = {mPi, mKs, mPi};
	std::shared_ptr<threeParticleMassGenerator> generator  = std::make_shared<threeParticleMassGenerator>(mD0, fs_masses, std::make_shared<kinematicSignature>(2));
	std::shared_ptr<efficiencyFunction>         efficiency = std::make_shared<BELLE_DtoKpipi_efficiency>();

	{
		std::vector<std::string> splitted = utils::splitString(resultFileName, '/');
		std::vector<std::string> parts    = utils::splitString(splitted[splitted.size()-1], '_');

		freeString = parts[2];
		inputLL    = atof(parts[3].c_str());

		std::string infoFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/BELLE_fit_results/fitInfos/BELLE_fit_";
		infoFileName += parts[2] + "_" + parts[4];

		std::string line;
		std::ifstream infoFile(infoFileName.c_str());
		while (std::getline(infoFile, line)) {
			if (line == "./Dnine.exe ") {
				continue;
			}
			if (line == "-prior ") {
				std::cout << "DostDrocess.cxx::main(): INFO: Using prior" << std::endl;
				prior = true;
			} else if(line == "-copy ") {
				std::cout << "DostDrocess.cxx::main(): INFO: Copy KpiRight to KpiWrong" << std::endl;
				copy = true;
			} else {
				int i = atoi(line.c_str());
				if (i == 0 and line != "0 ") {
					std::cout << "DostDrocess.cxx::main(): ERROR: Invalid option: '" << line << "'" << std::endl;
					return 1;
				}
				if (i < 0) {
					int index = -i-1;
					std::cout << "DostDrocess::main(...): INFO: Fixing wave '" << wave_names[index] << "' (index " << index << ") to zero" << std::endl;
					fixToZeroMap[index] = true;
				} else {
					freeMap[i] = true;
					std::cout << "DostDrocess::main(...): INFO: Freeing wave '" << wave_names[i] << "' (index " << i << ")" << std::endl;
				}
			}
		}
	}

	const std::vector<size_t> n_waves = getNwaves(freeMap);

	std::vector<std::shared_ptr<amplitude> > model    = get_model(freeMap, mD0, mPi, mKs);
	std::vector<std::shared_ptr<amplitude> > model_cp = get_model(freeMap, mD0, mPi, mKs, true);
	std::shared_ptr<amplitude> bg_amplitude           = get_bg_amplitude(fs_masses);

	const size_t n_model = model.size();

	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integral_points, generator, model, efficiency);
	std::shared_ptr<integrator> integral_cp = std::make_shared<integrator>(integral_points, generator, model_cp, efficiency);
	std::shared_ptr<integrator> integral_bg = std::make_shared<integrator>(integral_points, generator, std::vector<std::shared_ptr<amplitude> >(1,bg_amplitude), efficiency);

	{
		std::string integral_file_name = "integral_model_" + freeString;
		if (integral->loadIntegrals("./integralFiles/ps_"+integral_file_name+"_regular." + branchFileEnding,"./integralFiles/ac_"+integral_file_name+"_regular." + branchFileEnding)) {
			std::cout << "DostDrocess::main(...): INFO: Integral loaded" << std::endl;
		} else {
			std::cout << "DostDrocess::main(...): ERROR: Could not load integral" << std::endl;
			return 1;
		}

		std::string integral_cp_file_name =  "integral_cp_model_" + freeString;
		if (integral_cp->loadIntegrals("./integralFiles/ps_"+integral_cp_file_name+"_regular." + branchFileEnding,"./integralFiles/ac_"+integral_cp_file_name+"_regular." + branchFileEnding)) {
			std::cout << "DostDrocess::main(...): INFO: CP integral loaded" << std::endl;
		} else {
			std::cout << "DostDrocess::main(...): ERROR: Could not load CP integral" << std::endl;
			return 1;
		}

		std::string integral_bg_file_name = "integral_bg[" + bg_amplitude->name() + "]";
		if (integral_bg->loadIntegrals("./integralFiles/ps_"+integral_bg_file_name+"_regular." + branchFileEnding,"./integralFiles/ac_"+integral_bg_file_name+"_regular." + branchFileEnding)) {
			std::cout << "DostDrocess::main(...): INFO: BG integral loaded" << std::endl;
		} else {
			std::cout << "DostDrocess::main(...): ERROR: Could not load BG integral" << std::endl;
			return 1;
		}

	}

	for (std::shared_ptr<amplitude> a : model_cp) {
		model.push_back(a);
	}

	model.push_back(bg_amplitude);
	if (!integral->addIncoherentSector(integral_cp)) {
		std::cout << "DostDrocess::main(...): ERROR: Could not add CP sector to integral" << std::endl;
		return 1;
	}
	if (!integral->addIncoherentSector(integral_bg)) {
		std::cout << "DostDrocess::main(...): ERROR: Could not add bg sector to integral" << std::endl;
		return 1;
	}
	std::shared_ptr<logLikelihood_withPrior> ll = std::make_shared<logLikelihood_withPrior>(model, integral);

	const bool signalEvents = true;
	std::string dataFileName;
	if (signalEvents) {
		dataFileName = "./BELLE_data.root";
	} else {
		dataFileName = "./BELLE_bothSidebandsHigherMD.root";
	}
	std::vector<std::vector<double> > dataPoints = getBELLEevents(dataFileName, softpionSign);
	dataPoints = utils::sanitizeBELLEdataPoints(dataPoints, fs_masses);
	if (!ll->loadDataPoints(dataPoints, 25)) {
		std::cout << "DostDrocess::main(...): ERROR: Could not load data points" << std::endl;
		return 1;
	}
	std::cout << "DostDrocess::main(...): INFO: Finished loading data points to log-likelihood." << std::endl;

	if (!doTheFixingAndCopying(ll, n_waves, CP_scale_factor, copy, fixToZeroMap, true)) {
		std::cout << "DostDrocess::main(...): ERROR: Could not do the fixing and copying. Abort..." << std::endl;
		return 1;
	}
	if (prior) {
		ll-> setInterferencePriorStrength(1000000.);
	}

	double evalLL =  ll->eval(prodAmps);
	std::cout << "DostDrocess::main(...): INFO: Likelihood eval: " << evalLL << " difference to input is: " << inputLL-evalLL << std::endl;

	std::shared_ptr<efficiencyFunction> unitEfficiency = std::make_shared<threeParticlPerfectEfficiency>(std::make_shared<kinematicSignature>(2));
	std::string outFileNameDalitz = "dalitzPlotResults_"+freeString + "_" + branchFileEnding +".root";
	std::pair<bool, std::vector<double> > norms_ac = integral->getNormalizations(true);
	std::pair<bool, std::vector<double> > norms    = integral->getNormalizations(false);
	if (!norms.first or !norms_ac.first) {
		std::cout << "DostDrocess::main(...): ERROR: Could not get normalizations" << std::endl;
		return 1;
	}
//	for (size_t a = 0; a < norms.second.size(); ++a) {
//		norms.second[a] = 1./norms.second[a];//orms_ac.second[a];///norms.second[a];
//	}
	TFile* outFile = new TFile(outFileNameDalitz.c_str(), "RECREATE");
	TH2D dataHist  = makeDataDalitz(dataPoints, fs_masses);
	dataHist.Write();
	TH2D dataHistOrtho  = makeDataDalitz(dataPoints, fs_masses, 12,23);
	dataHistOrtho.Write();
	const bool singleWavePlots = false;
	if (singleWavePlots) {
		for (size_t a = 0; a < 2*n_model; ++a) {
			TH2D hist = makeDalizFromModel({prodAmps[a]}, {model[a]}, {norms.second[a]}, fs_masses, efficiency);
//			TH2D hist = makeDalizFromModel({prodAmps[a]}, {model[a]}, fs_masses, efficiency);
			hist.Write();
			TH2D histOrtho = makeDalizFromModel({prodAmps[a]}, {model[a]}, {norms.second[a]}, fs_masses, efficiency, 12,23);
//			TH2D histOrtho = makeDalizFromModel({prodAmps[a]}, {model[a]}, fs_masses, efficiency, 12,23);
			histOrtho.Write();
		}
	}
	size_t countWaves = 0;
	for (size_t w = 0; w < 9; ++w) {
		TH2D histWave = makeDalizFromModel(std::vector<std::complex<double> >(&prodAmps[countWaves], &prodAmps[countWaves + n_waves[w]]),
		                                   std::vector<std::shared_ptr<amplitude> >(&model[countWaves], &model[countWaves + n_waves[w]]),
		                                   std::vector<double>(&norms.second[countWaves], &norms.second[countWaves + n_waves[w]]),
		fs_masses, efficiency, 12, 13, "_"+wave_names[w]);
		histWave.Write();
		TH2D histWaveOrtho = makeDalizFromModel(std::vector<std::complex<double> >(&prodAmps[countWaves], &prodAmps[countWaves + n_waves[w]]),
		                                        std::vector<std::shared_ptr<amplitude> >(&model[countWaves], &model[countWaves + n_waves[w]]),
		                                        std::vector<double>(&norms.second[countWaves], &norms.second[countWaves + n_waves[w]]),
		fs_masses, efficiency, 12, 23, "_"+wave_names[w]);
		histWaveOrtho.Write();
		countWaves += n_waves[w];
	}

	TH2D histAll = makeDalizFromModel(std::vector<std::complex<double> >(&prodAmps[0], &prodAmps[n_model]),
	                                  std::vector<std::shared_ptr<amplitude> >(&model[0], &model[n_model]),
	                                  std::vector<double>(&norms.second[0], &norms.second[n_model]),
	fs_masses, efficiency);
	histAll.Write();
	TH2D histAllOrtho = makeDalizFromModel(std::vector<std::complex<double> >(&prodAmps[0], &prodAmps[n_model]),
	                                       std::vector<std::shared_ptr<amplitude> >(&model[0], &model[n_model]),
	                                       std::vector<double>(&norms.second[0], &norms.second[n_model]),
	fs_masses, efficiency, 12, 23);
	histAllOrtho.Write();
	TH2D histAllCP = makeDalizFromModel(std::vector<std::complex<double> >(&prodAmps[n_model], &prodAmps[2*n_model]),
	                                    std::vector<std::shared_ptr<amplitude> >(&model[n_model], &model[2*n_model]),
	                                    std::vector<double>(&norms.second[n_model], &norms.second[2*n_model]),
	fs_masses, efficiency, 12, 13, "_CP");
	histAllCP.Write();
	TH2D histAllOrthoCP = makeDalizFromModel(std::vector<std::complex<double> >(&prodAmps[n_model], &prodAmps[2*n_model]),
	                                         std::vector<std::shared_ptr<amplitude> >(&model[n_model], &model[2*n_model]),
	                                         std::vector<double>(&norms.second[n_model], &norms.second[2*n_model]),
	fs_masses, efficiency, 12, 23, "_CP");
	histAllOrthoCP.Write();
	TH2D histBg = makeDalizFromModel({prodAmps[2*n_model]}, {model[2*n_model]}, {norms.second[2*n_model]}, fs_masses, unitEfficiency);
//	TH2D histBg = makeDalizFromModel({prodAmps[2*n_model]}, {model[2*n_model]}, fs_masses, unitEfficiency);
	histBg.Write();
	TH2D histBgOrtho = makeDalizFromModel({prodAmps[2*n_model]}, {model[2*n_model]}, {norms.second[2*n_model]}, fs_masses, unitEfficiency, 12, 23);
//	TH2D histBgOrtho = makeDalizFromModel({prodAmps[2*n_model]}, {model[2*n_model]}, fs_masses, unitEfficiency, 12, 23);
	histBgOrtho.Write();
	std::shared_ptr<amplitude> constantAmpl = std::make_shared<constantAmplitude>(std::make_shared<kinematicSignature>(2));
	TH2D histAcceptance      = makeDalizFromModel({1.}, {constantAmpl}, {1.}, fs_masses, efficiency,12,13,"TimesEfficiency");
	histAcceptance.Write();
	TH2D histAcceptanceOrtho = makeDalizFromModel({1.}, {constantAmpl}, {1.}, fs_masses, efficiency,12,23,"TimesEfficiency");
	histAcceptanceOrtho.Write();

	outFile->Close();
	std::cout << "DostDrocess::main(...): INFO: Finished creating Dalitz plots... exit" << std::endl;
	return 0;
}
