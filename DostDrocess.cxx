#include "Dnine.h"
#include <sstream>

int main(int argc, char* argv[]) {

	std::cout << "DostDrocess::main(...): INFO: Branch file endings are: '" << branchFileEnding << "' '" << branchIntegralFileEnding << "'" << std::endl;

	if (argc < 2) {
		std::cout << "DostDrocess::main(...): ERROR: No input file given" << std::endl;
		return 1;
	}

// TODO:
	bool do_dalitz_plot = false;
	bool do_hessian     = false;
	bool do_covariance  = false;
//
	const int softpionSign = 0;

	std::string resultFileName = argv[1];
	for (int i = 2; i < argc; ++i) {
		std::string argString(argv[i]);
		if (argString == "-hessian") {
			std::cout << "DostDrocess::main(...): INFO: Create the hessian" << std::endl;
			do_hessian = true;
		} else if (argString == "-dalitz") {
			std::cout << "DostDrocess::main(...): INFO: Create Dalitz plots" << std::endl;
			do_dalitz_plot = true;
		} else if (argString == "-coma") {
			std::cout << "DostDrocess::main(...): INFO: Create covariance matrix" << std::endl;
			do_covariance = true;
		} else {
			std::cout << "DostDrocess::main(...): WARNING: Unknown command line argument '" << argv[i] << std::endl;
		}
	}

	if (do_hessian && do_covariance) {
		std::cout << "DostDrocess::main(...): INFO: Cannot do covariance matrix and hessian" << std::endl;
		return 1;
	}


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

	const double numLim = 1.e-5;

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
			if (line == "./Dnine.exe " or line == "Dnine.exe ") {
				continue;
			}
			if (line == "-prior ") {
				std::cout << "DostDrocess::main(...): INFO: Using prior" << std::endl;
				prior = true;
			} else if(line == "-copy ") {
				std::cout << "DostDrocess::main(...): INFO: Copy KpiRight to KpiWrong" << std::endl;
				copy = true;
			} else {
				int i = atoi(line.c_str());
				if (i == 0 and line != "0 ") {
					std::cout << "DostDrocess::main(...): ERROR: Invalid option: '" << line << "'" << std::endl;
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
		if (integral->loadIntegrals("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ps_"+integral_file_name+"_regular." + branchFileEnding,"/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ac_"+integral_file_name+"_regular." + branchFileEnding)) {
			std::cout << "DostDrocess::main(...): INFO: Integral loaded" << std::endl;
		} else {
			std::cout << "DostDrocess::main(...): ERROR: Could not load integral" << std::endl;
			return 1;
		}

		std::string integral_cp_file_name =  "integral_cp_model_" + freeString;
		if (integral_cp->loadIntegrals("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ps_"+integral_cp_file_name+"_regular." + branchFileEnding,"/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ac_"+integral_cp_file_name+"_regular." + branchFileEnding)) {
			std::cout << "DostDrocess::main(...): INFO: CP integral loaded" << std::endl;
		} else {
			std::cout << "DostDrocess::main(...): ERROR: Could not load CP integral" << std::endl;
			return 1;
		}

		std::string integral_bg_file_name = "integral_bg[" + bg_amplitude->name() + "]";
		if (integral_bg->loadIntegrals("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ps_"+integral_bg_file_name+"_regular." + branchFileEnding,"/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ac_"+integral_bg_file_name+"_regular." + branchFileEnding)) {
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
		dataFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/BELLE_data.root";
	} else {
		dataFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/BELLE_bothSidebandsHigherMD.root";
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
	if (pow(inputLL-evalLL, 2) > numLim) {
		std::cout << "DostDrocess::main(...): ERROR: Difference to input LL too big (> " << pow(numLim,.5) << ")" << std::endl;
		return 1;
	}

	if (do_hessian) {
		std::vector<std::string> splitted = utils::splitString(resultFileName, '/');
		std::vector<std::string> parts    = utils::splitString(splitted[splitted.size()-1], '_');

		std::string hessianFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/BELLE_fit_results/hessians/BELLE_parameterHessian";
		for (size_t i = 2; i < parts.size(); ++i) {
			hessianFileName += "_" + parts[i];
		}
		std::cout << "DostDrocess::main(...): INFO: Hessian file name is '" << hessianFileName << "'" << std::endl;
		std::ofstream outFile;		

		std::vector<double> realParameters = getParamsFromProdAmps(prodAmps, n_waves, true, copy, fixToZeroMap);
		std::vector<std::vector<double> > hessian = ll->makeFinalHessian(realParameters, ll->DDeval(prodAmps));

		outFile.open(hessianFileName.c_str());
		outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
		for (size_t i = 0; i < hessian.size(); ++i) {
			for (size_t j = 0; j < hessian.size(); ++j) {
				outFile << hessian[i][j] << " ";
			}
			outFile << std::endl;
		}
		outFile.close();
		std::cout << "DostDrocess::main(...): INFO: Hessian '" + hessianFileName + "' file created." << std::endl;
		std::cout << "now run '/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/sampleFromhessian.py " << hessianFileName << "'" << std::endl;

	}
	if (do_covariance) {
		std::ofstream outFile;		
		std::vector<std::string> splitted = utils::splitString(resultFileName, '/');
		std::vector<std::string> parts    = utils::splitString(splitted[splitted.size()-1], '_');

		std::string sampleFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/BELLE_fit_results/hessians/BELLE_sample";
		for (size_t i = 2; i < parts.size(); ++i) {
			sampleFileName += "_" + parts[i];
		}
		std::cout << "DostDrocess::main(...): INFO: Sample file name is '" << sampleFileName << "'" << std::endl;

		std::vector<double> sampledVector = utils::readRealValuesFromTextFile(sampleFileName);

		std::vector<double> realParameters = getParamsFromProdAmps(prodAmps, n_waves, true, copy, fixToZeroMap);
		size_t nSample = sampledVector.size()/realParameters.size();
		std::cout << "DostDrocess::main(...): INFO: Sample size is " << nSample << std::endl;

		std::vector<std::vector<double> > samplingDeltas = utils::reshape(sampledVector, realParameters.size(), nSample);
		std::vector<std::vector<double> > coma(2*prodAmps.size(), std::vector<double>(2*prodAmps.size(), 0.));

		for (size_t s = 0; s < nSample; ++s) {
			for (size_t p = 0; p < realParameters.size(); ++p ) {
				realParameters[p] += samplingDeltas[s][p];
			}
			std::vector<std::complex<double> > DprodAmps = ll->fullParamsToProdAmps(ll->getFullParameters(realParameters));
			for (size_t a1 = 0; a1 < DprodAmps.size(); ++a1) {
				for(size_t a2 = 0; a2 < DprodAmps.size(); ++a2) {
					std::complex<double> d1 = DprodAmps[a1] - prodAmps[a1];
					std::complex<double> d2 = DprodAmps[a2] - prodAmps[a2];

					coma[2*a1  ][2*a2  ] += d1.real() * d2.real() / nSample;
					coma[2*a1+1][2*a2  ] += d1.imag() * d2.real() / nSample;
					coma[2*a1  ][2*a2+1] += d1.real() * d2.imag() / nSample;
					coma[2*a1+1][2*a2+1] += d1.imag() * d2.imag() / nSample;
				}
			}
			for (size_t p = 0; p < realParameters.size(); ++p ) {
				realParameters[p] -= samplingDeltas[s][p];
			}
		}

		std::string comaFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/BELLE_fit_results/hessians/BELLE_COMA";
		for (size_t i = 2; i < parts.size(); ++i) {
			comaFileName += "_" + parts[i];
		}
		outFile.open(comaFileName.c_str());
		outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
		for (size_t p = 0; p < coma.size(); ++p ){
			for (size_t q = 0; q < coma.size(); ++q) {
				outFile << coma[p][q] << " ";
			}
			outFile << std::endl;
		}
		outFile.close();
		std::cout << "DostDrocess::main(...): INFO: COMA file '" << comaFileName << "' witten" << std::endl;
	}

	if(do_dalitz_plot) {

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
	}
	return 0;
}
