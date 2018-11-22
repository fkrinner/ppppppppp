#include"Dnine.h"

int main(int argc, char* argv[]) {

	TH2D ful = TH2D("full_dalitz", "fullDalitz", 100, 0.,2.2, 100 ,0.3 , 3.3);
	const size_t seed = size_t( time(NULL) );
	std::cout << "Dnine::main(...): INFO: Seed: " << seed << std::endl;
	srand(seed);

	bool copyRightToWrong = false;
	bool usePrior = false;

	double kin1max = 2.8;

	const int softpionSign = 0;

	const double f_sig  = 6.75192e-01;// +- 9.11090e-04
	const double f_rand = 8.95779e-01;// +- 2.14394e-03
	const double f_CP   = 0.492;

	const double total_signal_coefficient = f_sig + (1.-f_sig)*f_rand*(1.-f_CP);
	const double total_CP_coefficient    = (1.-f_sig)*(1.-f_rand)*f_CP;
	const double total_bg_coefficient    = (1.-f_sig)*f_rand;

	std::cout << "Signal coefficient: " << total_signal_coefficient << "; CP coefficient: " << total_CP_coefficient << "; BG coefficient: " << total_bg_coefficient << std::endl;

	std::vector<bool> free_map          = {false, false, false, false, false, false, false, false , false};
	std::vector<bool> fixToZeroList     = {false, false, false, false, false, false, false, false , false};
	for (int i = 1; i < argc; ++i) {
		std::string argvString = std::string(argv[i]);
		if (argvString == "-prior") {
			usePrior = true;
			std::cout << "Dnine::main(...): ERROR: Use inteference prior" << std::endl;
			continue;
		}
		if (argvString == "-copy") {
			copyRightToWrong = true;
			std::cout << "Dnine::main(...): ERROR: Copy parameters from KpiRight to KpiWrong" << std::endl;
			continue;
		}
		int waveIndex = atoi(argv[i]);
		if (waveIndex >= (int)free_map.size() || waveIndex < -9) {
			std::cout << "Dnine::main(...): ERROR: Invalid wave index to free or fix to zero: " << waveIndex << std::endl;
			return 1;
		}
		if (waveIndex < 0) {
			int indexToFix = - waveIndex - 1;
			std::cout << "Dnine::main(...): INFO: Fixing wave '" << wave_names[indexToFix] << "' (index " << indexToFix << ") to zero" << std::endl;
			fixToZeroList[indexToFix] = true;
		} else {
			std::cout << "Dnine::main(...): INFO: Freeing wave '" << wave_names[waveIndex] << "' (index " << waveIndex << ")" << std::endl;
			free_map[waveIndex] = true;
		}
	}
	std::string freeString = "";
	for (bool free : free_map) {
		freeString += std::to_string(free);
	}

	const std::vector<size_t> n_waves = getNwaves(free_map);
// // full fixed n_waves = {1,2,1,1,3,1,1,3,1,};


	const size_t integral_points = 6000*10*10*10*10; // Merged two 30000000 integrals for the model // Comment for ease bugfix
	const std::vector<double> fs_masses = {mPi, mKs, mPi};
	std::shared_ptr<threeParticleMassGenerator> generator = std::make_shared<threeParticleMassGenerator>(mD0, fs_masses, std::make_shared<kinematicSignature>(2));
	std::shared_ptr<efficiencyFunction> efficiency        = std::make_shared<BELLE_DtoKpipi_efficiency>();
	efficiency->setKin1max(kin1max);

	std::vector<std::shared_ptr<amplitude> > model       = get_model(free_map, mD0, mPi, mKs);
	std::vector<std::shared_ptr<amplitude> > fixed_model = get_model({false, false, false, false, false, false, false, false, false}, mD0, mPi, mKs);
	const size_t n_model = model.size();

	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integral_points, generator, model, efficiency);
	std::string integral_file_name = "integral_model_" + freeString; 

	std::shared_ptr<integrator> fixed_integral = std::make_shared<integrator>(integral_points, generator, fixed_model, efficiency);
	if (!fixed_integral->loadIntegrals("./integralFiles/ps_integral_model_000000000_regular.cut","./integralFiles/ac_integral_model_000000000_regular.cut")) {
		std::cout << "Dnine::main(...): ERROR: Could not load fixed model integrals" << std::endl;
		return 1;
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: './integralFiles/ps_integral_model_000000000_regular.cut' and './integralFiles/ac_integral_model_000000000_regular.cut'" <<std::endl;
	}

	if (!integral->loadIntegrals("./integralFiles/ps_"+integral_file_name+"_regular.cut","./integralFiles/ac_"+integral_file_name+"_regular.cut")) {
		if (!integral->integrate()) {
			std::cout << "Dnine::main(...): ERROR: Model integration failed" << std::endl;
			return 1;
		};
		integral->writeToFile("./integralFiles/ps_"+integral_file_name+".cut", false);
                integral->writeToFile("./integralFiles/ac_"+integral_file_name+".cut", true);
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: " << "'./integralFiles/ps_"+integral_file_name+"_regular.cut' and './integralFiles/ac_"+integral_file_name+"_regular.cut'" <<std::endl;
	}

	std::vector<std::shared_ptr<amplitude> > model_cp = get_model(free_map, mD0, mPi, mKs, true);
	std::shared_ptr<integrator> integral_cp = std::make_shared<integrator>(integral_points, generator, model_cp, efficiency);
	std::string integral_cp_file_name =  "integral_cp_model_" + freeString; 

	if (!integral_cp->loadIntegrals("./integralFiles/ps_"+integral_cp_file_name+"_regular.cut","./integralFiles/ac_"+integral_cp_file_name+"_regular.cut")) {
		if (!integral_cp->integrate()) {
			std::cout << "Dnine::main(...): ERROR: CP model integration failed" << std::endl;
			return 1;

		}
		integral_cp->writeToFile("./integralFiles/ps_"+integral_cp_file_name+".cut", false);
		integral_cp->writeToFile("./integralFiles/ac_"+integral_cp_file_name+".cut", true);
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: " << "'./integralFiles/ps_"+integral_cp_file_name+"_regular.cut' and './integralFiles/ac_"+integral_cp_file_name+"_regular.cut'" <<std::endl;
	}

	std::shared_ptr<amplitude> bg_amplitude = get_bg_amplitude();
	std::shared_ptr<integrator> integral_bg = std::make_shared<integrator>(integral_points, generator, std::vector<std::shared_ptr<amplitude> >(1,bg_amplitude), efficiency);

	std::string integral_bg_file_name = "integral_bg[" + bg_amplitude->name() + "]";

	if (!integral_bg->loadIntegrals("./integralFiles/ps_"+integral_bg_file_name+"_regular.cut","./integralFiles/ac_"+integral_bg_file_name+"_regular.cut")) {
		std::cout << "Dnine::main(...): WARNING: Could not load bg integral. Integrate" << std::endl;
		if (!integral_bg->integrate()) {
			std::cout << "Dnine::main(...): ERROR: Background integration failed" << std::endl;
			return 1;
		}
		integral_bg->writeToFile("./integralFiles/ps_"+integral_bg_file_name+".cut",false);
		integral_bg->writeToFile("./integralFiles/ac_"+integral_bg_file_name+".cut",true);
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: " << "'./integralFiles/ps_"+integral_bg_file_name+"_regular.cut' and './integralFiles/ac_"+integral_bg_file_name+"_regular.cut'" <<std::endl;
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

	std::shared_ptr<logLikelihood_withPrior> ll = std::make_shared<logLikelihood_withPrior>(model, integral);

	std::cout << "Dnine::main(...): INFO: Loading data points to log-likelihood..." << std::endl;
	if (!ll->loadDataPoints(dataPoints, 25)) {
		std::cout << "Dnine::main(...): ERROR: Could not load data points" << std::endl;
		return 1;
	}
	std::cout << "Dnine::main(...): INFO: Finished loading data points to log-likelihood." << std::endl;
	if (!ll->setExtended(true)) {
		std::cout << "Dnine::main(...): ERROR: Could not set extended" << std::endl;
		return 1;
	}

	double CP_scale_factor = pow(total_CP_coefficient/total_signal_coefficient, .5);


	if (!doTheFixingAndCopying(ll, n_waves, CP_scale_factor, copyRightToWrong, fixToZeroList, true)) {
		std::cout << "Dnine::main(...): ERROR: Could not do the fixing and copying. Abort..." << std::endl;
		return 1;
	}
// -- P A R A M E T E R   F I X I N G   A N D   C O P Y I N G   E N D S! ! ! 

//#define DO_ZERO_MODE_PRIORS
#ifdef DO_ZERO_MODE_PRIORS
	const double zeroModePriorStrength = 20000.;
	size_t countZero = 0;
	while (true) {
		std::string zeroFileName = "./zeroModeFiles/"+freeString+"_"+std::to_string(countZero)+".cut";
		std::vector<std::complex<double> > zeroMode = utils::readComplexValuesFromTextFile(zeroFileName, true);
		if (zeroMode.size() == ll->nAmpl()) {
			if (!ll->addPriorDirection(zeroModePriorStrength,zeroMode)) {
				std::cout << "Dnine::main(...): ERROR: Could not add zero mode from '" << zeroFileName << "'" << std::endl;
				return 1;
			} else {
				std::cout << "Dnine::main(...): INFO: Added prior for a zero mode with an intensity of: " << integral->totalIntensity(zeroMode, true) << "/" << integral->totalIntensity(zeroMode, false) << " (prior strength " << zeroModePriorStrength << ")" << std::endl;
			}
		} else if (zeroMode.size() == 0) {
			break;
		} else {
			std::cout << "Dnine::main(...): ERROR: Zero mode with unmatching size found: " << zeroMode.size() << std::endl;
			return 1;
		}
		++countZero;
	}
#endif//DO_ZERO_MODE_PRIORS

//#define DO_ALL_WAVE_PRIOR
#ifdef DO_ALL_WAVE_PRIOR
	const double allWavePriorStrength = .0001;
	std::cout << "Dnine::main(...): INFO: Add a prior with strength " << allWavePriorStrength << " for all waves" << std::endl;
	for (size_t a = 0; a < n_model; ++a) {
		std::vector<std::complex<double> > singleAmplShape(ll->nAmpl(), std::complex<double>(0.,0.));
		singleAmplShape[a] = std::complex<double>(1.,0.);
		if (!ll->addPriorDirection(allWavePriorStrength, singleAmplShape)) {
			std::cout << "Dnine::main(...): ERROR: Could not add single wave prior strength for wave " << a << std::endl;
			return 1;
		}
	}
#endif// DO_ALL_WAVE_PRIOR
	if (usePrior) {
		ll-> setInterferencePriorStrength(1000000.);
	}
//	ll->setNcallsPrint(20);
	std::cout << "Dnine::main(...): INFO: Start fitting" << std::endl;

// HERE ARE THE NFITS // // // // // // // // // // // // // // // // // // // /
	const size_t nFits = 1;
// // // // // // // // // // // // // // // // // // // // // // // // // // //

	double bestLL = std::numeric_limits<double>::infinity();
	std::vector<double> allLL(nFits,0.);
	std::vector<std::complex<double> > bestStart;
	std::vector<std::complex<double> > bestResult;
	
	for (size_t f = 0; f < nFits; ++f) {

//		std::vector<std::complex<double> > startValues;	
//		std::vector<std::complex<double> > BELLE_result_values = getStartValuesForBelleParameters(free_map, integral->getNormalizations(false).second, fixed_integral->getNormalizations(false).second);
//		randomizeStartValues(BELLE_result_values, .1, 1.9); // Add some randomization to maybe circumvent nlopt failure
//		for (size_t a = 0; a < BELLE_result_values.size(); ++a) {
//			BELLE_result_values[a] = std::complex<double>(utils::random2(), utils::random2());
//		}
//		for (size_t a = 0; a < n_model; ++a) {
//			BELLE_result_values.push_back(CP_scale_factor * BELLE_result_values[a]);
//		}
//		BELLE_result_values.push_back(std::complex<double>(pow(nData * total_bg_coefficient/integral_bg->totalIntensity(std::vector<std::complex<double> >(1,std::complex<double>(1.,0.)), true), .5),0.));

		std::complex<double> bg_amplitude(pow(nData * total_bg_coefficient/integral_bg->totalIntensity(std::vector<std::complex<double> >(1,std::complex<double>(1.,0.)), true), .5), 0.);
		std::vector<std::complex<double> > startValues = getRandomizedStartValues( n_waves, CP_scale_factor, true, copyRightToWrong, bg_amplitude);

		{ // Normalize to number of events... keep the BG coefficient
			std::complex<double> BG_coeff = startValues[2*n_model];
			std::vector<std::complex<double> > BG_only_vec(startValues.size(), std::complex<double>(0.,0.));
			BG_only_vec[2*n_model] = BG_coeff;
			startValues[2*n_model] = std::complex<double>(0.,0.);
			double scaleFactor = pow(integral->totalIntensity(startValues, true)/(nData - integral->totalIntensity(BG_only_vec, true)), -.5);
			for (size_t a = 0; a < 2*n_model; ++a) {
				startValues[a] *= scaleFactor;
			}
			startValues[2*n_model] = BG_coeff;
		}
		std::vector<double> realParameters = getParamsFromProdAmps(startValues, n_waves, true, copyRightToWrong, fixToZeroList);

		std::pair<double, std::vector<std::complex<double> > > retVal;

		std::cout << "Dnine::main(...): INFO: Start values gotten... Fit!" << std::endl;

		std::cout << "Dnine::main(...): INFO: Number of parameters is " << ll->getNpar() << " (should match " << realParameters.size() << ")" << std::endl;
		ll->resetNcalls();
		try {
			retVal = ll->fitNlopt(realParameters);
			std::cout << "Dnine::main(...): INFO: Success after " << ll->nCalls() << " function calls" << std::endl;
		} catch(...) {
			continue;
		}

		allLL[f] = retVal.first;
		if (retVal.first < bestLL) {
			bestLL     = retVal.first;
			bestResult = retVal.second;
			bestStart  = startValues;
		}
		if (std::isinf(bestLL)) {
			bestStart = startValues;
		}
	} 
	if (std::isinf(bestLL)) {
		std::cout << "Dnine::main(...): ERROR: Not a single fit coverged. Call the fitter again to reproduce the error, then abort" << std::endl;
		std::vector<double> realParameters = getParamsFromProdAmps(bestStart, n_waves, true, copyRightToWrong, fixToZeroList);
		ll->fitNlopt(realParameters);
		std::cout << "Dnine::main(...): ERROR: Abort!" << std::endl;
		return 1;
	}

	std::string outFileName = "./BELLE_startValues/BELLE_startValues" + freeString;
	outFileName += "_"+std::to_string(seed)+".cut";
	std::ofstream outFile;
	outFile.open(outFileName.c_str());
	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	const double startValueIntentity = integral->totalIntensity(bestStart, true);
	for (std::complex<double>& amp : bestStart) {
		outFile << amp << std::endl;
	}
	outFile.close();

	outFileName = "./BELLE_fit_results/BELLE_fit_" + freeString;
	outFileName += "_"+std::to_string(bestLL)+"_"+std::to_string(seed)+".cut";
	outFile.open(outFileName.c_str());
	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (std::complex<double>& amp : bestResult) {
		outFile << amp << std::endl;
	}
	outFile.close();

	outFileName = "./BELLE_fit_results/worseLikelihoods/corresponding_BELLE_fit_" + freeString;
	outFileName += "_"+std::to_string(seed)+".cut";
	outFile.open(outFileName.c_str());
	outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (double ll : allLL) {
		outFile << ll << std::endl;
	}
	outFile.close();

	outFileName = "./BELLE_fit_results/fitInfos/BELLE_fit_" + freeString;
	outFileName += "_"+std::to_string(seed)+".cut";
	outFile.open(outFileName.c_str());
	for (int i = 0; i < argc; ++i) {
		outFile << argv[i] << " " << std::endl;
	}
	outFile.close();


//#define WRITE_HESSIANS
#ifdef WRITE_HESSIANS
	if (writeHessians) {

		outFileName = "./BELLE_fit_results/hessians/BELLE_hessian_" + freeString;
		outFileName += "_"+std::to_string(bestLL)+"_"+std::to_string(seed)+".cut";
		std::vector<std::vector<double> > hessian = ll->makeFinalHessian(ll->prodAmpsToFullParams(bestResult), ll->logLikelihood::DDeval(bestResult));
		outFile.open(outFileName.c_str());
		std::vector<std::pair<size_t,double> > fixedParams = ll->getFixedParameters();
		const size_t dim = hessian.size() + fixedParams.size();
		size_t countI = 0;
		for (size_t i = 0; i < dim; ++i) {
			bool iIsFixed = false;
			for (std::pair<size_t,double> fix: fixedParams) {
				if (i == fix.first) {
					iIsFixed = true;
					break;
				}
			}
			size_t countJ = 0;
			for (size_t j = 0; j < dim; ++j) {
				bool jIsFixed = false;
				for (std::pair<size_t,double> fix: fixedParams) {
					if (j == fix.first) {
						jIsFixed = true;
						break;
					}
				}
				if (iIsFixed || jIsFixed) {
					outFile << "0. ";
	//				if (jIsFixed) {
	//					++countJ;
	//				}
				} else {
					outFile << hessian[countI][countJ] << " ";
					++countJ;
				}
			}
			if (!iIsFixed) {
				++countI;
			}
			outFile << std::endl;
		}
		outFile.close();
	}
#endif//WRITE_HESSIANS

/// Getting interference of the result
	double norm = 0.;
	for (size_t a = 0; a < n_model; ++a) {
		norm += std::norm(bestResult[a]);
	}
	norm = pow(norm,.5);
	std::vector<std::complex<double> > modelSector(bestResult.size(), std::complex<double>(0.,0.));
	for (size_t a = 0; a < n_model; ++a) {
		modelSector[a] = bestResult[a]/norm;
	}
	std::cout << "-----------------------------------------------------------------" << std::endl;
	std::cout << "The number of described events is: " << integral->totalIntensity(bestResult, true) << " (should be " << nData << "; start values had " << startValueIntentity << ")" << std::endl;
	std::cout << "The intensity after interference is: " << integral->totalIntensity(modelSector, false) << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;

	bool makeAmplitudeDalitzPlots = utils::updateBestLLfile("bestLL_"+freeString, bestLL, std::to_string(seed));
	if (makeAmplitudeDalitzPlots) {
		std::shared_ptr<efficiencyFunction> unitEfficiency = std::make_shared<threeParticlPerfectEfficiency>(std::make_shared<kinematicSignature>(2));

		std::string outFileNameDalitz = "dalitzPlotResults_"+freeString	+".root";
		std::pair<bool, std::vector<double> > norms_ac = integral->getNormalizations(true);
		std::pair<bool, std::vector<double> > norms    = integral->getNormalizations(false);
		if (!norms.first or !norms_ac.first) {
			std::cout << "Dnine::main(...): ERROR: Could not get normalizations" << std::endl;
			return 1;
		}
		for (size_t a = 0; a < norms.second.size(); ++a) {
			norms.second[a] = 1./norms.second[a];;//orms_ac.second[a];///norms.second[a];
		}
		TFile* outFile = new TFile(outFileNameDalitz.c_str(), "RECREATE");
		TH2D dataHist  = makeDataDalitz(dataPoints, fs_masses);
		dataHist.Write();
		TH2D dataHistOrtho  = makeDataDalitz(dataPoints, fs_masses, 12,23);
		dataHistOrtho.Write();
		const bool singleWavePlots = false;
		if (singleWavePlots) {
			for (size_t a = 0; a < 2*n_model; ++a) {
				TH2D hist = makeDalizFromModel({bestResult[a]}, {model[a]}, {norms.second[a]}, fs_masses, efficiency);
	//			TH2D hist = makeDalizFromModel({bestResult[a]}, {model[a]}, fs_masses, efficiency);
				hist.Write();
				TH2D histOrtho = makeDalizFromModel({bestResult[a]}, {model[a]}, {norms.second[a]}, fs_masses, efficiency, 12,23);
	//			TH2D histOrtho = makeDalizFromModel({bestResult[a]}, {model[a]}, fs_masses, efficiency, 12,23);
				histOrtho.Write();
			}
		}
		TH2D histAll = makeDalizFromModel(std::vector<std::complex<double> >(&bestResult[0], &bestResult[n_model]), 
		                                  std::vector<std::shared_ptr<amplitude> >(&model[0], &model[n_model]), 
		                                  std::vector<double>(&norms.second[0], &norms.second[n_model]), 
		fs_masses, efficiency);
		histAll.Write();
		TH2D histAllOrtho = makeDalizFromModel(std::vector<std::complex<double> >(&bestResult[0], &bestResult[n_model]), 
		                                       std::vector<std::shared_ptr<amplitude> >(&model[0], &model[n_model]), 
		                                       std::vector<double>(&norms.second[0], &norms.second[n_model]), 
		fs_masses, efficiency, 12, 23);
		histAllOrtho.Write();
		TH2D histAllCP = makeDalizFromModel(std::vector<std::complex<double> >(&bestResult[n_model], &bestResult[2*n_model]), 
		                                    std::vector<std::shared_ptr<amplitude> >(&model[n_model], &model[2*n_model]), 
		                                    std::vector<double>(&norms.second[n_model], &norms.second[2*n_model]),
		fs_masses, efficiency, 12, 13, "_CP");
		histAllCP.Write();
		TH2D histAllOrthoCP = makeDalizFromModel(std::vector<std::complex<double> >(&bestResult[n_model], &bestResult[2*n_model]), 
		                                         std::vector<std::shared_ptr<amplitude> >(&model[n_model], &model[2*n_model]), 
		                                         std::vector<double>(&norms.second[n_model], &norms.second[2*n_model]),
		fs_masses, efficiency, 12, 23, "_CP");
		histAllOrthoCP.Write();

		TH2D histBg = makeDalizFromModel({bestResult[2*n_model]}, {model[2*n_model]}, {norms.second[2*n_model]}, fs_masses, unitEfficiency);
//		TH2D histBg = makeDalizFromModel({bestResult[2*n_model]}, {model[2*n_model]}, fs_masses, unitEfficiency);
		histBg.Write();
		TH2D histBgOrtho = makeDalizFromModel({bestResult[2*n_model]}, {model[2*n_model]}, {norms.second[2*n_model]}, fs_masses, unitEfficiency, 12, 23);
//		TH2D histBgOrtho = makeDalizFromModel({bestResult[2*n_model]}, {model[2*n_model]}, fs_masses, unitEfficiency, 12, 23);
		histBgOrtho.Write();
		outFile->Close();
		std::cout << "Dnine::main(...): INFO: Finished creating Dalitz plots... exit" << std::endl;
		return 0;
	}
	return 0;
}
