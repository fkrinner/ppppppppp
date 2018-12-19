#include "Dnine.h"
#include <sstream>

int main() {

	std::cout << "D-GenerateMC::main(...): INFO: Branch file endings are: '" << branchFileEnding << "' '" << branchIntegralFileEnding << "'" << std::endl;

	const size_t      integral_points = 6000*10*10*10*10; // Merged two 30000000 integrals for the model // Comment for ease bugfix
	std::vector<bool> freeMap         = {false, false, false, false, false, false, false, false, false};
	std::vector<bool> fixToZeroMap    = {false, false, false, false, false, false, false, false, false};

	const double f_sig                    = 6.75192e-01;// +- 9.11090e-04
	const double f_rand                   = 8.95779e-01;// +- 2.14394e-03
	const double f_CP                     = 0.492;
	const double total_signal_coefficient = f_sig + (1.-f_sig)*f_rand*(1.-f_CP);
	const double total_CP_coefficient     = (1.-f_sig)*(1.-f_rand)*f_CP;
	const double total_bg_coefficient     = (1.-f_sig)*f_rand;
	const double CP_scale_factor          = pow(total_CP_coefficient/total_signal_coefficient, .5);

	std::string freeString = "000000000";

	const std::vector<double>                   fs_masses  = {mPi, mKs, mPi};
	std::shared_ptr<threeParticleMassGenerator> generator  = std::make_shared<threeParticleMassGenerator>(mD0, fs_masses, std::make_shared<kinematicSignature>(2));
	std::shared_ptr<efficiencyFunction>         efficiency = std::make_shared<BELLE_DtoKpipi_efficiency>();

	const std::vector<size_t> n_waves = getNwaves(freeMap);

	std::vector<std::shared_ptr<amplitude> > model    = get_model(freeMap, mD0, mPi, mKs);

	const size_t n_model = model.size();

	std::vector<std::shared_ptr<amplitude> > model_cp = get_model(freeMap, mD0, mPi, mKs, true);
	std::shared_ptr<amplitude> bg_amplitude           = get_bg_amplitude(fs_masses);

	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integral_points, generator, model, efficiency);
	std::shared_ptr<integrator> integral_cp = std::make_shared<integrator>(integral_points, generator, model_cp, efficiency);
	std::shared_ptr<integrator> integral_bg = std::make_shared<integrator>(integral_points, generator, std::vector<std::shared_ptr<amplitude> >(1,bg_amplitude), efficiency);

	{
		std::string integral_file_name = "integral_model_" + freeString;
		if (integral->loadIntegrals("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ps_"+integral_file_name+"_regular." + branchFileEnding,"/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ac_"+integral_file_name+"_regular." + branchFileEnding)) {
			std::cout << "D-GenerateMC::main(...): INFO: Integral loaded" << std::endl;
		} else {
			std::cout << "D-GenerateMC::main(...): ERROR: Could not load integral" << std::endl;
			return 1;
		}

		std::string integral_cp_file_name =  "integral_cp_model_" + freeString;
		if (integral_cp->loadIntegrals("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ps_"+integral_cp_file_name+"_regular." + branchFileEnding,"/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ac_"+integral_cp_file_name+"_regular." + branchFileEnding)) {
			std::cout << "D-GenerateMC::main(...): INFO: CP integral loaded" << std::endl;
		} else {
			std::cout << "D-GenerateMC::main(...): ERROR: Could not load CP integral" << std::endl;
			return 1;
		}

		std::string integral_bg_file_name = "integral_bg[" + bg_amplitude->name() + "]";
		if (integral_bg->loadIntegrals("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ps_"+integral_bg_file_name+"_regular." + branchFileEnding,"/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/integralFiles/ac_"+integral_bg_file_name+"_regular." + branchFileEnding)) {
			std::cout << "D-GenerateMC::main(...): INFO: BG integral loaded" << std::endl;
		} else {
			std::cout << "D-GenerateMC::main(...): ERROR: Could not load BG integral" << std::endl;
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
	std::vector<std::complex<double> > prodAmps = getBelleProdAmps();
	const size_t nA = prodAmps.size();
	for (size_t a = 0; a < nA; ++a) {
		prodAmps.push_back(prodAmps[a] * CP_scale_factor);
	}
	prodAmps.push_back(std::complex<double>(0.,0.));
	double intens = integral->totalIntensity(prodAmps, true);
	prodAmps[prodAmps.size() - 1] = std::complex<double>(pow(intens*total_bg_coefficient/(total_CP_coefficient+total_signal_coefficient),.5),0.);

	std::vector<double> norms = integral->getNormalizations().second;

	std::vector<std::shared_ptr<amplitude> > generationAmplitudes;
	generationAmplitudes.push_back(std::make_shared<modelAmplitude>(std::vector<std::complex<double> >       (&prodAmps[0],   &prodAmps[n_model]), 
	                                                                std::vector<std::shared_ptr<amplitude> > (&model[0],      &model[n_model]),
	                                                                std::vector<double>                      (&norms[0],      &norms[n_model]), "model"));
	generationAmplitudes.push_back(std::make_shared<modelAmplitude>(std::vector<std::complex<double> >       (&prodAmps[n_model],   &prodAmps[2*n_model]), 
	                                                                std::vector<std::shared_ptr<amplitude> > (&model[n_model],      &model[2*n_model]),
	                                                                std::vector<double>                      (&norms[n_model],      &norms[2*n_model]), "CP_model"));
	generationAmplitudes.push_back(std::make_shared<modelAmplitude>(std::vector<std::complex<double> >       (&prodAmps[2*n_model],   &prodAmps[2*n_model+1]), 
	                                                                std::vector<std::shared_ptr<amplitude> > (&model[2*n_model],      &model[2*n_model+1]),
	                                                                std::vector<double>                      (&norms[2*n_model],      &norms[2*n_model+1]),"BG"));

	modelGenerator mg(generationAmplitudes, generator, efficiency);

	const size_t nGenerate = 1000000;
	std::vector<std::vector<double> > dataPoints = mg.generateDataPoints(nGenerate,nGenerate);

	std::string outFileName = "MC_data";
	std::ofstream outFile;
	outFile.open(outFileName.c_str());
	for (const std::vector<double>& point : dataPoints) {
		outFile << point[0] << " " << point[1] << " " << point[2] << std::endl;
	}
	outFile.close();
	return 0;
}
