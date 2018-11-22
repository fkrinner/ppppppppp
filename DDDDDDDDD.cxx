#include"Dnine.h"

int main() {

	TH2D ful = TH2D("full_dalitz", "fullDalitz", 100, 0.,2.2, 100 ,0.3 , 3.3);
	const size_t seed = size_t( time(NULL) );
	std::cout << "Dnine::main(...): INFO: Seed: " << seed << std::endl;
	srand(seed);

	const double f_sig  = 6.75192e-01;// +- 9.11090e-04
	const double f_rand = 8.95779e-01;// +- 2.14394e-03
	const double f_CP   = 0.492;

	const double total_signal_coefficient = f_sig + (1.-f_sig)*f_rand*(1.-f_CP);
	const double total_CP_coefficient    = (1.-f_sig)*(1.-f_rand)*f_CP;
	const double total_bg_coefficient    = (1.-f_sig)*f_rand;

	const size_t nEvents = 1081486;

	std::cout << "Signal coefficient: " << total_signal_coefficient << "; CP coefficient: " << total_CP_coefficient << "; BG coefficient: " << total_bg_coefficient << std::endl;

	std::vector<std::string> wave_names = {"KpiRightS", "KpiRightP", "KpiRightD", "KpiWrongS", "KpiWrongP", "KpiWrongD", "piPiS", "piPiP", "piPiD"};

	const size_t integral_points = 6000*10*10*10*10; // Merged two 30000000 integrals for the model // Comment for ease bugfix
	const std::vector<double> fs_masses = {mPi, mKs, mPi};

	std::shared_ptr<threeParticleMassGenerator> generator = std::make_shared<threeParticleMassGenerator>(mD0, fs_masses, std::make_shared<kinematicSignature>(2));
	std::shared_ptr<efficiencyFunction> efficiency        = std::make_shared<BELLE_DtoKpipi_efficiency>();

	std::vector<std::shared_ptr<amplitude> > fixed_model = get_model({false, false, false, false, false, false, false, false, false}, mD0, mPi, mKs);

	std::shared_ptr<integrator> fixed_integral = std::make_shared<integrator>(integral_points, generator, fixed_model, efficiency);
	if (!fixed_integral->loadIntegrals("./integralFiles/ps_integral_model_000000000_regular.dat","./integralFiles/ac_integral_model_000000000_regular.dat")) {
		std::cout << "Dnine::main(...): ERROR: Could not load fixed model integrals" << std::endl;
		return 1;
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: './integralFiles/ps_integral_model_000000000_regular.dat' and './integralFiles/ac_integral_model_000000000_regular.dat'" <<std::endl;
	}

	std::vector<std::shared_ptr<amplitude> > model_cp = get_model({false, false, false, false, false, false, false, false, false}, mD0, mPi, mKs, true);
	std::shared_ptr<integrator> integral_cp = std::make_shared<integrator>(integral_points, generator, model_cp, efficiency);
	std::string integral_cp_file_name =  "integral_cp_model_000000000"; 

	if (!integral_cp->loadIntegrals("./integralFiles/ps_"+integral_cp_file_name+"_regular.dat","./integralFiles/ac_"+integral_cp_file_name+"_regular.dat")) {
		if (!integral_cp->integrate()) {
			std::cout << "Dnine::main(...): ERROR: CP model integration failed" << std::endl;
			return 1;

		}
		integral_cp->writeToFile("./integralFiles/ps_"+integral_cp_file_name+".dat", false);
		integral_cp->writeToFile("./integralFiles/ac_"+integral_cp_file_name+".dat", true);
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: " << "'./integralFiles/ps_"+integral_cp_file_name+"_regular.dat' and './integralFiles/ac_"+integral_cp_file_name+"_regular.dat'" <<std::endl;
	}

	std::shared_ptr<amplitude> bg_amplitude = get_bg_amplitude();
	std::shared_ptr<integrator> integral_bg = std::make_shared<integrator>(integral_points, generator, std::vector<std::shared_ptr<amplitude> >(1,bg_amplitude), efficiency);

	std::string integral_bg_file_name = "integral_bg[" + bg_amplitude->name() + "]";

	if (!integral_bg->loadIntegrals("./integralFiles/ps_"+integral_bg_file_name+"_regular.dat","./integralFiles/ac_"+integral_bg_file_name+"_regular.dat")) { // Use two times ps, since the BG parameterization is already acceptance correctezd...
		std::cout << "Dnine::main(...): WARNING: Could not load bg integral. Integrate" << std::endl;
		if (!integral_bg->integrate()) {
			std::cout << "Dnine::main(...): ERROR: Background integration failed" << std::endl;
			return 1;
		}
		integral_bg->writeToFile("./integralFiles/ps_"+integral_bg_file_name+".dat",false);
		integral_bg->writeToFile("./integralFiles/ac_"+integral_bg_file_name+".dat",true);
	} else {
		std::cout << "Dnine::main(...): INFO: Loaded: " << "'./integralFiles/ps_"+integral_bg_file_name+"_regular.dat' and './integralFiles/ac_"+integral_bg_file_name+"_regular.dat'" <<std::endl;
	}

	std::vector<std::complex<double> > BELLE_transition_amplitudes = getBelleProdAmps();
	std::vector<std::complex<double> > scaled_transition_amplitudes(0);

	double CP_scale_factor = pow(total_CP_coefficient/total_signal_coefficient,.5);
	for (std::complex<double> A : BELLE_transition_amplitudes) {
		scaled_transition_amplitudes.push_back(A * CP_scale_factor);
	}
	std::complex<double> BG_transition_amplitude(pow(total_bg_coefficient*(fixed_integral->totalIntensity(BELLE_transition_amplitudes) + integral_cp->totalIntensity(scaled_transition_amplitudes)),.5),0.);

	std::cout << "The signal intensity is: " << fixed_integral->totalIntensity(BELLE_transition_amplitudes) << "; the CP intensity is: " << integral_cp->totalIntensity(scaled_transition_amplitudes) << "; the BG intensity is: " << std::norm(BG_transition_amplitude) << std::endl;

	std::shared_ptr<amplitude> amplitudeOfModel   = std::make_shared<modelAmplitude>(BELLE_transition_amplitudes , fixed_model, fixed_integral->getNormalizations(false).second, "modelAmplitude");
	std::shared_ptr<amplitude> amplitudeOfCPModel = std::make_shared<modelAmplitude>(scaled_transition_amplitudes, model_cp   , integral_cp->getNormalizations(false).second, "CPmodelAmplitude");
	std::shared_ptr<amplitude> backgroundModel    = std::make_shared<modelAmplitude>(std::vector<std::complex<double> > (1,BG_transition_amplitude) , std::vector<std::shared_ptr<amplitude> > (1,bg_amplitude), integral_bg->getNormalizations(false).second, "BGmodelAmplitude");

	modelGenerator generator_of_the_model({amplitudeOfModel, amplitudeOfCPModel, backgroundModel}, generator, efficiency);

	std::vector<std::vector<double> > MC_data_points = generator_of_the_model.generateDataPoints(nEvents, nEvents);
	writeTextDataFile("./Monte-Carlo_data_set.dat", MC_data_points);

	return 0;
}
