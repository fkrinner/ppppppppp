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
#include"modelAmplitude.h"
#include"modelGenerator.h"
#include"utils.h"
#include"logLikelihood.h"

int main() {
	size_t seed = size_t( time(NULL) );
	std::cout << "Seed: " << seed << std::endl;
//	size_t seed = 13041988 + 13434 + 1;
	srand(seed);
//	utils::opening();

	const double mPi =  .13957018;
//	const double mD0 = 1.86484;
	const double mDp = 1.8696200;

	const size_t nTries         = 10;
	const size_t nPoints        = 1000000;
	const size_t integralPoints = 1000000;

	const bool bose             = true;
	const bool extended         = true;

	const bool writeGenerated   = false;
	const bool writeAllAttempts = false;

	const bool useF0   = true;
	const bool useRho  = true;
	const bool useF2   = true;

	const bool freeF0  = true;
	const bool freeRho = true;
	const bool freeF2  = false;

	std::ofstream outFile;

	std::vector<double> fsMasses = {mPi, mPi, mPi};
	std::shared_ptr<threeParticleMassGenerator> gen = std::make_shared<threeParticleMassGenerator>(mDp, fsMasses);

	std::shared_ptr<constant> cnst = std::make_shared<constant>();

	std::shared_ptr<simpleBW> f0  = std::make_shared<simpleBW>( .98, .1 );
//	std::shared_ptr<simpleBW> f0  = std::make_shared<simpleBW>(1.4 , .1 );
	std::shared_ptr<simpleBW> rho = std::make_shared<simpleBW>( .77, .16);
	std::shared_ptr<simpleBW> f2  = std::make_shared<simpleBW>(1.27, .2 );

	std::shared_ptr<zeroMode0pp> zero0pp = std::make_shared<zeroMode0pp>(mDp*mDp, mPi*mPi);
	std::shared_ptr<zeroMode1mm> zero1mm = std::make_shared<zeroMode1mm>(mDp*mDp, mPi*mPi);

	std::shared_ptr<sameMassZeroS> Sangle = std::make_shared<sameMassZeroS>(mPi);
	std::shared_ptr<sameMassOneP>  Pangle = std::make_shared<sameMassOneP>(mPi);
	std::shared_ptr<sameMassTwoD>  Dangle = std::make_shared<sameMassTwoD>(mPi);

	std::shared_ptr<sameMassZeroSnonRelativistic> SangleNR = std::make_shared<sameMassZeroSnonRelativistic>(mPi);
	std::shared_ptr<sameMassOnePnonRelativistic>  PangleNR = std::make_shared<sameMassOnePnonRelativistic>(mPi);
	std::shared_ptr<sameMassTwoDnonRelativistic>  DangleNR = std::make_shared<sameMassTwoDnonRelativistic>(mPi);

	std::shared_ptr<ratioOfDependences> ratioS = std::make_shared<ratioOfDependences>(Sangle, SangleNR);
	std::shared_ptr<ratioOfDependences> ratioP = std::make_shared<ratioOfDependences>(Pangle, PangleNR);
	std::shared_ptr<ratioOfDependences> ratioD = std::make_shared<ratioOfDependences>(Dangle, DangleNR);

	std::shared_ptr<threeParticleIsobaricAmplitude> Swave = std::make_shared<threeParticleIsobaricAmplitude>(bose, "0mp0ppPiS", f0,  SangleNR);
	std::shared_ptr<threeParticleIsobaricAmplitude> Pwave = std::make_shared<threeParticleIsobaricAmplitude>(bose, "0mp1mmPiP", rho, PangleNR);
	std::shared_ptr<threeParticleIsobaricAmplitude> Dwave = std::make_shared<threeParticleIsobaricAmplitude>(bose, "0mp2ppPiP", f2,  DangleNR);

	std::shared_ptr<threeParticleIsobaricAmplitude> Szero = std::make_shared<threeParticleIsobaricAmplitude>(bose, "Szero", zero0pp, Sangle);
	std::shared_ptr<threeParticleIsobaricAmplitude> Pzero = std::make_shared<threeParticleIsobaricAmplitude>(bose, "Pzero", zero1mm, Pangle);


	std::vector<std::shared_ptr<amplitude> > amplitudes;
	if (useF2) {
		amplitudes.push_back(Dwave);
	}
	if (useRho) {
		amplitudes.push_back(Pwave);
	}
	if (useF0) {
		amplitudes.push_back(Swave);
	}

	const bool freedIsobar = freeF0 or freeRho or freeF2;

	size_t nAmpl = amplitudes.size();
	std::shared_ptr<integrator> integral = std::make_shared<integrator>(integralPoints, gen, amplitudes);
	std::cout << "Starting integration" << std::endl;
	integral->integrate();
	std::cout << "Finished integration" << std::endl;

	std::vector<double> normalizations;
	std::vector<std::complex<double> > transitionAmplitudes;
	for (size_t a = 0; a < nAmpl; ++a) {
		std::pair<bool, std::complex<double> > diag = integral->element(a,a);
		if (not diag.first) {
			std::cerr << "Could not get diagonal element" << std::endl;
			return 1;
		}
		normalizations.push_back(1./pow(diag.second.real(), .5));
		transitionAmplitudes.push_back(std::complex<double>(utils::random2(), utils::random2()));
	}
	std::vector<std::complex<double> > taToSet = {std::complex<double>( 1.,0.),std::complex<double>( 0.,1.), std::complex<double>(-1.,0.)};

	for (size_t i = 0; i < std::min(taToSet.size(), transitionAmplitudes.size()); ++i) {
		transitionAmplitudes[i] = taToSet[i];
	}

	std::shared_ptr<modelAmplitude> model = std::make_shared<modelAmplitude>(transitionAmplitudes, amplitudes, normalizations);

	modelGenerator generator(model, gen);
	std::cout << "Starting generation" << std::endl;
	std::vector<std::vector<double> > generatedPoints = generator.generateDataPoints(nPoints, nPoints); // Won't loose anything, if the burn-in is as long as the sample
	std::cout << "Finished generation" << std::endl;

	if (writeGenerated) {
		const double s23base   = mDp*mDp + 3*mPi*mPi;
		outFile.open("./generated.dat");
		for (std::vector<double> vv : generatedPoints) {
			const double s12 = vv[1];
			const double s13 = vv[2];
			const double s23 = s23base - s12 - s13;
			outFile << s12 << " " << s13 << " " << s23 << std::endl;
		}		
		outFile.close();
	}

	const double width = .02;
	
	std::vector<double> binning = {.278, .3, .3 + width};
	while (binning[binning.size()-1] < 1.760) {
		binning.push_back(binning[binning.size()-1] + width);
	}

	std::vector<double> binningF0 = { .278,  .320,  .360,  .400,  .440,  .480,  .520,  .560,  .600,
                             .640,  .680,  .720,  .760,  .800,  .840,  .880,  .920,  .930,
                             .940,  .950,  .960,  .970,  .980,  .990, 1.000, 1.010, 1.020,
                            1.030, 1.040, 1.050, 1.060, 1.070, 1.080, 1.120, 1.160,
                            1.200, 1.240, 1.280, 1.320, 1.360, 1.400, 1.440, 1.480,
                            1.520, 1.560, 1.600, 1.640, 1.680, 1.720, 1.760};

//	binningF0 = binning;

	for (size_t i = 0; i< binningF0.size(); ++i) {
		binningF0[i] = binningF0[i]*binningF0[i];
	}

	std::vector<double> binningRho = {0.278,0.32, 0.36, 0.4,  0.44, 0.48, 0.52, 0.56, 0.6,  0.64, 0.68, 0.7,  0.72,
                             0.74, 0.76, 0.78, 0.8,  0.82, 0.84, 0.86, 0.88, 0.9,  0.92, 0.96, 1.0,  1.04,
                             1.08, 1.12, 1.16, 1.2,  1.24, 1.28, 1.32, 1.36, 1.4,  1.44, 1.48, 1.52, 1.56,
                             1.6,  1.64, 1.68, 1.72, 1.76};

//	binningRho = binning;

	for (size_t i = 0; i< binningRho.size(); ++i) {
		binningRho[i] = binningRho[i]*binningRho[i];
	}

	std::vector<double> binningF2 = {0.278,0.32, 0.36, 0.4,  0.44, 0.48, 0.52, 0.56, 0.6,  0.64, 0.68, 
					0.72, 0.76, 0.8,  0.84, 0.88, 0.92, 0.96, 1.0,  1.04, 1.08, 1.12, 
					1.16, 1.18, 1.2,  1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 
					1.38, 1.4,  1.44, 1.48, 1.52, 1.56, 1.6,  1.64, 1.68, 1.72, 1.76 };

//	binningF2 = binning;
	for (size_t i = 0; i < binningF2.size(); ++i) {
		binningF2[i] = binningF2[i]*binningF2[i];
	}

	if (freedIsobar) {
		if (useF0 and freeF0) {
			outFile.open("./binningF0.dat");
			for (const double& b : binningF0) {
				outFile << b << " ";
			}
			outFile.close();
		}
		if (useRho and freeRho) {
			outFile.open("./binningRho.dat");
			for (const double& b : binningRho) {
				outFile << b << " ";
			}
			outFile.close();
		}
		if (useF2 and freeF2) {
			outFile.open("./binningF2.dat");
			for (const double& b : binningF2) {
				outFile << b << " ";
			}
			outFile.close();
		}
		outFile.open("./nDp.dat");
		if (useF0) {
			if (freeF0) {
				outFile << binningF0.size()-1;
			} else {
				outFile << 1;
			}
		} else {
			outFile << 0;
		}
		if (useRho) {
			if (freeRho) {
				outFile << " " << binningRho.size()-1;
			} else {
				outFile << " " << 1;
			}
		} else {
			outFile << " " << 0;
		}
		if (useF2) {
			if (freeF2) {
				outFile << " " << binningF2.size()-1;
			} else {
				outFile << " 1";
			}
		} else {
			outFile << " 0";
		}
		outFile.close();
	}

	std::vector<std::shared_ptr<amplitude> > amplitudesFit = amplitudes;
	std::shared_ptr<integrator> integralFit                = integral;

	if (freedIsobar) {
//		amplitudesFit = {Swave, Pwave};
		amplitudesFit = {};
		if (useF2) {
			if (freeF2) {
				for (size_t b = 0; b < binningF2.size() - 1; ++b) {
					std::shared_ptr<stepLike> step = std::make_shared<stepLike>(binningF2[b], binningF2[b+1]);
					std::shared_ptr<threeParticleIsobaricAmplitude> stepWave = std::make_shared<threeParticleIsobaricAmplitude>(bose, std::string("0mp2pp[") + std::to_string(b) + std::string("]PiD"), step, DangleNR);
					amplitudesFit.push_back(stepWave);
				}
			} else {
				amplitudesFit.push_back(Dwave);
			}
		}
		if (useF0) {
			if (freeF0) {
				for (size_t b = 0; b < binningF0.size() - 1; ++b) {
					std::shared_ptr<stepLike> step  = std::make_shared<stepLike>(binningF0[b],binningF0[b+1]);
					std::shared_ptr<threeParticleIsobaricAmplitude> stepWave = std::make_shared<threeParticleIsobaricAmplitude>(bose, std::string("0mp0pp[") + std::to_string(b) + std::string("]PiS"), step, SangleNR);	
					amplitudesFit.push_back(stepWave);
				}
			} else {
				amplitudesFit.push_back(Swave);
			}
		}
		if (useRho) {
			if (freeRho) {
				for (size_t b = 0; b < binningRho.size() - 1; ++b) {
					std::shared_ptr<stepLike> step  = std::make_shared<stepLike>(binningRho[b],binningRho[b+1]);
					std::shared_ptr<threeParticleIsobaricAmplitude> stepWave = std::make_shared<threeParticleIsobaricAmplitude>(bose, std::string("0mp1mm[") + std::to_string(b) + std::string("]PiP"), step, PangleNR);	
					amplitudesFit.push_back(stepWave);
				}
			} else { 
				amplitudesFit.push_back(Pwave);
			}
		}
		integralFit = std::make_shared<integrator>(integralPoints, gen, amplitudesFit);
		std::cout << "Starting integration 2" << std::endl;
		integralFit->integrate();
		std::cout << "Finished integration 2" << std::endl;
	}



	logLikelihood ll(amplitudesFit, integralFit);
	integralFit->writeToFile("./integral.dat");
//	std::cout << "artificial stopping point after writing of the integrals" << std::endl;
//	return 0;
	ll.setFixFirstPhase(true);
	std::cout << "Start loading data points" << std::endl;
	ll.loadDataPoints(generatedPoints);
	std::cout << "Finished loading data points" << std::endl;


	ll.setExtended(extended);

	double avgAmpl = pow(nPoints, .5);
	double bestLike = 1./0.;
	std::vector<std::complex<double> > bestVals;
	std::cout << "Start fitting" << std::endl;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	std::vector<std::complex<double> >pa(ll.nAmpl());
	if (writeAllAttempts) {
		outFile.open("./allAttempts.dat");
	}
	for (size_t _ = 0; _ < nTries; ++_) {
		std::cout << "At: #" << _ << std::endl;
		for (size_t a = 0; a < ll.nAmpl(); ++a) {
			pa[a] = avgAmpl * std::complex<double>(utils::random2(), utils::random2());
		}	
		try {
			std::pair<double, std::vector<std::complex<double> > > retVal = ll.fitNlopt(pa);
			if (retVal.first < bestLike) {
				std::cout << "try " << _ << " bestLike improved from " << bestLike << " to " << retVal.first << std::endl;
				bestLike = retVal.first;
				bestVals = retVal.second;
			}
			if (writeAllAttempts) {
				outFile << retVal.first;
				for (std::complex<double> p : retVal.second) {
					outFile << " " << p.real() << " " << p.imag();
				}
				outFile << _ << " " << std::endl;
			}
		} catch(std::runtime_error err) {
			std::cout << "Caught an error... continue" << std::endl;
			continue;
		}
	}
	if (writeAllAttempts) {
		outFile.close();
	}

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout << "Finished fitting " << std::ctime(&end_time) << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
	std::vector<std::complex<double> > doppel(1,transitionAmplitudes[0] * bestVals[1] / transitionAmplitudes[1] / bestVals[0]);
	if (not freedIsobar) {
		doppel.push_back(transitionAmplitudes[1] * bestVals[2] / transitionAmplitudes[2] / bestVals[1]);
	}
	std::cout << "Doppel:";
	for (std::complex<double> d : doppel){
		std::cout << " " << d;
	}
	std::cout << std::endl;

	double intens = integralFit->totalIntensity(bestVals);
	std::cout << intens << " events described, " << nPoints << " points generated" << std::endl;
	std::cout << "The best likelihood is " << bestLike << std::endl;
	std::vector<std::vector<double> > hessian = ll.DDeval(bestVals);

	std::string outFileName = std::string("./hessian_") + std::to_string(seed) + std::string(".dat");
	outFile.open(outFileName.c_str());
	for (size_t i = 0; i < 2*ll.nAmpl(); ++i) {
		for (size_t j = 0; j < 2*ll.nAmpl(); ++j) {
			outFile << hessian[i][j] << " ";
		}
		outFile << std::endl;
	}
	outFile.close();

	if (freedIsobar) {
		outFileName = std::string("./amplitudes_") + std::to_string(seed) + std::string(".dat");
		outFile.open(outFileName.c_str());
		for (size_t a = 0; a < ll.nAmpl(); ++a) {
			std::complex<double> amp = bestVals[a];
			outFile << amp << std::endl;
		}
		outFile.close();
		outFileName = std::string("./amplitudes_") + std::to_string(seed) + std::string(".gnu");
		outFile.open(outFileName.c_str());
		for (size_t a = 0; a < ll.nAmpl(); ++a) {
			std::complex<double> amp = bestVals[a];
			outFile << a << " " << norm(amp) << " " << amp.real() << " " << amp.imag() << std::endl;
		}
		outFile.close();
	}
	outFileName = "./bestLikeFile";
	outFile.open(outFileName.c_str(), std::ios_base::app);
	outFile << seed << " " << bestLike <<std::endl;
	outFile.close();
	return 0;
}
/* Testcode for hessians...
	double delta = 0.00001;
	std::complex<double> deltaR(delta,0.);
	std::complex<double> deltaI(0.,delta);
	double v = ll.eval(predictedFreedIsobar);
	std::vector<double> val = ll.Deval(predictedFreedIsobar);
	std::vector<std::vector<double> > hess = ll.DDeval(predictedFreedIsobar);
	for(size_t i = 0; i < ll.nAmpl(); ++i) {
		predictedFreedIsobar[i] += deltaR;
		double dr = (ll.eval(predictedFreedIsobar) - v)/delta;
		predictedFreedIsobar[i] += deltaI - deltaR;
		double di = (ll.eval(predictedFreedIsobar) - v)/delta;
		predictedFreedIsobar[i] -= deltaI;
		std::cout << "DR: " << dr << " " << val[2*i] << std::endl;
		std::cout << "DI: " << di << " " << val[2*i+1] << std::endl;
	}
	for(size_t i = 0; i < ll.nAmpl(); ++i) {
		predictedFreedIsobar[i] += deltaR;
		std::vector<double> DR = ll.Deval(predictedFreedIsobar);
		predictedFreedIsobar[i] += deltaI - deltaR;
		std::vector<double> DI = ll.Deval(predictedFreedIsobar);
		predictedFreedIsobar[i] -= deltaI;
		for (size_t j = 0; j < ll.nAmpl(); ++j) {
			double D = (DR[2*j] - val[2*j])/delta;
			std::cout << "DRDR " << D << " " << hess[2*i][2*j] << std::endl;
			D = (DR[2*j+1] - val[2*j+1])/delta;
			std::cout << "DIDR " << D << " " << hess[2*i][2*j+1] << std::endl;
			D = (DI[2*j] - val[2*j])/delta;
			std::cout << "DRDI " << D << " " << hess[2*i+1][2*j] << std::endl;
			D = (DI[2*j+1] - val[2*j+1])/delta;
			std::cout << "DIDI " << D << " " << hess[2*i+1][2*j+1] << std::endl;
			std::cout << "-------------------------------------" << std::endl;
		}
		std::cout << "=====================================" << std::endl;
	}
	return 0;
*/

/* Code for writing amplitudes
	const double sMin  = 0.;
	const double sMax  = mDp*mDp;
	const size_t nBins = 1000;
	const double width = (sMax - sMin)/(nBins);
	outFile.open("zeroModeAmpls.dat");
	for (size_t i = 0; i < nBins; ++i) {
		double sI = sMin + (.5+i)*width;
		for (size_t j = 0; j < nBins; ++j) {
			double sJ = sMin + (.5*j)*width;
			std::vector<double> kin = {mDp*mDp, sI, sJ};
			if (not gen->isValidPoint(kin)) {
				continue;
			}
			double zm0pp = Szero->eval(kin).real();
			double zm1mm = Pzero->eval(kin).real();
			outFile << sI << " " << sJ << " " << zm0pp <<" "<< zm1mm << std::endl;
		}
	}
	outFile.close();
	return 0;
*/

