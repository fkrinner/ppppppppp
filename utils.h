#ifndef UTILS__SLITU
#define UTILS__SLITU
#include<random>
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<complex>
#include<iomanip>
#include<memory>

#include"logLikelihood.h"

namespace utils {
	void opening() {
		std::cout << "This is ppppppppp, a\n Program for\n Pseudoscalar to\n Pseudoscalar,\n Pseudoscalar,\n Pseudoscalar\n Processes for\n Pushing our\n Paper\n Publication" << std::endl;
	}

	double random() {
		double retVal = ((double) rand() / (RAND_MAX));
		return retVal;
	}

	double random2() { 
		return 2*random() -1.;	
	}

	double breakupMomentumSquared(const double s, const double s1, const double s2) {
		return (s*s + s1*s1 + s2*s2 - 2*(s*s1 + s*s2 + s1*s2))/4/s;
	}

	double degToRad(double deg) {
		return deg * M_PI / 180.;
	}

	bool checkComplexDouble() {
		std::vector<std::complex<double> > vals = {std::complex<double>(1.,2.), std::complex<double>(3.,4.)};
		double* val = (double*)&vals[0];
		bool works = true;
		if (val[0] != 1.) {
			works = false;
		} 
		if (val[1] != 2.) {
			works = false;
		}
		if (val[2] != 3.) {
			works = false;
		}
		if (val[3] != 4.) {
			works = false;
		}
		if (!works) {
			std::cout << "utils::checkComplexDouble(): ERROR: Complex double may not be reinterpreted as re, im array" << std::endl;
		}
		return works;
	}

	std::vector<double> readRealValuesFromTextFile(const std::string& inFileName, bool allowZeroLength = false) {
		std::vector<double> retVal;
		double v;
		std::ifstream fin(inFileName.c_str(), std::ifstream::in);
		while (fin>>v) {
			retVal.push_back(v);
		}
		if (retVal.size() == 0 and allowZeroLength) {
			std::cout << "utils::readRealValuesFromTextFile(...): ERROR: No values could be read" << std::endl;
			throw;
		}
		return retVal;
	}

	std::vector<std::vector<double> > reshape(const std::vector<double>& inVector, size_t dimX, size_t dimY) {
		std::vector<std::vector<double> > retVal(dimY, std::vector<double>(dimX, 0.));
		size_t countX = 0;
		size_t countY = 0;
		for (const double& val : inVector) {
			retVal[countY][countX] = val;
			++countX;
			if (countX == dimX) {
				countX = 0;
				++countY;
			}
		}
		if (countX != 0) {
			std::cout << "utils::reshape(...): ERROR: Did not finish at line end" << std::endl;
			throw;
		}
		if (countY != dimY) {
			std::cout << "utils::reshape(...): ERROR: Did not find specified number of lines" << std::endl;
			throw;
		}
		return retVal;	
	}

	std::vector<std::complex<double> > readComplexValuesFromTextFile(const std::string& inFileName, bool allowZeroLength = false) {
		std::vector<std::complex<double> > retVal;
		std::complex<double> c;
		std::ifstream fin(inFileName.c_str(), std::ifstream::in);
		while (fin>>c) {
			retVal.push_back(c);
		}
		if (retVal.size() == 0 and !allowZeroLength) {
			std::cout << "utils::readComplexValuesFromTextFile(...): ERROR: No values could be read" << std::endl;
			throw;
		}
		return retVal;
	}

	bool updateBestLLfile(const std::string& llFileName, double newValue, const std::string& additionalInfo = "") { // Returns "true", <value> is better than the value in <llFileName> or the file doe not yet exist
		double val;
		bool updateFile = false;
		{
			std::ifstream fin(llFileName.c_str(), std::ifstream::in);
			if (fin >> val) {
				if (newValue < val) {
					updateFile = true;
				}
			} else {
				updateFile = true;
			}
		}
		if (updateFile) {
			std::ofstream outFile;
			outFile.open(llFileName.c_str());
			outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1);
			outFile << newValue << std::endl;
			outFile << additionalInfo;
			outFile.close();
		}
		return updateFile;
	}

	std::pair<bool, std::vector<std::vector<std::complex<double> > > > readComplexMatrixFromTextFile(const std::string& inFileName, size_t dim) {
		std::vector<std::vector<std::complex<double> > > retVal(dim, std::vector<std::complex<double> >(dim));
		size_t col = 0;
		size_t lin = 0;
		std::complex<double> c;
		std::ifstream fin(inFileName.c_str(), std::ifstream::in);
		while (fin>>c) {
			retVal[lin][col] = c;
			col += 1;
			if (col == dim) {
				col  = 0;
				lin += 1;
			}
		}
		if (col != 0 || lin != dim) {
			std::cout << "utils::readComplexMatrixFromTextFile(...): ERROR: Mismatch in matrix dimensions: ("<<col<<","<<lin<<")= (col,lin) != (0,dim) = (0," <<dim<<")"<<std::endl;
			return std::pair<bool, std::vector<std::vector<std::complex<double> > > >(false, std::vector<std::vector<std::complex<double> > >());
		}
		return std::pair<bool, std::vector<std::vector<std::complex<double> > > >(true, retVal);
	}

	bool isSanePoint(const std::vector<double>& point, const std::vector<double>& fsMasses) {
		const double mPi  = fsMasses[0];
		const double mKs  = fsMasses[1];

		const double p1pK = (point[1] - mPi*mPi - mKs*mKs)/2;
		const double Epi  = pow(point[2],.5)/2.; // Half the isobar mass in the isobar rest frame
		const double EKs  = (point[0] - point[2] - mKs*mKs)/4/Epi; // 4*Epi = 2*mPiPi = 2*pow(kin[2],.5)
		const double pPi  = pow(Epi*Epi - mPi*mPi, .5);
		const double pKs  = pow(EKs*EKs - mKs*mKs, .5);
		const double cosT = (Epi*EKs - p1pK)/pPi/pKs;
		
		if (cosT < -1. or cosT > 1.) {
			return false;
		}

		const double masssum = fsMasses[0]*fsMasses[0] + fsMasses[1]*fsMasses[1] + fsMasses[2]*fsMasses[2];
		if (point[0] - point[1] - point[2] + masssum < 0.) {
			return false;
		}
		if (point[2] < pow(fsMasses[0] + fsMasses[2], 2)) {
			return false;
		}
		if (point[2] > pow(pow(point[0],.5) - fsMasses[1],2)) {
			return false;
		}

		return true;
	}

	std::vector<std::vector<double> > sanitizeBELLEdataPoints(const std::vector<std::vector<double> >& inData, const std::vector<double>& fsMasses) {
		if (fsMasses.size() != 3) {
			std::cout << "utils::sanitizeBELLEdataPoints(...): ERROR: fsMasses has to have size 3" << std::endl;
			throw;
		}
		const size_t nIn = inData.size();
		if (nIn == 0) {
			return std::vector<std::vector<double> >();
		}
		const size_t dim = inData[0].size();
		size_t count = 0;
		std::vector<std::vector<double> > retVal(nIn, std::vector<double>(dim, 0.));
		for (const std::vector<double> & point : inData) {
			if (isSanePoint(point, fsMasses)) {
				retVal[count][0] = point[0];
				retVal[count][1] = point[1];
				retVal[count][2] = point[2];
				++count;
			}
		}
		retVal.resize(count);
		std::cout << "utils::sanitizeDataPoints(...): INFO: " << count << " events of " << nIn << " events are sane" << std::endl;

		return retVal;
	}


	void checkDerivatives(std::shared_ptr<logLikelihoodBase> ll, const std::vector<std::complex<double> > &ampl, double delta, bool doSecond) {
		std::vector<std::complex<double> > pa = ampl;
		const double evl = ll->eval(pa);
		const std::vector<double> Devl = ll->Deval(pa);
		for (size_t a = 0; a < pa.size(); ++a) {
			double dd = delta * pow(std::norm(pa[a]),.5);
			pa[a] += std::complex<double>(dd,0.);
			std::cout << Devl[2*a  ] << " || " << (ll->eval(pa) - evl)/dd << std::endl;
			pa[a] += std::complex<double>(-dd,dd);
			std::cout << Devl[2*a+1] << " || " << (ll->eval(pa) - evl)/dd << std::endl;
			pa[a] += std::complex<double>(0.,-dd);
		}
		if (doSecond) {
			const std::vector<std::vector<double> > DDevl = ll->DDeval(pa);
			std::cout << "utils::checkDerivatives(...): INFO: Second derivative gotten. Checking numerically now..." << std::endl;
			std::vector<double> D;
			for (size_t a = 0; a < pa.size(); ++a) {
				double dd = delta * pow(std::norm(pa[a]),.5);
				pa[a] += std::complex<double>(dd,0.);
				D = ll->Deval(pa);
				for (size_t b = 0; b < pa.size(); ++b) {
					if (a == b) {
						std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
					}
					std::cout << DDevl[2*a  ][2*b  ] << " ++ " << (D[2*b  ] - Devl[2*b  ])/dd << std::endl;
					std::cout << DDevl[2*a  ][2*b+1] << " ++ " << (D[2*b+1] - Devl[2*b+1])/dd << std::endl;
					if (a == b) {
						std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
					}
				}
				pa[a] += std::complex<double>(-dd,dd);
				D = ll->Deval(pa);
				for (size_t b = 0; b < pa.size(); ++b) {
					if (a == b) {
						std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
					}
					std::cout << DDevl[2*a+1][2*b  ] << " ++ " << (D[2*b  ] - Devl[2*b  ])/dd << std::endl;
					std::cout << DDevl[2*a+1][2*b+1] << " ++ " << (D[2*b+1] - Devl[2*b+1])/dd << std::endl;
					if (a == b) {
						std::cout << " ! - ! - ! - ! - ! - ! - !" << std::endl;
					}
				}
				pa[a] += std::complex<double>(0.,-dd);
			}
		}
	}

	std::vector<std::string> splitString(const std::string& inputString, const char& splitChar) {
		std::vector<std::string> retVal;
		std::string part = "";
		for (const char& c : inputString) {
			if (c == splitChar) {
				if (part != "") {
					retVal.push_back(part);
					part = "";
				}
			} else {
				part += c;
			}
		}
		if (part != "") {
			retVal.push_back(part);
		}
		return retVal;
	}
}
#endif//UTILS__SLITU

