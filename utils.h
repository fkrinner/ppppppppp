#ifndef UTILS__SLITU
#define UTILS__SLITU
#include<random>
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<complex>
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

	std::vector<std::vector<double> > sanitizeDataPoints(const std::vector<std::vector<double> >& inData, const std::vector<double>& fsMasses) {
		if (fsMasses.size() != 3) {
			std::cout << "utils::sanitizeDataPoints(...): ERROR: fsMasses has to have size 3" << std::endl;
			throw;
		}
		const size_t nIn = inData.size();
		if (nIn == 0) {
			return std::vector<std::vector<double> >();
		}
		const size_t dim = inData[0].size();
		const double masssum = fsMasses[0]*fsMasses[0] + fsMasses[1]*fsMasses[1] + fsMasses[2]*fsMasses[2];
		size_t count = 0;
		std::vector<std::vector<double> > retVal(nIn, std::vector<double>(dim, 0.));
		for (const std::vector<double> & point : inData) {
			if (point[0] - point[1] - point[2] + masssum > 0.) {
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
}
#endif//UTILS__SLITU

