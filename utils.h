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

	std::vector<std::vector<std::complex<double> > > readComplexMatrixFromTextFile(const std::string& inFileName, size_t dim) {
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
			std::cout << "Mismatch in matrix dimensions: ("<<col<<","<<lin<<")= (col,lin) != (0,dim) = (0," <<dim<<")"<<std::endl;
			throw;
		}
		return retVal;
	}
}
#endif//UTILS__SLITU

