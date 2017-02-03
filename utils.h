#ifndef UTILS__SLITU
#define UTILS__SLITU
#include<random>
#include<iostream>
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
}
#endif//UTILS__SLITU

