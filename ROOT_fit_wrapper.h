#ifndef ROOT_FIT_WRAPPER__
#define ROOT_FIT_WRAPPER__
#include<memory>
#include<complex>
#include<string>
#include"logLikelihood.h"

class ROOT_fit_wrapper {
	public:
		ROOT_fit_wrapper(std::shared_ptr<logLikelihoodBase> ll);

		static double DoEval(const double* xx);
//		static std::vector<double> DoGradient(const double* xx);

		std::pair<double, std::vector<std::complex<double> > > fit(const std::vector<std::complex<double> >& startVals);

		static std::shared_ptr<logLikelihoodBase> __ll;
	protected:
		size_t                             _minimizerStrategy;
		double                             _initStepSize;
		std::shared_ptr<logLikelihoodBase> _ll;
};
#endif//ROOT_FIT_WRAPPER__
