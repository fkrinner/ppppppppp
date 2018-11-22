#include "ROOT_fit_wrapper.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include<string>

std::shared_ptr<logLikelihoodBase> ROOT_fit_wrapper::__ll = std::shared_ptr<logLikelihoodBase>(nullptr);

ROOT_fit_wrapper::ROOT_fit_wrapper(std::shared_ptr<logLikelihoodBase> ll) : _minimizerStrategy(0), _initStepSize(1.), _ll(ll) {} 

double ROOT_fit_wrapper::DoEval(const double* xx) {
	std::vector<double> grad(0);
	return __ll->nloptCall(std::vector<double>(xx, xx + __ll->getNpar()), grad);
}

//std::vector<double> ROOT_fit_wrapper::DoGradient(const double* xx) {
//	std::vector<double> grad(__ll->getNpar(),0.);
//	__ll->nloptCall(std::vector<double>(xx, xx + __ll->getNpar()), grad);
//	return grad;
//}

std::pair<double, std::vector<std::complex<double> > > ROOT_fit_wrapper::fit(const std::vector<std::complex<double> >& startVals) {
	__ll = _ll;
	ROOT::Math::Functor f(&ROOT_fit_wrapper::DoEval, __ll->getNpar());
	std::shared_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
	min->SetStrategy(_minimizerStrategy);
	min->SetFunction(f);
	const std::vector<double> startPars = _ll->cutGradient(_ll->prodAmpsToFullParams(startVals)); // // // // Use this, since no better option now... (will change the start parameters a bit)
	for (size_t p = 0; p < _ll->getNpar(); ++p) {
		std::string parName = "par_" + std::to_string(p);
		double stepSize = _initStepSize;
		min->SetVariable(p, parName, startPars[p], stepSize);
	}
	min->Minimize();
	const std::vector<double> result(min->X(),min->X() + _ll->getNpar());
	std::vector<std::complex<double> > finalProdAmps = _ll->fullParamsToProdAmps(_ll->getFullParameters(result));
	return std::pair<double, std::vector<std::complex<double> > > (_ll->eval(finalProdAmps) ,finalProdAmps);
}
