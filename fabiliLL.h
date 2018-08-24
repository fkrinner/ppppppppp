#ifndef FABILILL__
#define FABILILL__
#include"logLikelihood.h"
#include"fabiliFunction.h"

class fabiliLL : public logLikelihood, public fabiliFunction {
	public:
		fabiliLL(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral);

		double scalarEval(const std::vector<double>& parameters) const override;
		evalType eval(const std::vector<double>& parameters) const override;
};
#endif//FABILILL__
