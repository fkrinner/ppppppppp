#ifndef FABILILL_OPM__
#define FABILILL_OMP__
#include"logLikelihood.h"
#include"fabiliFunction.h"

class fabiliLL_openMP : public logLikelihood, public fabiliFunction {
	public:
		fabiliLL_openMP(std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral);

		double scalarEval(const std::vector<double>& parameters) const override;
		evalType eval(const std::vector<double>& parameters) const override;
};
#endif//FABILILL_OMP__
