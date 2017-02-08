#ifndef LOGLIKELIHOOD__
#define LOGLIKELIHOOD__
#include<vector>
#include<memory>
#include<complex>
#include"amplitude.h"
#include"integrator.h"
#include"kinematicSignature.h"
class logLikelihood {
	public:
		logLikelihood (std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral);

		std::pair<double, std::vector<std::complex<double> > > fitNlopt   (std::vector<std::complex<double> >& parameters);
		std::pair<double, std::vector<std::complex<double> > > fitROOT    (std::vector<std::complex<double> >& parameters);
		double                                                 nloptCall  (const std::vector<double> &x, std::vector<double> &grad) const;
		double                                                 operator() (const double* prodAmps)                                  const;
		double                                                 eval       (std::vector<std::complex<double> >& prodAmps)            const;
		std::vector<double>                                    Deval      (std::vector<std::complex<double> >& prodAmps)            const;
		std::vector<std::vector<double> >                      DDeval     (std::vector<std::complex<double> >& prodAmps)            const;

		bool                                                   loadDataPoints (const std::vector<std::vector<double> >& dataPoints);

		bool                    setFixFirstPhase (bool flag)                {_fixFirstPhase = flag; return true;}
		bool                    fixFirstPhase    ()                   const {return _fixFirstPhase;}
		size_t                  getNpar          ()                   const;
		size_t                  nAmpl            ()                   const {return _nAmpl;}

	private:
		bool                                             _fixFirstPhase;
		std::shared_ptr<kinematicSignature>              _kinSignature;
		mutable size_t                                   _nCalls;
		size_t                                           _nCallsPrint;
		size_t                                           _nAmpl;
		size_t                                           _nPoints;
		std::shared_ptr<integrator>                      _integral;
		std::vector<std::shared_ptr<amplitude> >         _amplitudes;
		std::vector<std::vector<std::complex<double> > > _points;
};
#endif//LOGLIKELIHOOD__
