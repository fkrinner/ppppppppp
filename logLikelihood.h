#ifndef LOGLIKELIHOOD__
#define LOGLIKELIHOOD__
#include<vector>
#include<memory>
#include<complex>
#include<iostream>
#include"amplitude.h"
#include"integrator.h"
#include"kinematicSignature.h"

class logLikelihoodBase {
	public:
		logLikelihoodBase (size_t nAmpl, std::shared_ptr<kinematicSignature> kinSignature, std::shared_ptr<integrator> integral);

		std::pair<double, std::vector<std::complex<double> > > fitNlopt   (std::vector<std::complex<double> >& parameters);
//		std::pair<double, std::vector<std::complex<double> > > fitROOT    (std::vector<std::complex<double> >& parameters);
		double                                                 nloptCall  (const std::vector<double> &x, std::vector<double> &grad) const;

		virtual double                                         eval       (std::vector<std::complex<double> >& prodAmps)            const 
                                                                       {std::cout << "Base class eval called (" << prodAmps.size() << " parameters), returning 0." << std::endl; return 0.;}
		virtual std::vector<double>                            Deval      (std::vector<std::complex<double> >& prodAmps)            const
                                                                       {std::cout << "Base class Deval called (" << prodAmps.size() << " parameters), returning empty vector" << std::endl; return std::vector<double>();}
		virtual std::vector<std::vector<double> >              DDeval     (std::vector<std::complex<double> >& prodAmps)            const
                                                                       {std::cout << "Base class DDeval called (" << prodAmps.size() << " parameters), returning empty vector" << std::endl; return std::vector<std::vector<double> >();}
		virtual bool                                           loadDataPoints (const std::vector<std::vector<double> >& dataPoints) 
                                                                       {std::cout << "No loadDataPoints(...) method specified. Do not load " << dataPoints.size() << " points" << std::endl; return false;}

		bool                    setFixFirstPhase (bool flag)                {_fixFirstPhase = flag; return true;}
		bool                    fixFirstPhase    ()                   const {return _fixFirstPhase;}
		bool                    setExtended      (bool flag)                {_extended = flag; return true;}
		size_t                  getNpar          ()                   const;
		size_t                  nAmpl            ()                   const {return _nAmpl;}
		size_t                  nCalls           ()                   const {return _nCalls;}
	protected:
		bool                                             _fixFirstPhase;
		bool                                             _extended;
		std::shared_ptr<kinematicSignature>              _kinSignature;
		mutable size_t                                   _nCalls;
		size_t                                           _nCallsPrint;
		size_t                                           _nAmpl;
		size_t                                           _nPoints;
		std::shared_ptr<integrator>                      _integral;
};

class logLikelihood : public logLikelihoodBase {
	public:
		logLikelihood (std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral);

		double                                                 eval       (std::vector<std::complex<double> >& prodAmps)            const override;
		std::vector<double>                                    Deval      (std::vector<std::complex<double> >& prodAmps)            const override;
		std::vector<std::vector<double> >                      DDeval     (std::vector<std::complex<double> >& prodAmps)            const override;

		size_t                                                 getSector(size_t a) const;
		bool                                                   loadDataPoints (const std::vector<std::vector<double> >& dataPoints) override;
		bool                                                   setCoherenceBorders(std::vector<size_t> borders);

		bool                    setNstore        (size_t nStore);
		std::pair<std::vector<double>, std::vector<double> > getStoredPoints() const;
	protected:
		size_t                                           _nSect;
		std::vector<std::shared_ptr<amplitude> >         _amplitudes;
		std::vector<std::vector<std::complex<double> > > _points;
		std::vector<size_t>                              _amplitudeCoherenceBorders;
		std::vector<std::vector<size_t> >                _contributingWaves;

		size_t                                           _nStore;
		mutable size_t                                   _storageLocation;
		mutable std::vector<std::complex<double> >       _lastPar;
		mutable std::vector<double>                      _storeEvals;
		mutable std::vector<double>                      _storeSteps;
};

class logLikelihoodAllFree : public logLikelihoodBase {
	public:
		logLikelihoodAllFree (std::vector<double> binning, std::vector<std::shared_ptr<angularDependence> > freedAmplitudes, std::shared_ptr<integrator> integral);

		double                                                 eval       (std::vector<std::complex<double> >& prodAmps)            const override;
		std::vector<double>                                    Deval      (std::vector<std::complex<double> >& prodAmps)            const override;
		std::vector<std::vector<double> >                      DDeval     (std::vector<std::complex<double> >& prodAmps)            const override;

		bool                                                   loadDataPoints (const std::vector<std::vector<double> >& dataPoints) override;
		std::pair<bool, std::pair<size_t, size_t> >            findBin(const std::vector<double>& point) const;

	protected:
		size_t                                                         _nBins;
		std::vector<double>                                            _binning;
		std::vector<std::vector<size_t> >                              _eventsPerBin;
		std::vector<std::vector<std::vector<std::complex<double> > > > _amplitudesInBin;
};
#endif//LOGLIKELIHOOD__
