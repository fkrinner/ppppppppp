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
	
		std::pair<double, std::vector<std::complex<double> > > fitNlopt   (const std::vector<double>& parameters);
		std::pair<double, std::vector<std::complex<double> > > fitNlopt   (const std::vector<std::complex<double> >& parameters);
		double                                                 nloptCall  (const std::vector<double> &x, std::vector<double> &grad) const;

		virtual double                                         eval       (const std::vector<std::complex<double> >& prodAmps)            const 
                                                                       {std::cout << "Base class eval called (" << prodAmps.size() << " parameters), returning 0." << std::endl; return 0.;}
		virtual std::vector<double>                            Deval      (const std::vector<std::complex<double> >& prodAmps)            const
                                                                       {std::cout << "Base class Deval called (" << prodAmps.size() << " parameters), returning empty vector" << std::endl; return std::vector<double>();}
		virtual std::vector<std::vector<double> >              DDeval     (const std::vector<std::complex<double> >& prodAmps)            const
                                                                       {std::cout << "Base class DDeval called (" << prodAmps.size() << " parameters), returning empty vector" << std::endl; return std::vector<std::vector<double> >();}
		virtual bool                                           loadDataPoints (const std::vector<std::vector<double> >& dataPoints, size_t maxNperEvent = 9) 
                                                                       {std::cout << "No loadDataPoints(...) method specified. Do not load " << dataPoints.size() << " points (maxNperEvent = " << maxNperEvent << ")" << std::endl; return false;}

		std::vector<double>                makeFinalGradient(const std::vector<double>& params, const std::vector<double>& fullGradient) const;
		std::vector<std::vector<double> >  makeFinalHessian(const std::vector<double>& params, const std::vector<std::vector<double> >& fullHessian) const;

		std::vector<double>                getFullParameters(const std::vector<double>& params) const;
		std::vector<double>                cutGradient( const std::vector<double>& params) const;

		std::vector<std::complex<double> > fullParamsToProdAmps(const std::vector<double>& params                 ) const;
		std::vector<double>                prodAmpsToFullParams(const std::vector<std::complex<double> >& prodAmps) const;


		std::vector<std::vector<double> >  DDconstrainedProdAmps(const std::vector<double>& params) const;
		std::vector<std::vector<double> >  DprodAmpsNumerical(const std::vector<double>& params, double delta = 1.e-5) const;

		bool                    setNcallsPrint    (size_t nCallPrint)        {_nCallsPrint = nCallPrint; return true;}
		bool                    setExtended       (bool flag)                {_extended = flag; return true;}
		bool                    fixParameter      (size_t nPar, double value);
		bool                    addCopyParameter  (size_t nPar, size_t copyFrom, int nScale = -1, double scaler = 1.);

		size_t                  getNparTot        ()                   const;
		size_t                  getNpar           ()                   const;
		size_t                  nAmpl             ()                   const {return _nAmpl;}
		size_t                  nCalls            ()                   const {return _nCalls;}
		bool                    resetNcalls       ()                         {_nCalls = 0; return true;}

		const std::vector<std::pair<size_t,double> >              getFixedParameters() const {return _fixedParameters;}
		const std::vector<std::tuple<size_t,size_t,int, double> > getCopyParameters () const {return _copyParameters ;}

// //		bool                                                setNumericalGradientCheck(bool flag) {_ngc = flag; return true;} // DELETE
		std::vector<double>                                 getStoreParams() const {return _storeParams;}                    // DELETE
	protected:
		bool                                                _extended;
		std::shared_ptr<kinematicSignature>                 _kinSignature;
		mutable size_t                                      _nCalls;
		size_t                                              _nCallsPrint;
		size_t                                              _nAmpl;
		size_t                                              _nPoints;
		size_t                                              _nScale; // Number of scale parameters
		double                                              _numDelta;
		std::shared_ptr<integrator>                         _integral;
		std::vector<std::pair<size_t,double> >              _fixedParameters;
		std::vector<std::tuple<size_t,size_t,int, double> > _copyParameters;
// //
// //		bool                                                _ngc;         // DELETE
		mutable std::vector<double>                         _storeParams; // DELETE
		
};

class logLikelihood : public logLikelihoodBase {
	public:
		logLikelihood (std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral);

		virtual double                                         eval       (const std::vector<std::complex<double> >& prodAmps)            const override;
		virtual std::vector<double>                            Deval      (const std::vector<std::complex<double> >& prodAmps)            const override;
		virtual std::vector<std::vector<double> >              DDeval     (const std::vector<std::complex<double> >& prodAmps)            const override;

		virtual bool                                           loadDataPoints (const std::vector<std::vector<double> >& dataPoints, size_t maxNperEvent = 9) override;

		size_t                                                 getSector(size_t a) const;
		bool                                                   setCoherenceBorders(std::vector<size_t> borders);

		bool                    setNstore        (size_t nStore);
		bool                    setMaxII         (double maxII) {_maximumIncoherentIntensity = maxII; return true;}
		std::pair<std::vector<double>, std::vector<double> > getStoredPoints() const;		
	protected:
		size_t                                           _nSect;
		double                                           _maximumIncoherentIntensity;
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

class logLikelihood_withPrior : public logLikelihood {
	public: 
		logLikelihood_withPrior (std::vector<std::shared_ptr<amplitude> > amplitudes, std::shared_ptr<integrator> integral);

		virtual double                                         eval       (const std::vector<std::complex<double> >& prodAmps) const override;
		virtual std::vector<double>                            Deval      (const std::vector<std::complex<double> >& prodAmps) const override;
		virtual std::vector<std::vector<double> >              DDeval     (const std::vector<std::complex<double> >& prodAmps) const override;

		bool addPriorDirection(double strength, const std::vector<std::complex<double> > direction);

		bool setInterferencePriorStrength(double strength);

		double                    interferencePriorFunc(double coherent, double incoherent) const;
		std::pair<double, double> DinterferencePriorFunc(double coherent, double incoherent) const;
		std::vector<double>       DDinterferencePriorFunc(double coherent, double incoherent) const;

	protected:
		double                                           _interferencePriorStrength;

		size_t                                           _nPrior;
		std::vector<double>                              _priorStrengths;
		std::vector<std::vector<std::complex<double> > > _priorDirections;

};

class logLikelihoodAllFree : public logLikelihoodBase {
	public:
		logLikelihoodAllFree (std::vector<double> binning, std::vector<std::shared_ptr<angularDependence> > freedAmplitudes, std::shared_ptr<integrator> integral);

		double                                                 eval       (const std::vector<std::complex<double> >& prodAmps)            const override;
		std::vector<double>                                    Deval      (const std::vector<std::complex<double> >& prodAmps)            const override;
		std::vector<std::vector<double> >                      DDeval     (const std::vector<std::complex<double> >& prodAmps)            const override;

		bool                                                   loadDataPoints (const std::vector<std::vector<double> >& dataPoints, size_t maxNperEvent = 9) override;
		std::pair<bool, std::pair<size_t, size_t> >            findBin(const std::vector<double>& point) const;

	protected:
		size_t                                                         _nBins;
		std::vector<double>                                            _binning;
		std::vector<std::vector<size_t> >                              _eventsPerBin;
		std::vector<std::vector<std::vector<std::complex<double> > > > _amplitudesInBin;
};
#endif//LOGLIKELIHOOD__
