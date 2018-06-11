#ifndef INTEGRATOR__
#define INTEGRATOR__
#include"amplitude.h"
#include"generator.h"
#include"kinematicSignature.h"
#include"efficiencyFunction.h"
#include<string>
#include<vector>
#include<memory>
#include<complex>

class integrator {
	public:
		integrator(size_t integralPoints, std::shared_ptr<generator> pointGenerator, const std::vector<std::shared_ptr<amplitude> >& amplitudes, std::shared_ptr<efficiencyFunction>& efficiency);

		bool                                             integrate         ();
		bool                                             loadIntegrals(const std::string& psFileName, const std::string& accFileName);
		std::vector<std::vector<std::complex<double> > > getIntegralMatrix (bool accCorr = false)                     const;
		std::pair<bool, std::complex<double> >           element           (size_t i, size_t j, bool accCorr = false) const;

		double                            totalIntensity   (const std::vector<std::complex<double> >& prodAmpl, bool accCorr = false) const;
		std::vector<double>               DtotalIntensity  (const std::vector<std::complex<double> >& prodAmpl, bool accCorr = false) const;
		std::vector<std::vector<double> > DDtotalIntensity (const std::vector<std::complex<double> >& prodAmpl, bool accCorr = false) const;

		bool                                isIntegrated ()                     const {return _isIntegrated;}
		bool                                setNpoints   (size_t n);
		bool                                writeToFile  (std::string fileName, bool accCorr = false) const;
		size_t                              nAmpl        ()                     const {return _nAmpl;}
		std::pair<bool, size_t>             getNpoints   ()                     const; 
		std::pair<bool, std::string>        getWaveName  (size_t i)             const;
		std::shared_ptr<kinematicSignature> kinSignature ()                     const {return _kinSignature;}
	protected:
		bool                                             _isIntegrated;
		std::shared_ptr<kinematicSignature>              _kinSignature;
		size_t                                           _nAmpl;
		size_t                                           _nPoints;
		std::vector<std::shared_ptr<amplitude> >         _amplitudes;
		std::shared_ptr<generator>                       _generator;
		std::shared_ptr<efficiencyFunction>              _efficiency;
		std::vector<std::vector<std::complex<double> > > _integralMatrix;
		std::vector<std::vector<std::complex<double> > > _accCorrIntegralMatrix;
};
#endif//INTEGRATOR__

