#ifndef AMPLITUDE__
#define AMPLITUDE__
#include<string>
#include<complex>
#include<vector>
#include<memory>
#include"massShape.h"
#include"kinematicSignature.h"
#include"angularDependence.h"

class amplitude {
	public:
		amplitude (std::shared_ptr<kinematicSignature> kinSignature, std::string name);

		virtual std::complex<double> eval (const std::vector<double>& kin) const;
		double  intens (const std::vector<double>& kin) const {return std::norm(eval(kin));}

		std::string                         name         () const {return _name;}
		size_t                              nKin         () const {return _kinSignature->nKin();}
		std::shared_ptr<kinematicSignature> kinSignature () const {return _kinSignature;}


	protected:
		std::shared_ptr<kinematicSignature> _kinSignature;
		std::string                         _name;

};

class threeParticleIsobaricAmplitude : public amplitude {
	public:
		threeParticleIsobaricAmplitude(bool boseSymmetrize, std::string name, std::shared_ptr<massShape> shape, std::shared_ptr<angularDependence> angDep);
		
		std::complex<double> eval(const std::vector<double>& kin) const override;
		std::complex<double> evalSingleSymmTerm(const std::vector<double>& kin) const;

	private:
		bool                               _bose;
		std::shared_ptr<massShape>         _massShape;
		std::shared_ptr<angularDependence> _angDep;
};

class threeParticlaIsobaricAmplitudeNoBose : public amplitude {
	public:
		threeParticlaIsobaricAmplitudeNoBose(size_t isobarIndex, std::string name, std::shared_ptr<massShape> shape, std::shared_ptr<angularDependence> angDep, std::vector<double> fsMasses);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		size_t                             _isobarIndex;
		double                             _sumFSmasses;
		std::shared_ptr<massShape>         _massShape;
		std::shared_ptr<angularDependence> _angDep;
	
};

class dalitzMonomialAmplitude : public amplitude {
	public:
		dalitzMonomialAmplitude(std::shared_ptr<kinematicSignature> kinSignature, double exponent1, double exponent2);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		double _exponent1;
		double _exponent2;

};

class dalitzPolynomialAmplitude : public amplitude {
	public:
		dalitzPolynomialAmplitude(std::shared_ptr<kinematicSignature> kinSignature, const std::string configurationFile, double xMin, double xMax, double yMin, double yMax);

		std::complex<double> eval(const std::vector<double>& kin) const override;
		std::pair<double, double> getXY(const std::vector<double>& kin) const;

	private:
		size_t              _nTerms;
		double              _xMin;
		double              _xMax;
		double              _xWidth;
		double              _yMin;
		double              _yMax;
		double              _yWidth;
		std::vector<size_t> _xExponents;
		std::vector<size_t> _yExponents;
		std::vector<double> _coefficients;
};

class mCosTintensPolynomial : public amplitude {
	public:
		mCosTintensPolynomial(std::shared_ptr<kinematicSignature> kinSignature, const std::string configurationFile, double motherMass, std::vector<double> fsMasses, size_t isobarCombination = 12);

		std::complex<double> eval(const std::vector<double>& kin) const override;
		std::pair<double, double> getMcosT(const std::vector<double>& kin) const;
		std::pair<double, double> getXY(const std::pair<double, double> mCosT) const;

		bool setCosTLimits(const std::pair<double,double> newLimits);
		bool setMlimits(const std::pair<double,double> newLimits);
	private:
		size_t                   _isobarCombination;
		size_t                   _nTerms;
		double                   _mWidth;
		double                   _cosTwidth;
		std::pair<double,double> _mLimits;
		std::pair<double,double> _cosTlimits;
		std::vector<size_t>      _xExponents;
		std::vector<size_t>      _yExponents;
		std::vector<double>      _coefficients;
		std::vector<double>      _fsMasses;
};

class lookupAmplitudeIntens : public amplitude {
	public:
		lookupAmplitudeIntens(std::shared_ptr<kinematicSignature> kinSignature, const std::string& name, double sMinX, double widthX, double sMinY, double widthY, const std::vector<std::vector<double> > & intensities);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	protected:
		size_t _nX;
		size_t _nY;

		double _sMinX;
		double _widthX;

		double _sMinY;
		double _widthY;

		std::vector<std::vector<double> > _data;
};
#endif//AMPLITUDE__
