#ifndef MASSSHAPE__
#define MASSSHAPE__
#include<string>
#include<complex>
#include<vector>
class massShape {
	public:
		massShape (std::string name, std::vector<double> pars, std::vector<std::string> parNames = std::vector<std::string>());

		virtual std::complex<double> eval (double s) const;

		bool                         setParameter (size_t n, double val);
		std::pair<bool, double>      getParameter (size_t n)                   const;
		bool                         setParName   (size_t n, std::string name);
		std::pair<bool, std::string> getParName   (size_t n)                   const;

		std::string                  name         () const {return _name;}

	protected:
		size_t                   _nPar;
		std::string              _name;
		std::vector<double>      _parameters;
		std::vector<std::string> _parameterNames;
};

class simpleBW : public massShape {
	public:
		simpleBW (double mass, double width);

		std::complex<double> eval (double s) const override;
};

class stepLike : public massShape {
	public:
		stepLike (double sMin, double sMax);

		std::complex<double> eval (double s) const override;
};

class constant : public massShape {
	public:
		constant ();

		std::complex<double> eval (double s) const override;
};

class zeroMode0pp : public massShape {
	public:
		zeroMode0pp (double s, double m2);

		std::complex<double> eval (double s) const override;
};

class zeroMode1mm : public massShape {
	public:
		zeroMode1mm (double s, double m2);

		std::complex<double> eval (double s) const override;
};

class polynomialMassShape : public massShape {
	public:
		polynomialMassShape(std::vector<std::complex<double> > coefficients, double baseExponent = 1.);

		std::complex<double> eval (double s) const override;
	protected:
		size_t _polDeg;
		double _baseExponent;	
};

class BELLEbreitWigner : public massShape {
	public:
		BELLEbreitWigner(std::string name, double mass, double width, size_t spin, double motherMass, double bachelorMass, double daughterMass1, double daughterMass2);

		std::complex<double> eval (double s) const override;
	protected:
		size_t _spin;

		double _motherMass;
		double _bachelorMass;
		double _daughterMass1;
		double _daughterMass2;

		double _Rr;		
		double _RD;
};

class BELLE_LASS_KpiS : public massShape {
	public:
		BELLE_LASS_KpiS(const std::vector<double>& parameters, double mPi, double mKs, bool usePS = false);
		
		std::complex<double> eval(double s12) const override;
	protected:
		bool   _usePS;
		double _mPi;
		double _mKs;
	
};

class BELLE_LASS : public massShape {
	public:
		BELLE_LASS(const std::vector<double>& parameters, double mPi, double mK);

		std::complex<double> eval(double s12) const override;

	protected:
		double _mPi;
		double _mK;

};

class BELLE_KMatrix : public massShape {
	public:
		BELLE_KMatrix(const std::vector<double>& parameters);

		std::complex<double> eval(double s12) const override;
};
#endif//MASSSHAPE__
