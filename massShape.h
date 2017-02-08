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
#endif//MASSSHAPE__
