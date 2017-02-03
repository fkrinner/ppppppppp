#ifndef ANGULARDEPENDENCE__
#define ANGULARDEPENDENCE__
#include<string>
#include<vector>
#include<complex>
#include"kinematicSignature.h"

class angularDependence {
	public:
		angularDependence(kinematicSignature kinSignature, std::string name);

		virtual std::complex<double> eval(const std::vector<double>& kin) const;

		std::string name() const {return _name;}	
		size_t      nKin() const {return _kinSignature.nKin();}
		kinematicSignature kinSignature() const {return _kinSignature;}
	protected:
		kinematicSignature _kinSignature;
		std::string        _name;
};

class sameMassZeroS : public angularDependence {
	public:
		sameMassZeroS(double fsMass);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};

class sameMassOneP : public angularDependence {
	public:
		sameMassOneP(double fsMass);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};

class sameMassTwoD : public angularDependence {
	public:
		sameMassTwoD(double fsMass);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};
#endif// ANGULARDEPENDENCE__
