#ifndef ANGULARDEPENDENCE__
#define ANGULARDEPENDENCE__
#include<string>
#include<vector>
#include<complex>
#include<memory>
#include"kinematicSignature.h"

class angularDependence {
	public:
		angularDependence (std::shared_ptr<kinematicSignature> kinSignature, std::string name);

		virtual std::complex<double> eval (const std::vector<double>& kin) const;

		std::string                         name()         const {return _name;}	
		size_t                              nKin()         const {return _kinSignature->nKin();}
		std::shared_ptr<kinematicSignature> kinSignature() const {return _kinSignature;}
	protected:
		std::shared_ptr<kinematicSignature> _kinSignature;
		std::string                         _name;
};

class sameMassZeroS : public angularDependence {
	public:
		sameMassZeroS (double fsMass);

		std::complex<double> eval (const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};

class sameMassOneP : public angularDependence {
	public:
		sameMassOneP (double fsMass);

		std::complex<double> eval (const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};

class sameMassTwoD : public angularDependence {
	public:
		sameMassTwoD (double fsMass);

		std::complex<double> eval (const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};

class sameMassZeroSnonRelativistic : public angularDependence {
	public:
		sameMassZeroSnonRelativistic (double fsMass);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};

class sameMassOnePnonRelativistic : public angularDependence {
	public:
		sameMassOnePnonRelativistic (double fsMass);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};

class sameMassTwoDnonRelativistic : public angularDependence {
	public:
		sameMassTwoDnonRelativistic (double fsMass);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		double _fsMass;
};

class arbitraryMass_S_nonRelativistic : public angularDependence {
	public:
		arbitraryMass_S_nonRelativistic(size_t isobarIndex, std::vector<double> fsMasses);

		std::complex<double> eval(const std::vector<double>& kin) const override;

	private:
		size_t              _isobarIndex;
		std::vector<double> _fsMasses;
};

class arbitraryMass_P_nonRelativistic : public angularDependence {
	public:
		arbitraryMass_P_nonRelativistic(size_t isobarIndex, std::vector<double> fsMasses);

		std::complex<double> eval(const std::vector<double>& kin) const override;

	private:
		size_t              _isobarIndex;
		std::vector<double> _fsMasses;
};

class ratioOfDependences : public angularDependence {
	public:
		ratioOfDependences (std::shared_ptr<angularDependence> numerator, std::shared_ptr<angularDependence> denominator);

		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		std::shared_ptr<angularDependence> _numerator;
		std::shared_ptr<angularDependence> _denominator;
};

class BELLE_S : public angularDependence {
	public: 
		BELLE_S(size_t isobarIndex);
		std::complex<double> eval(const std::vector<double>& kin) const override;
	private:
		size_t _isobarIndex;
};

class BELLE_P : public angularDependence {
	public: 
		BELLE_P(size_t isobarIndex);
		std::complex<double> eval(const std::vector<double>& kin) const override;

		bool setFSmasses(const std::vector<double> newMasses);
	private:
		size_t _isobarIndex;
		std::vector<double> _fsMasses;
};

class BELLE_D : public angularDependence {
	public: 
		BELLE_D(size_t isobarIndex);
		std::complex<double> eval(const std::vector<double>& kin) const override;

		bool setFSmasses(const std::vector<double> newMasses);
	private:
		size_t _isobarIndex;
		std::vector<double> _fsMasses;
};
#endif// ANGULARDEPENDENCE__
