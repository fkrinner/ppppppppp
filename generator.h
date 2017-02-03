#ifndef GENERATOR__
#define GENERATOR__
#include<string>
#include<vector>
#include"kinematicSignature.h"
class generator {
	public:
		generator();

		virtual std::vector<double> generate() const;

		size_t nKin() const {return _kinSignature.nKin();};
		kinematicSignature kinSignature() const {return _kinSignature;}

		bool setMaxFail(size_t n) {_maxFail = n; return true;};
	protected:
		kinematicSignature _kinSignature;
		size_t             _maxFail;
		mutable size_t     _failCount;
};

class threeParticleMassGenerator : public generator {
	public:
		threeParticleMassGenerator(double initialMass,const std::vector<double>& fsMasses);
		std::vector<double> generate() const override;
		bool isValidPoint(const std::vector<double>& kin) const;
	protected:
		double              _initialMass;
		std::vector<double> _fsMasses;
		
};
#endif//GENERATOR__
