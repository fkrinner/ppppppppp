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

		std::string                         name         () const {return _name;}
		size_t                              nKin         () const {return _kinSignature->nKin();}
		std::shared_ptr<kinematicSignature> kinSignature () const { return _kinSignature;}
		

	protected:
		std::shared_ptr<kinematicSignature> _kinSignature;
		std::string                         _name;

};

class threeParticleIsobaricAmplitude : public amplitude {
	public:
		threeParticleIsobaricAmplitude(bool boseSymmetrize, std::string name, std::shared_ptr<massShape> shape, std::shared_ptr<angularDependence> angDep);
		
		std::complex<double> eval(const std::vector<double>& kin) const;
		std::complex<double> evalSingleSymmTerm(const std::vector<double>& kin) const;

	private:
		bool                               _bose;
		std::shared_ptr<massShape>         _massShape;
		std::shared_ptr<angularDependence> _angDep;
};
#endif//AMPLITUDE__
