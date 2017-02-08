#ifndef MODELAMPLITUDE__
#define MODELAMPLITUDE__
#include<vector>
#include<complex>
#include<string>
#include<memory>
#include"amplitude.h"
#include"kinematicSignature.h"
class modelAmplitude {
	public:
		modelAmplitude (std::vector<std::complex<double> > transitionAmplitudes, std::vector<std::shared_ptr<amplitude> > amplitudes, std::vector<double> normalizations);

		double               intens (const std::vector<double>& kin) const {return std::norm(ampl(kin));}
		std::complex<double> ampl   (const std::vector<double>& kin) const;

		std::shared_ptr<kinematicSignature> kinSignature           ()                                   const {return _kinSignature;}
		bool                                setTransitionAmplitude (size_t n, std::complex<double> amp);
	protected:
		size_t                                   _nAmpl;
		std::shared_ptr<kinematicSignature>      _kinSignature;
		std::vector<std::complex<double> >       _transitionAmplitudes;
		std::vector<std::shared_ptr<amplitude> > _amplitudes;
		std::vector<double>                      _normalizations;
};
#endif// MODELAMPLITUDE__
