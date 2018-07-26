#ifndef MODELAMPLITUDE__
#define MODELAMPLITUDE__
#include<vector>
#include<complex>
#include<string>
#include<memory>
#include"amplitude.h"
#include"kinematicSignature.h"
class modelAmplitude : public amplitude {
	public:
		modelAmplitude (std::vector<std::complex<double> > transitionAmplitudes, std::vector<std::shared_ptr<amplitude> > amplitudes, std::vector<double> normalizations, std::string name);
		std::complex<double> eval   (const std::vector<double>& kin) const override;

		bool                                setTransitionAmplitude (size_t n, std::complex<double> amp);
	protected:
		size_t                                   _nAmpl;
		std::vector<std::complex<double> >       _transitionAmplitudes;
		std::vector<std::shared_ptr<amplitude> > _amplitudes;
		std::vector<double>                      _normalizations;
};
#endif// MODELAMPLITUDE__
