#ifndef MODELGENERATOR__
#define MODELGENERATOR__
#include<memory>
#include<vector>
#include"modelAmplitude.h"
#include"generator.h"
#include"kinematicSignature.h"
#include"utils.h"
class modelGenerator {
	public:
		modelGenerator(std::shared_ptr<modelAmplitude> modelAmplitude, std::shared_ptr<generator> generator);

		std::pair<std::pair<size_t, double>, std::vector<std::vector<double > > > burnIn(size_t nPoints) const;
		std::vector<std::vector<double> > generateDataPoints(size_t nPoints, size_t nBurnIn = 1000) const;

		std::pair<double, std::vector<double> > getSinglePoint() const;

		bool setMaxFail(size_t n) {_maxFail = n; return true;};
	protected:
		size_t                          _maxFail;
		mutable size_t                  _failCount;
		kinematicSignature              _kinSignature;
		std::shared_ptr<generator>      _generator;
		std::shared_ptr<modelAmplitude> _model;
};
#endif// MODELGENERATOR__
