#include"modelGenerator.h"
#include<iostream>
modelGenerator::modelGenerator(std::shared_ptr<amplitude> model, std::shared_ptr<generator> gen):
	_maxFail(1000), _failCount(0), _kinSignature(model->kinSignature()), _generator(gen), _model(model) {
	if (not(*_kinSignature == *(_generator->kinSignature()))) {
		std::cerr << "modelGenerator::modelGenerator(...): ERROR: Kinematic signatures do not match" << std::endl;
	}
}

std::pair<std::pair<size_t, double>, std::vector<std::vector<double > > > modelGenerator::burnIn(size_t nPoints) const {
	double maxWeight = 0.;
	std::vector<double> weights(nPoints);
	std::vector<std::vector<double> > points(nPoints);
	for (size_t p = 0; p < nPoints; ++p) {
		std::vector<double> kin = _generator->generate();
		double intens = _model->intens(kin);
		weights[p] = intens;
		points[p]  = kin;
		if (intens > maxWeight) {
			maxWeight = intens;
		}
	}
	std::vector<std::vector<double> > deWeighted(nPoints);
	size_t remaining = 0;
	for (size_t p = 0; p < nPoints; ++p) {
		if (weights[p] > utils::random() * maxWeight) {
			deWeighted[remaining] = points[p];
			++remaining;
		}
	}
	return std::pair<std::pair<size_t, double>, std::vector<std::vector<double> > >(std::pair<size_t, double>(remaining, maxWeight), deWeighted);
}

std::vector<std::vector<double> > modelGenerator::generateDataPoints(size_t nPoints, size_t nBurnIn) const {
	double maxWeight = 0.;
	size_t found     = 0;
	std::vector<std::vector<double> > retVal(nPoints);
	{
		std::pair<std::pair<size_t, double>, std::vector<std::vector<double > > > burn = burnIn(nBurnIn);
		found     = std::min(burn.first.first, nPoints); 
		maxWeight = burn.first.second;
		for (size_t i = 0; i < found; ++i) {
			retVal[i] = burn.second[i];
		}
	}
	while (found < nPoints) {
		std::vector<double> kin = _generator->generate();
		double intens = _model->intens(kin);
		if (intens > utils::random() * maxWeight) {
			retVal[found] = kin;
			++found;
			_failCount = 0;
		} else {
			++_failCount;
			if (_failCount > _maxFail) {
				std::cerr << "modelGenerator::generateDataPoints(...): Over " << _maxFail << " failed attempts. Aborting." << std::endl;
				throw;
			}
		}
		if (intens > maxWeight) {
			std::vector<std::vector<double> > newVal(nPoints);
			size_t count = 0;
			for (size_t i = 0; i < found; ++i) {
				if (maxWeight > utils::random() * intens) {
					newVal[count] = retVal[i];
					++count;
				}
			}
//			std::cout << (double)count/found << " vs " << maxWeight/intens << std::endl;
			found     = count;
			maxWeight = intens;
			retVal    = newVal;
		}		
	}
	return retVal;
}

std::pair<double, std::vector<double> > modelGenerator::getSinglePoint() const { 
	std::vector<double> kin = _generator->generate();
	double intens = _model->intens(kin);
	return std::pair<double, std::vector<double> > (intens, kin);
}
