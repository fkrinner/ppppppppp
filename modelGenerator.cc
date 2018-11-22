#include"modelGenerator.h"
#include<iostream>

modelGenerator::modelGenerator(std::shared_ptr<amplitude> model, std::shared_ptr<generator> gen, std::shared_ptr<efficiencyFunction> efficiency) : 
	modelGenerator::modelGenerator(std::vector<std::shared_ptr<amplitude> >(1, model), gen, efficiency) {};

modelGenerator::modelGenerator(std::vector<std::shared_ptr<amplitude> > model, std::shared_ptr<generator> gen, std::shared_ptr<efficiencyFunction> efficiency):
	_maxFail(1000), _failCount(0), _kinSignature(gen->kinSignature()), _generator(gen), _model(model), _efficiency(efficiency) {
	if (model.size() == 0) {
		std::cout << "modelGenerator::modelGenerator(...): ERROR: No amplitude given" << std::endl;
		throw;
	}
	if (not(*_kinSignature == *(_model[0]->kinSignature()))) {
		std::cout << "modelGenerator::modelGenerator(...): ERROR: Kinematic signatures does not match with generator" << std::endl;
		throw;
	}
	if (_efficiency) {
		if (not(*_kinSignature == *(_efficiency->kinSignature()))) {
			std::cout << "modelGenerator::modelGenerator(...): ERROR: Kinematic signatures does not match with efficiency" << std::endl;
			throw;
		}
	}
}

std::pair<std::pair<size_t, double>, std::vector<std::vector<double > > > modelGenerator::burnIn(size_t nPoints) const {
	double maxWeight = 0.;
	std::vector<double> weights(nPoints);
	std::vector<std::vector<double> > points(nPoints);
	for (size_t p = 0; p < nPoints; ++p) {
		std::vector<double> kin = _generator->generate();
		double intens = 0;
		for (const std::shared_ptr<amplitude> ampl: _model) {
			intens += ampl->intens(kin);
		}
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
		double intens = 0.;
		for (std::shared_ptr<amplitude> ampl : _model) {
			intens += ampl->intens(kin);
		}
		if (intens > utils::random() * maxWeight) {
			retVal[found] = kin;
			++found;
			_failCount = 0;
		} else {
			++_failCount;
			if (_failCount > _maxFail) {
				std::cout << "modelGenerator::generateDataPoints(...): Over " << _maxFail << " failed attempts. Aborting." << std::endl;
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
	double intens = 0.;
	for (std::shared_ptr<amplitude> ampl : _model) {
		intens += ampl->intens(kin);
	}
	if (_efficiency) {
		intens *= _efficiency->eval(kin);
	}
	return std::pair<double, std::vector<double> > (intens, kin);
}
