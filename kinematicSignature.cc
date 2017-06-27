#include"kinematicSignature.h"
#include<iostream>

kinematicSignature::kinematicSignature(size_t identifier) : _identifier(identifier) {};

bool kinematicSignature::operator==(const kinematicSignature& other) const {
	return _identifier == other._identifier;
}

std::vector<std::vector<double> > kinematicSignature::getBoseSymmetrizedKinematics(const std::vector<double>& kin) const {
	if (kin.size() != nKin()) {
		std::cerr << "kinematicSignature::getBoseSymmetrizedKinematics(...): ERROR: Number of kinematic variables does not match. Returning empty vetor." << std::endl;
		return std::vector<std::vector<double> >();
	}
	switch(_identifier) {
		case 0: {return std::vector<std::vector<double> >();}
		case 1: {std::vector<std::vector<double> > retVal(2); retVal[0] = {kin[0], kin[1], kin[2]}; retVal[1] = {kin[0], kin[2], kin[1]}; return retVal;}
		case 2: {std::vector<std::vector<double> > retVal(1); retVal[0] = {kin[0], kin[1], kin[2]}; return retVal;}
	}
	std::cerr << "kinematicSignature::getBoseSymmetrizedKinematics(): ERROR: Unknown identifier: " << _identifier << ". Retruning empty vector" <<  std::endl;
	return std::vector<std::vector<double> >();
}

size_t kinematicSignature::nKin() const {
	switch (_identifier) {
		case 0: return 0;
		case 1: return 3;
		case 2: return 3;
	}
	std::cerr << "kinematicSignature::nKin(): ERROR: Unknown identifier: " << _identifier << std::endl;
	return 0;
}

std::vector<size_t> kinematicSignature::isobarMassIndices() const { 
	switch (_identifier) {
		case 0: return std::vector<size_t>();
		case 1: return std::vector<size_t>(1,1);
		case 2: return std::vector<size_t>();
	}
	std::cerr << "kinematicSignature::nKin(): ERROR: Unknown identifier: " << _identifier << std::endl;
	return std::vector<size_t>();
}

void kinematicSignature::print() const {
	switch (_identifier) {
		case 0: std::cout << "Empty kinematic signature" << std::endl; return;
		case 1: std::cout << "Signature for one initial state going into three final state particles.\nThe kinematic variables are {s, s_{12}, s_{13}\nParticles 2 and 3 are assumed to be identical" <<  std::endl; return;
		case 2: std::cout << "Signature for one initial state going into three final state particles.\nThe kinematic variables are {s, s_{12}, s_{13}\nNo particles are identical" << std::endl; return;
	}
	std::cerr << "kinematicSignature::print(): ERROR: Unknown identifier: " << _identifier << std::endl;
}
