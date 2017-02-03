#ifndef KINEMATICSIGNATURE__
#define KINEMATICSIGNATURE__
#include<string>
#include<vector>
class kinematicSignature {
	public:
		kinematicSignature() : kinematicSignature(0) {};
		kinematicSignature(size_t identifier);
		
		bool operator==(const kinematicSignature& other) const;

		std::vector<std::vector<double> > getBoseSymmetrizedKinematics(const std::vector<double>& kin) const;

		size_t              nKin() const;
		std::vector<size_t> isobarMassIndices() const;
		
		void print() const;		
	protected:
		size_t _identifier;

};
#endif// KINEMATICSIGNATURE__
