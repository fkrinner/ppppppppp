#ifndef EFFICIENCYFUNCTION__
#define EFFICIENCYFUNCTION__
#include<string>
#include<vector>
#include<memory>
#include"kinematicSignature.h"
class efficiencyFunction {
	public:
		efficiencyFunction ();

		size_t                              nKin        ()         const {return _kinSignature->nKin();};
		std::shared_ptr<kinematicSignature> kinSignature()         const {return _kinSignature;}

		virtual double eval(const std::vector<double>& kin) const;
	protected:
		std::shared_ptr<kinematicSignature> _kinSignature;
};

class threeParticlPerfectEfficiency : public efficiencyFunction {
	public:
		threeParticlPerfectEfficiency (std::shared_ptr<kinematicSignature> kinSig = std::make_shared<kinematicSignature>(1));

		double eval(const std::vector<double>& kin) const override;
		
};

class BELLE_DtoKpipi_efficiency : public efficiencyFunction {
	public: // Kineamtic variable are {m_D^2, m_{Kpi(RS}^2, m_{pipi}^2}
		BELLE_DtoKpipi_efficiency ();

		double eval(const std::vector<double>& kin) const override;

		bool setMpipiCosTCut(const double minM2Pi, const double maxAbsCosT) {_minM2Pisquared=minM2Pi*minM2Pi; _maxAbsCosT = maxAbsCosT; return true;}
	private:
		double _minM2Pisquared;  // events with m2Pi <  _minM2Pi will be cut, if their |cosT| is bigger than _maxAbsCost (store the square, however)
		double _maxAbsCosT;      // events with bigger |cosT| > _maxAbsCost will be cut, if their m2Pi is smaller than _minM2Pi
		
};
#endif//EFFICIENCYFUNCTION__
