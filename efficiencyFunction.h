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

		virtual double call(const std::vector<double>& kin)          const;
	protected:
		std::shared_ptr<kinematicSignature> _kinSignature;
};

class threeParticlPerfectEfficiency : public efficiencyFunction {
	public:
		threeParticlPerfectEfficiency (std::shared_ptr<kinematicSignature> kinSig = std::make_shared<kinematicSignature>(1));

		virtual double call(const std::vector<double>& kin)          const;
		
};

class BELLE_DtoKpipi_efficiency : public efficiencyFunction {
	public: // Kineamtic variable are {m_D^2, m_{Kpi(RS}^2, m_{pipi}^2}
		BELLE_DtoKpipi_efficiency ();

		virtual double call(const std::vector<double>& kin)          const;
};
#endif//EFFICIENCYFUNCTION__
