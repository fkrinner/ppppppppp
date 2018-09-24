#include"massShape.h"
#include"utils.h"
#include<iostream>

massShape::massShape(std::string name, std::vector<double> pars, std::vector<std::string> parNames):
	_nPar(pars.size()), _name(name), _parameters(pars) {
	if (parNames.size() == _nPar) {
		_parameterNames = parNames;
	} else {
		if (parNames.size() > 0) {
			std::cerr << "massShape::massShape(...): ERROR: Number of parameter names does not match and is not zero." << std::endl;
		}
		_parameterNames = std::vector<std::string>();
		for (size_t i = 0; i < _nPar; ++i) {
			_parameterNames.push_back(std::string("unname_parameter_")+std::to_string(i));
		}
	}
}

std::complex<double> massShape::eval(double m) const {
	std::cerr << "massShape::eval(double): ERROR: Evaluating the bas class at mass m = " << m << ", returning zero." << std::endl;
	return std::complex<double>(0.,0.);
}

bool massShape::setParameter(size_t n, double val) {
	if (n < _nPar) {
		_parameters[n] = val;
		return true;
	}
	return false;
}
	
std::pair<bool, double> massShape::getParameter(size_t n) const {
	if (n<_nPar){
		return std::pair<bool, double>(true, _parameters[n]);	
	}
	return std::pair<bool, double>(false, 0.);
}

bool massShape::setParName(size_t n, std::string name) {
	if(n<_nPar) {
		_parameterNames[n] = name;
		return true;
	}
	return false;
}

std::pair<bool, std::string> massShape::getParName(size_t n) const {
	if (n<_nPar) {
		return std::pair<bool, std::string>(true, _parameterNames[n]);
	}
	return std::pair<bool, std::string>(false, "");
}
//------------------------------------------------------------------------------
simpleBW::simpleBW(double mass, double width):
	massShape("simpleBW", {mass, width}, {"m0","G0"}) {}

std::complex<double> simpleBW::eval(double s) const {
	double num = _parameters[0] * _parameters[1];
	std::complex<double> den = std::complex<double>(_parameters[0]*_parameters[0] - s, -_parameters[0]*_parameters[1]);
	return num/den;
}
//------------------------------------------------------------------------------
stepLike::stepLike(double sMin, double sMax) :
	massShape("stepLike", {sMin, sMax}, {"sMin","sMax"}) {}

std::complex<double> stepLike::eval(double s) const {
	if (s > _parameters[0] and s <= _parameters[1]) {
		return std::complex<double>(1.,0.);
	}
	return std::complex<double>(0.,0.);
}
//------------------------------------------------------------------------------
constant::constant() :
	massShape("constant", {}, {}) {};

std::complex<double> constant::eval(double s) const {
	return std::complex<double>(s/s,0.);
}
//------------------------------------------------------------------------------
zeroMode0pp::zeroMode0pp(double s, double m2) :
	massShape("zeroMode0pp", {s, m2}, {"s", "m2"}) {};

std::complex<double> zeroMode0pp::eval(double s12) const {
	double retVal = (_parameters[0] - 3*s12 + 3*_parameters[1])/8;
//	std::cout << "Called zeroMode0pp with " << s12 << " giving " << retVal << std::endl;
	return std::complex<double>(retVal, 0.);
}
//------------------------------------------------------------------------------
zeroMode1mm::zeroMode1mm(double s, double m2) :
	massShape("zeroMode0pp", {s, m2}, {"s", "m2"}) {};

std::complex<double> zeroMode1mm::eval(double s12) const {
	double retVal = _parameters[0]/(_parameters[0] + s12 - _parameters[1]);
//	std::cout << "Called zeroMode1mm with " << s12 << " giving " << retVal << std::endl;
	return std::complex<double>(retVal, 0.);
}
//------------------------------------------------------------------------------
polynomialMassShape::polynomialMassShape(std::vector<std::complex<double> > coefficients, double baseExponent) :
	massShape(std::string("polynomialMassShape_deg") + std::to_string(coefficients.size()-1), {}, {}), _polDeg(0), _baseExponent(baseExponent) {

	if (!utils::checkComplexDouble()) {
		std::cout << "polynomialMassShape::polynomialMassShape(...): ERROR: std::complex<double>* is not double real, double imag array" << std::endl;
		throw;
	}
	for (size_t c = 0; c < coefficients.size(); ++c) {
		_parameters.push_back(coefficients[c].real());
		_parameterNames.push_back(std::string("c_") + std::to_string(c) + std::string("_r"));
		_parameters.push_back(coefficients[c].imag());
		_parameterNames.push_back(std::string("c_") + std::to_string(c) + std::string("_i"));
	}
	_nPar   = 2*coefficients.size();
	_polDeg = coefficients.size();
}

std::complex<double> polynomialMassShape::eval(double s12) const {
	std::complex<double> retVal(0.,0.);
	std::complex<double>* coeffs = (std::complex<double>*)&_parameters[0];
	for (size_t c = 0; c < _polDeg; ++c) {
		retVal += coeffs[c] * pow(s12, _baseExponent*c);
	}
	return retVal;
}
//------------------------------------------------------------------------------
BELLEbreitWigner::BELLEbreitWigner(std::string name, double mass, double width, size_t spin, double motherMass, double bachelorMass, double daughterMass1, double daughterMass2) :
	massShape("BELLEbw_"+name, {mass, width}, {name+"_mass", name+"_width"}), _spin(spin), _motherMass(motherMass), _bachelorMass(bachelorMass), _daughterMass1(daughterMass1), _daughterMass2(daughterMass2), _Rr(1.5), _RD(5.)
{
	if (_motherMass - _bachelorMass < _parameters[0]) {
		std::cout  << "BELLEbreitWigner::BELLEbreitWigner(...): ERROR: On shell resonance mass of '" << _name << "' too heavy for decay of mother particle: " << _motherMass << " -> " << _parameters[0] << " + " << _bachelorMass << std::endl;
		throw;
	}
	if (_daughterMass1 + _daughterMass2 > _parameters[0]) {
		std::cout << "BELLEbreitWigner::BELLEbreitWigner(...): ERROR: On shell resonance mass of '" << _name << "' too light for decay into daughter particles: " << _parameters[0] << " -> " << _daughterMass1 << " + " << _daughterMass2 << std::endl;
		throw;
	}	

	if (_spin > 2) {
		std::cout << "BELLEbreitWigner::BELLEbreitWigner(...): ERROR: Spin > 2 not supportet yet" << std::endl;
		throw;
	}
}

std::complex<double> BELLEbreitWigner::eval(double s12) const {

	const double m12   = pow(s12, .5);

	if ((m12 < _daughterMass1 + _daughterMass2) || (m12 > _motherMass - _bachelorMass)) {
		return std::complex<double>(0.,0.);
	}

	const double S          = _motherMass*_motherMass;
	const double sr         = _parameters[0]*_parameters[0];
	const double sDaughter1 = _daughterMass1*_daughterMass1;
	const double sDaughter2 = _daughterMass2*_daughterMass2;
	const double sBatch     = _bachelorMass*_bachelorMass;

	const double pr  = pow(pow(sr - sDaughter1 - sDaughter2, 2) - 4*sDaughter1*sDaughter2, .5)/2/_parameters[0];
	const double pAB = pow(pow(s12 - sDaughter1 - sDaughter2, 2) - 4*sDaughter1*sDaughter2, .5)/2/m12;

	const double pD   = pow(pow(S - sr - sBatch, 2) - 4*sr*sBatch, .5)/2/_motherMass;
	const double pABC = pow(pow(S - s12 - sBatch, 2) - 4*s12*sBatch, .5)/2/_motherMass;

	double Fr = 1.;
	double FD = 1.;
	if (_spin == 1) {
		Fr = pow((pow(_Rr*pr,2)+1)/(pow(_Rr*pAB,2)+1),.5);
		FD = pow((pow(_RD*pD,2)+1)/(pow(_RD*pABC,2)+1),.5);
	} else if (_spin == 2) {
		const double xr =  _Rr*_Rr*pr *pr;
		const double xAB = _Rr*_Rr*pAB*pAB;
		Fr = pow((pow(xr-3.,2) + 9*xr)/(pow(xAB-3.,2) + 9*xAB),.5);

		const double xD   = _RD*_RD*pD  *pD;
		const double xABC = _RD*_RD*pABC*pABC;
		FD = pow((pow(xD-3.,2) + 9*xD)/(pow(xABC-3.,2) + 9*xABC),.5);
	}
	
	const double Gamma = _parameters[1]* _parameters[0]/m12 * Fr*Fr * pow(pAB/pr, 2*_spin+1);
	
	std::complex<double> retVal = std::complex<double>(Fr*FD,0.)/std::complex<double>(_parameters[0]*_parameters[0] - s12, - _parameters[0]*Gamma);
//	if (std::isnan(retVal.real()) || std::isnan(retVal.imag())) {
//		std::cout << "s12 = " << s12 << "; sDaughter1 = " << sDaughter1 << "; sDaughter2 = " << sDaughter2 << std::endl;
//		std::cout << "pr = " << pr << "; pAB = " << pAB << std::endl;
//		std::cout << "Fr = " << Fr << "; FD = " << FD << "; Gamma = " << Gamma << std::endl;
//	}
	return retVal;
}
//------------------------------------------------------------------------------
BELLE_LASS_KpiS::BELLE_LASS_KpiS(const std::vector<double>& parameters, double mPi, double mKs, bool usePS) : massShape("BELLE_LASS_KpiS", parameters, {"a", "r", "M0", "G0", "phiF", "phiR", "phiRsin", "F", "R", "MMax"}), _usePS(usePS), _mPi(mPi), _mKs(mKs) {}

std::complex<double> BELLE_LASS_KpiS::eval(double s12) const {
	double M  = pow(s12, .5);
	std::complex<double> retVal(0.,0.);
	if (M > _mPi + _mKs && M < _parameters[9]){
		const double q  = pow(utils::breakupMomentumSquared(s12, _mKs*_mKs, _mPi*_mPi), .5);
		const double q0 = pow(utils::breakupMomentumSquared(_parameters[2]*_parameters[2], _mKs*_mKs, _mPi*_mPi), .5);

		const double G = _parameters[3] * _parameters[2]/M * q/q0;
		const double deltaR = atan(_parameters[2] * G/(_parameters[2]*_parameters[2] - s12));
		const double deltaF = atan(2. * _parameters[0] * q / (2. + _parameters[0] * _parameters[1] * q * q));

		retVal += _parameters[7] * sin(deltaF + _parameters[4]) * exp(std::complex<double>(0., deltaF + _parameters[4])) + 
		          _parameters[8] * sin(deltaR + _parameters[6]) * exp(std::complex<double>(0., deltaR + _parameters[5])) * exp(std::complex<double>(0., 2* ( deltaF + _parameters[4])));
		if (_usePS) {
			retVal *= .5*M/q;
		}

/* Copied the code form LAURA++, but this version has much less parameters... check again when Stefan is back
		const double rVal     = _parameters[1];
		const double aVal     = _parameters[0];
		const double resMass  = _parameters[2];
		const double resWidth = _parameters[3];
		const double qRatio = q/q0;
		const double totWidth = resWidth*qRatio*(resMass/M);
		const double massSqTerm = resMass*resMass - s12;
		std::complex<double> altern(massSqTerm, resMass*totWidth);
		altern *= (resMass*resMass*resWidth/q0)/(massSqTerm*massSqTerm + resMass*resMass*totWidth*totWidth);

		const double tandeltaB = (2.0*aVal*q)/(2.0 + aVal*rVal*q*q);
		const double tanSq = tandeltaB*tandeltaB;
		const double cos2PhaseShift = (1.0 - tanSq) / (1.0 + tanSq);
		const double sin2PhaseShift = 2.0*tandeltaB / (1.0 + tanSq);
		std::complex<double> phaseShift(cos2PhaseShift, sin2PhaseShift);

		altern *= phaseShift;
		const double qcotdeltaB = 1.0/aVal + (rVal*q*q)/2.0;
		std::complex<double> bkgAmplitude(qcotdeltaB, q);
		bkgAmplitude *= M/(qcotdeltaB*qcotdeltaB + q*q);

		altern += bkgAmplitude;

		std::cout << retVal << " " << altern << " but check me again, when Stafen is here, since LAURA++ does not use phiFR, FRm" << std::endl;
*/
	}
	return retVal;
}
//------------------------------------------------------------------------------
BELLE_LASS::BELLE_LASS(const std::vector<double>& parameters, double mPi, double mK) : massShape("BELLE_LASS", parameters, {"a", "r", "M0", "G0", "phiF", "phiR", "F", "R"}), _mPi(mPi), _mK(mK) {}

std::complex<double> BELLE_LASS::eval(double s12) const {

//	double pi180inv = 1.0/EvtConst::radToDegrees;
//	double ampl =_amp;
//	double theta=_phase;

//	double s = mab2;

	const double _a    = _parameters[0]; //_LASS_a;
	const double _r    = _parameters[1]; //_LASS_r;
	const double _R    = _parameters[7]; //_LASS_R;
	const double _phiR = _parameters[5]; //_LASS_phi_R;
	const double _F    = _parameters[6]; //_LASS_F;
	const double _phiF = _parameters[4]; //_LASS_phi_F;

	// T = R sin(deltaR) exp[i (deltaR+phiR)] exp[i2 (deltaB+phiB)] + B sin(deltaB+phiB) exp[i (deltaB+phiB)]
	//   = R exp[i (phiR+2phiB)] sin(deltaR) exp[i deltaR] exp[i2 deltaB] + B exp[i phiB] { cos(phiB) + cot(deltaB) sin(phiB) } sin(deltaB) exp[i deltaB]
	//   = R exp[i (phiR+2phiB) m0*Gamma(m)/(m0*m0-m*m-i*m0*Gamma(m)) exp[i2 deltaB] + B exp[i phiB] { cos(phiB) + cot(deltaB) sin(phiB) m/[q*cot(deltaB)-i*q]
	// The propagator is T/rho, where rho = 2 q/sqrt(s) is the two-body phase space

//	const double gamma = _parameters[3]; //_gammaR;
//	const double bwm   = _parameters[2]; //_massR;

	const double mR     = _parameters[2]; //_massR;
	const double gammaR = _parameters[3]; //_gammaR;

	double fR=1.0; // K*0(1430) has spin zero
	int power=1; // power is 1 for spin zero

	const double mAB=pow(s12,.5); //   (_p4_d1+_p4_d2).mass();

	const double mA=_mPi; //_p4_d1.mass();
	const double mB=_mK; //_p4_d2.mass();

	const double pAB=pow( (((mAB*mAB-mA*mA-mB*mB)*(mAB*mAB-mA*mA-mB*mB)/4.0) -
			  mA*mA*mB*mB)/(mAB*mAB),.5);
	const double q=pAB;

	const double pR=pow( (((mR*mR-mA*mA-mB*mB)*(mR*mR-mA*mA-mB*mB)/4.0) -
			  mA*mA*mB*mB)/(mR*mR),.5);

	if (isnan(pAB)) {
		return std::complex<double>(0.,0.);
	}

	// compute running width g
	const double g = gammaR*pow(pAB/pR,power)*(mR/mAB)*fR*fR;

	const std::complex<double> propagator_relativistic_BreitWigner = 1./std::complex<double>(mR*mR - mAB*mAB,-mR*g);

	// non-resonant phase shift
	const double cot_deltaF = 1.0/(_a*q) + 0.5*_r*q;
	const double qcot_deltaF = 1.0/_a + 0.5*_r*q*q;

	// calculate resonant part
	const std::complex<double> expi2deltaF = std::complex<double>(qcot_deltaF, q)/ std::complex<double>(qcot_deltaF, -q);

	const std::complex<double> resonant_term_T = _R * std::complex<double>(cos(_phiR + 2 * _phiF), sin(_phiR + 2 * _phiF)) * propagator_relativistic_BreitWigner * mR * gammaR * mR / pR * expi2deltaF;

	// calculate non-resonant part
	const std::complex<double>  non_resonant_term_F = _F * std::complex<double>(cos(_phiF), sin(_phiF)) * (cos(_phiF) + cot_deltaF * sin(_phiF)) * sqrt(s12) / std::complex<double>(qcot_deltaF, -q);

	// sum up non-resonant and resonant terms
	const std::complex<double> LASS_contribution = non_resonant_term_F + resonant_term_T;

//	// convert std::complex<double> to TComplex
//	TComplex LASS_contribution_TComplex (LASS_contribution._rpart, LASS_contribution._ipart);

//	TComplex matrixEl = ampl * TComplex(cos(theta*pi180inv), sin(theta*pi180inv)) * LASS_contribution_TComplex;
	return LASS_contribution;
}
//------------------------------------------------------------------------------
BELLE_KMatrix::BELLE_KMatrix(const std::vector<double>& parameters) : massShape("BELLE_KMatrix", parameters, 
	{"Kmatrix_beta1_Amplitude","Kmatrix_beta1_Phase",
	 "Kmatrix_beta2_Amplitude","Kmatrix_beta2_Phase",
	 "Kmatrix_beta3_Amplitude","Kmatrix_beta3_Phase",
	 "Kmatrix_beta4_Amplitude","Kmatrix_beta4_Phase",
	 "Kmatrix_beta5_Amplitude","Kmatrix_beta5_Phase",
	 "Kmatrix_f_prod_11_Amplitude", "Kmatrix_f_prod_11_Phase",
	 "Kmatrix_f_prod_12_Amplitude", "Kmatrix_f_prod_12_Phase",
	 "Kmatrix_f_prod_13_Amplitude", "Kmatrix_f_prod_13_Phase",
	 "Kmatrix_f_prod_14_Amplitude", "Kmatrix_f_prod_14_Phase",
	 "Kmatrix_f_prod_15_Amplitude", "Kmatrix_f_prod_15_Phase",
	 "Kmatrix_s_prod_0"}) {}

std::complex<double> BELLE_KMatrix::eval(double s12) const {
	double s = s12;

	//Define the complex coupling constants
	//The convention is as follow
	//i=0 --> pi+ pi-
	//i=1 --> KK
	//i=2 --> 4pi
	//i=3 --> eta eta
	//i=4 --> eta eta'
	std::vector<std::vector<double> > g(5,std::vector<double>(5,0.)); // Coupling constants. The first index is the pole index. The second index is the decay channel

	//pi+pi- channel
	g[0][0]=0.22889;
	g[1][0]=0.94128;
	g[2][0]=0.36856;
	g[3][0]=0.33650;
	g[4][0]=0.18171;
	//K+K- channel
	g[0][1]=-0.55377;
	g[1][1]=0.55095;
	g[2][1]=0.23888;
	g[3][1]=0.40907;
	g[4][1]=-0.17558;
	//4pi channel
	g[0][2]=0;
	g[1][2]=0;
	g[2][2]=0.55639;
	g[3][2]=0.85679;
	g[4][2]=-0.79658;
	//eta eta channel
	g[0][3]=-0.39899;
	g[1][3]=0.39065;
	g[2][3]=0.18340;
	g[3][3]=0.19906;
	g[4][3]=-0.00355;
	//eta eta' channel
	g[0][4]=-0.34639;
	g[1][4]=0.31503;
	g[2][4]=0.18681;
	g[3][4]=-0.00984;
	g[4][4]=0.22358;

	// Pole masses
	std::vector<double> ma(5,0.);   // Pole masses. The unit is in GeV

	ma[0]=0.651;
	ma[1]=1.20360;
	ma[2]=1.55817;
	ma[3]=1.21000;
	ma[4]=1.82206;

	// scattering data
	// double s_0_scatt=-3.92637;
	// double s_A=1.0;
	// double s_A0=-0.15;

	//Now define the K-matrix pole
	std::complex<double> n11,n12,n13,n14,n15,n21,n22,n23,n24,n25,n31,n32,n33,n34,n35,n41,n42,n43,n44,n45,n51,n52,n53,n54,n55;
	double  rho1sq,rho2sq,rho4sq,rho5sq;//,rho3sq
	std::complex<double> rho1,rho2,rho3,rho4,rho5;
	std::vector<std::complex<double> > rho(5,std::complex<double>(0.,0.));
	std::complex<double> pole,SVT,Alder;
	std::complex<double> det;
	std::vector<std::vector<std::complex<double> > >i(5,std::vector<std::complex<double> >(5, std::complex<double>(0.,0.)));  //inverse of the (I-iKp) matrix
	std::vector<std::vector<double> > f(5, std::vector<double>(5,0.));

	//Initalize the mass of the resonance
	double mpi   = 0.13957;
	double mK    = 0.493677;     //using charged K value
	double meta  = 0.54775;    //using PDG value
	double metap = 0.95778;   //using PDG value

	//Initialize the matrix to value zero
	std::vector<std::vector<std::complex<double> > > K(5, std::vector<std::complex<double> >(5, std::complex<double>(0.,0.)));
//	for(Int_t k=0;k<5;k++) {
//		for(Int_t l=0;l<5;l++) {
//			i[k][l]=std::complex<double>(0,0);
//			K[k][l]=std::complex<double>(0,0);
//			f[k][l]=0;
////			f_scatt[k][l]=0;
//		}
//		rho[k]=0;
//	}

	//Input the _f[i][j] scattering data
	double s_scatt=-3.92637;
	double sa=1.0;
	double sa_0=-0.15;

	f[0][0]=0.23399;  // f^scatt
	f[0][1]=0.15044;
	f[0][2]=-0.20545;
	f[0][3]=0.32825;
	f[0][4]=0.35412;

	f[1][0]=f[0][1];
	f[2][0]=f[0][2];
	f[3][0]=f[0][3];
	f[4][0]=f[0][4];

	//Construct the phase-space factor
	//For eta-eta' there is no difference term
	rho1sq=(1.0-(pow((mpi+mpi),2)/s));   //pi+ pi- phase factor
	if(rho1sq >=0) {
		rho1=std::complex<double>(sqrt(rho1sq),0);
	} else {
		rho1=std::complex<double>(0,sqrt(-rho1sq));
	}
	rho[0]=rho1;

	rho2sq=(1.0-(pow((mK+mK),2)/s));
	if(rho2sq >=0) {
		rho2=std::complex<double>(sqrt(rho2sq),0);
	} else {
		rho2=std::complex<double>(0,sqrt(-rho2sq));
	}
	rho[1]=rho2;
	//using the A&S 4pi phase space Factor:
	rho3=std::complex<double>(0,0);
	//Shit, not continue
	if(s<=1) {
		double real = 1.2274+0.00370909/(s*s) - (0.111203)/(s) - 6.39017*s +16.8358*s*s - 21.8845*s*s*s + 11.3153*s*s*s*s;
		double cont32=sqrt(1.0-(16.0*mpi*mpi));
		rho3=std::complex<double>(cont32*real,0);
	} else {
		rho3=std::complex<double>(sqrt(1.0-(16.0*mpi*mpi/s)),0);
	}
	rho[2]=rho3;
	//
	rho4sq=(1.0-(pow((meta+meta),2)/s));
	if(rho4sq>=0) {
		rho4=std::complex<double>(sqrt(rho4sq),0);
	} else {
		rho4=std::complex<double>(0,sqrt(-rho4sq));
	}
	rho[3]=rho4;
	//
	rho5sq=(1.0-(pow((meta+metap),2)/s));
	if(rho5sq >=0) {
		rho5=std::complex<double>(sqrt(rho5sq),0);
	} else {
		rho5=std::complex<double>(0,sqrt(-rho5sq));
	}
	rho[4]=rho5;

	//sum the pole
	//equation (3) in the E791 K-matrix paper
	for(size_t k = 0; k < 5 ; ++k) {
		for(size_t l = 0; l < 5; ++l) {
			for (size_t pole_index = 0; pole_index < 5; ++pole_index) {
				double A = g[pole_index][k]*g[pole_index][l];
				double B = ma[pole_index]*ma[pole_index]-s;
				K[k][l]  = K[k][l]+std::complex<double>(A/B,0);
			}
		}
	}

	//add the SVT part
	for(size_t k = 0; k < 5; ++k) {
		for(size_t l = 0; l < 5; ++l) {
			double C=f[k][l]*(1.0-s_scatt);
			double D=(s-s_scatt);
			K[k][l]=K[k][l]+std::complex<double>(C/D,0);
		}
	}

	//Include the Alder zero term:
	for(size_t k = 0; k < 5; ++k) {
		for(size_t l = 0; l < 5; ++l) {
			double E=(s-(sa*mpi*mpi*0.5))*(1.0-sa_0);
			double F=(s-sa_0);
			K[k][l]=K[k][l]*std::complex<double>(E/F,0);
		}
	}

	n11=std::complex<double>(1,0)-std::complex<double>(0,1)*K[0][0]*rho[0];
	n12=std::complex<double>(0,0)-std::complex<double>(0,1)*K[0][1]*rho[1];
	n13=std::complex<double>(0,0)-std::complex<double>(0,1)*K[0][2]*rho[2];
	n14=std::complex<double>(0,0)-std::complex<double>(0,1)*K[0][3]*rho[3];
	n15=std::complex<double>(0,0)-std::complex<double>(0,1)*K[0][4]*rho[4];

	n21=std::complex<double>(0,0)-std::complex<double>(0,1)*K[1][0]*rho[0];
	n22=std::complex<double>(1,0)-std::complex<double>(0,1)*K[1][1]*rho[1];
	n23=std::complex<double>(0,0)-std::complex<double>(0,1)*K[1][2]*rho[2];
	n24=std::complex<double>(0,0)-std::complex<double>(0,1)*K[1][3]*rho[3];
	n25=std::complex<double>(0,0)-std::complex<double>(0,1)*K[1][4]*rho[4];

	n31=std::complex<double>(0,0)-std::complex<double>(0,1)*K[2][0]*rho[0];
	n32=std::complex<double>(0,0)-std::complex<double>(0,1)*K[2][1]*rho[1];
	n33=std::complex<double>(1,0)-std::complex<double>(0,1)*K[2][2]*rho[2];
	n34=std::complex<double>(0,0)-std::complex<double>(0,1)*K[2][3]*rho[3];
	n35=std::complex<double>(0,0)-std::complex<double>(0,1)*K[2][4]*rho[4];

	n41=std::complex<double>(0,0)-std::complex<double>(0,1)*K[3][0]*rho[0];
	n42=std::complex<double>(0,0)-std::complex<double>(0,1)*K[3][1]*rho[1];
	n43=std::complex<double>(0,0)-std::complex<double>(0,1)*K[3][2]*rho[2];
	n44=std::complex<double>(1,0)-std::complex<double>(0,1)*K[3][3]*rho[3];
	n45=std::complex<double>(0,0)-std::complex<double>(0,1)*K[3][4]*rho[4];

	n51=std::complex<double>(0,0)-std::complex<double>(0,1)*K[4][0]*rho[0];
	n52=std::complex<double>(0,0)-std::complex<double>(0,1)*K[4][1]*rho[1];
	n53=std::complex<double>(0,0)-std::complex<double>(0,1)*K[4][2]*rho[2];
	n54=std::complex<double>(0,0)-std::complex<double>(0,1)*K[4][3]*rho[3];
	n55=std::complex<double>(1,0)-std::complex<double>(0,1)*K[4][4]*rho[4];

	  //Compute determinant
	det = (n15*n24*n33*n42*n51 - n14*n25*n33*n42*n51 - n15*n23*n34*n42*n51 +
		 n13*n25*n34*n42*n51 + n14*n23*n35*n42*n51 - n13*n24*n35*n42*n51 -
		 n15*n24*n32*n43*n51 + n14*n25*n32*n43*n51 + n15*n22*n34*n43*n51 -
		 n12*n25*n34*n43*n51 - n14*n22*n35*n43*n51 + n12*n24*n35*n43*n51 +
		 n15*n23*n32*n44*n51 - n13*n25*n32*n44*n51 - n15*n22*n33*n44*n51 +
		 n12*n25*n33*n44*n51 + n13*n22*n35*n44*n51 - n12*n23*n35*n44*n51 -
		 n14*n23*n32*n45*n51 + n13*n24*n32*n45*n51 + n14*n22*n33*n45*n51 -
		 n12*n24*n33*n45*n51 - n13*n22*n34*n45*n51 + n12*n23*n34*n45*n51 -
		 n15*n24*n33*n41*n52 + n14*n25*n33*n41*n52 + n15*n23*n34*n41*n52 -
		 n13*n25*n34*n41*n52 - n14*n23*n35*n41*n52 + n13*n24*n35*n41*n52 +
		 n15*n24*n31*n43*n52 - n14*n25*n31*n43*n52 - n15*n21*n34*n43*n52 +
		 n11*n25*n34*n43*n52 + n14*n21*n35*n43*n52 - n11*n24*n35*n43*n52 -
		 n15*n23*n31*n44*n52 + n13*n25*n31*n44*n52 + n15*n21*n33*n44*n52 -
		 n11*n25*n33*n44*n52 - n13*n21*n35*n44*n52 + n11*n23*n35*n44*n52 +
		 n14*n23*n31*n45*n52 - n13*n24*n31*n45*n52 - n14*n21*n33*n45*n52 +
		 n11*n24*n33*n45*n52 + n13*n21*n34*n45*n52 - n11*n23*n34*n45*n52 +
		 n15*n24*n32*n41*n53 - n14*n25*n32*n41*n53 - n15*n22*n34*n41*n53 +
		 n12*n25*n34*n41*n53 + n14*n22*n35*n41*n53 - n12*n24*n35*n41*n53 -
		 n15*n24*n31*n42*n53 + n14*n25*n31*n42*n53 + n15*n21*n34*n42*n53 -
		 n11*n25*n34*n42*n53 - n14*n21*n35*n42*n53 + n11*n24*n35*n42*n53 +
		 n15*n22*n31*n44*n53 - n12*n25*n31*n44*n53 - n15*n21*n32*n44*n53 +
		 n11*n25*n32*n44*n53 + n12*n21*n35*n44*n53 - n11*n22*n35*n44*n53 -
		 n14*n22*n31*n45*n53 + n12*n24*n31*n45*n53 + n14*n21*n32*n45*n53 -
		 n11*n24*n32*n45*n53 - n12*n21*n34*n45*n53 + n11*n22*n34*n45*n53 -
		 n15*n23*n32*n41*n54 + n13*n25*n32*n41*n54 + n15*n22*n33*n41*n54 -
		 n12*n25*n33*n41*n54 - n13*n22*n35*n41*n54 + n12*n23*n35*n41*n54 +
		 n15*n23*n31*n42*n54 - n13*n25*n31*n42*n54 - n15*n21*n33*n42*n54 +
		 n11*n25*n33*n42*n54 + n13*n21*n35*n42*n54 - n11*n23*n35*n42*n54 -
		 n15*n22*n31*n43*n54 + n12*n25*n31*n43*n54 + n15*n21*n32*n43*n54 -
		 n11*n25*n32*n43*n54 - n12*n21*n35*n43*n54 + n11*n22*n35*n43*n54 +
		 n13*n22*n31*n45*n54 - n12*n23*n31*n45*n54 - n13*n21*n32*n45*n54 +
		 n11*n23*n32*n45*n54 + n12*n21*n33*n45*n54 - n11*n22*n33*n45*n54 +
		 n14*n23*n32*n41*n55 - n13*n24*n32*n41*n55 - n14*n22*n33*n41*n55 +
		 n12*n24*n33*n41*n55 + n13*n22*n34*n41*n55 - n12*n23*n34*n41*n55 -
		 n14*n23*n31*n42*n55 + n13*n24*n31*n42*n55 + n14*n21*n33*n42*n55 -
		 n11*n24*n33*n42*n55 - n13*n21*n34*n42*n55 + n11*n23*n34*n42*n55 +
		 n14*n22*n31*n43*n55 - n12*n24*n31*n43*n55 - n14*n21*n32*n43*n55 +
		 n11*n24*n32*n43*n55 + n12*n21*n34*n43*n55 - n11*n22*n34*n43*n55 -
		 n13*n22*n31*n44*n55 + n12*n23*n31*n44*n55 + n13*n21*n32*n44*n55 -
		 n11*n23*n32*n44*n55 - n12*n21*n33*n44*n55 + n11*n22*n33*n44*n55);

	  //The 1st row of the inverse matrix. This matrix is {(I-iKp)^-1}_0j
	i[0][0]   = (n25*n34*n43*n52 -
		     n24*n35*n43*n52 - n25*n33*n44*n52 + n23*n35*n44*n52 +
		     n24*n33*n45*n52 - n23*n34*n45*n52 - n25*n34*n42*n53 +
		     n24*n35*n42*n53 + n25*n32*n44*n53 - n22*n35*n44*n53 -
		     n24*n32*n45*n53 + n22*n34*n45*n53 + n25*n33*n42*n54 -
		     n23*n35*n42*n54 - n25*n32*n43*n54 + n22*n35*n43*n54 +
		     n23*n32*n45*n54 - n22*n33*n45*n54 - n24*n33*n42*n55 +
		     n23*n34*n42*n55 + n24*n32*n43*n55 - n22*n34*n43*n55 -
		     n23*n32*n44*n55 + n22*n33*n44*n55)/det;

	i[0][1]   = (-n15*n34*n43*n52 +
		     n14*n35*n43*n52 + n15*n33*n44*n52 - n13*n35*n44*n52 -
		     n14*n33*n45*n52 + n13*n34*n45*n52 + n15*n34*n42*n53 -
		     n14*n35*n42*n53 - n15*n32*n44*n53 + n12*n35*n44*n53 +
		     n14*n32*n45*n53 - n12*n34*n45*n53 - n15*n33*n42*n54 +
		     n13*n35*n42*n54 + n15*n32*n43*n54 - n12*n35*n43*n54 -
		     n13*n32*n45*n54 + n12*n33*n45*n54 + n14*n33*n42*n55 -
		     n13*n34*n42*n55 - n14*n32*n43*n55 + n12*n34*n43*n55 +
		     n13*n32*n44*n55 - n12*n33*n44*n55)/det;

	i[0][2]   = (n15*n24*n43*n52 -
		     n14*n25*n43*n52 - n15*n23*n44*n52 + n13*n25*n44*n52 +
		     n14*n23*n45*n52 - n13*n24*n45*n52 - n15*n24*n42*n53 +
		     n14*n25*n42*n53 + n15*n22*n44*n53 - n12*n25*n44*n53 -
		     n14*n22*n45*n53 + n12*n24*n45*n53 + n15*n23*n42*n54 -
		     n13*n25*n42*n54 - n15*n22*n43*n54 + n12*n25*n43*n54 +
		     n13*n22*n45*n54 - n12*n23*n45*n54 - n14*n23*n42*n55 +
		     n13*n24*n42*n55 + n14*n22*n43*n55 - n12*n24*n43*n55 -
		     n13*n22*n44*n55 + n12*n23*n44*n55)/det;

	i[0][3]   = (-n15*n24*n33*n52 +
		     n14*n25*n33*n52 + n15*n23*n34*n52 - n13*n25*n34*n52 -
		     n14*n23*n35*n52 + n13*n24*n35*n52 + n15*n24*n32*n53 -
		     n14*n25*n32*n53 - n15*n22*n34*n53 + n12*n25*n34*n53 +
		     n14*n22*n35*n53 - n12*n24*n35*n53 - n15*n23*n32*n54 +
		     n13*n25*n32*n54 + n15*n22*n33*n54 - n12*n25*n33*n54 -
		     n13*n22*n35*n54 + n12*n23*n35*n54 + n14*n23*n32*n55 -
		     n13*n24*n32*n55 - n14*n22*n33*n55 + n12*n24*n33*n55 +
		     n13*n22*n34*n55 - n12*n23*n34*n55)/det;

	i[0][4]   = (n15*n24*n33*n42 -
		     n14*n25*n33*n42 - n15*n23*n34*n42 + n13*n25*n34*n42 +
		     n14*n23*n35*n42 - n13*n24*n35*n42 - n15*n24*n32*n43 +
		     n14*n25*n32*n43 + n15*n22*n34*n43 - n12*n25*n34*n43 -
		     n14*n22*n35*n43 + n12*n24*n35*n43 + n15*n23*n32*n44 -
		     n13*n25*n32*n44 - n15*n22*n33*n44 + n12*n25*n33*n44 +
		     n13*n22*n35*n44 - n12*n23*n35*n44 - n14*n23*n32*n45 +
		     n13*n24*n32*n45 + n14*n22*n33*n45 - n12*n24*n33*n45 -
		     n13*n22*n34*n45 + n12*n23*n34*n45)/det;

	// convert betas and f_prods to TComplex

	std::vector<std::complex<double> > _beta(5, std::complex<double>(0.,0.));

	_beta[0] = _parameters[0] * std::complex<double>(cos(_parameters[1]),sin(_parameters[1]));
	_beta[1] = _parameters[2] * std::complex<double>(cos(_parameters[3]),sin(_parameters[3]));
	_beta[2] = _parameters[4] * std::complex<double>(cos(_parameters[5]),sin(_parameters[5]));
	_beta[3] = _parameters[6] * std::complex<double>(cos(_parameters[7]),sin(_parameters[7]));
	_beta[4] = _parameters[8] * std::complex<double>(cos(_parameters[9]),sin(_parameters[9]));

	std::complex<double> _fr11prod = _parameters[10] * std::complex<double>(cos(_parameters[11]),sin(_parameters[11]));
	std::complex<double> _fr12prod = _parameters[12] * std::complex<double>(cos(_parameters[13]),sin(_parameters[13]));
	std::complex<double> _fr13prod = _parameters[14] * std::complex<double>(cos(_parameters[15]),sin(_parameters[15]));
	std::complex<double> _fr14prod = _parameters[16] * std::complex<double>(cos(_parameters[17]),sin(_parameters[17]));
	std::complex<double> _fr15prod = _parameters[18] * std::complex<double>(cos(_parameters[19]),sin(_parameters[19]));

	double _s0prod = _parameters[20];

	std::vector<std::complex<double> > U1j(5, std::complex<double>(0.,0.));
	for(size_t j = 0; j < 5; ++j) {
		U1j[j] = i[0][j];
	}

	// compute product of inverse matrix times production vector, split production vector into two pieces
	std::complex<double> value0(0.,0.);
	std::complex<double> value1(0.,0.);

	  // compute inverse_matrix times first part of production vector, sum all the poles
	for(size_t l = 0; l < 5; ++l) {
		for (size_t pole_index = 0; pole_index < 5; ++pole_index) {
			std::complex<double> A = _beta[pole_index]*g[pole_index][l];
			double B               = ma[pole_index]*ma[pole_index]-s;
			value0                += U1j[l] * A / B;
		}
	}

	// compute inverse_matrix times second part of production vector
	value1 += U1j[0]*_fr11prod;
	value1 += U1j[1]*_fr12prod;
	value1 += U1j[2]*_fr13prod;
	value1 += U1j[3]*_fr14prod;
	value1 += U1j[4]*_fr15prod;
	value1 *= (1-_s0prod)/(s-_s0prod);

	// Computes the final F0 vector by adding the two parts
	std::complex<double> value = value0 + value1;

//	TComplex F_0_TComplex (value._rpart, value._ipart);

	return value;

}

