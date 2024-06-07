#include "Wilson_THDM.h"



void THDM_LO_Strategy::init(double scale, WilsonSet& C_match) {

	Parameters* sm = Parameters::GetInstance();
	Parameters* mod = nullptr;

	if (is_thdm) {
		mod = Parameters::GetInstance(2);
		lu=(*mod)("YU", 22);
		ld=(*mod)("YD", 22);
		Logger::getInstance()->info("YU : " + std::to_string(lu));
		Logger::getInstance()->info("37 : " + std::to_string((*mod)("MASS",37)));
	}
	else {
		mod = Parameters::GetInstance(1);
	}
    double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale, "running", "pole"); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale); //mass bottom 6 (at pole)

    double sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.); //1 = param-> gp and 2 = param->g2

	Logger::getInstance()->info("g1 : " + std::to_string((*sm)("GAUGE",1)));
	Logger::getInstance()->info("g2 : " + std::to_string((*sm)("GAUGE",2)));
	Logger::getInstance()->info("mt : " + std::to_string((*sm)("MASS",6)));
    double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*mod)("MASS",37),2.); // param->mass_H (25)
	Logger::getInstance()->info("xt : " + std::to_string(xt));
    complex_t C7H_0=1./3.*lu*lu*F7_1(yt) - lu*ld*F7_2(yt);
	complex_t C8H_0=1./3.*lu*lu*F8_1(yt) - lu*ld*F8_2(yt);

	complex_t C9H_0=(1.-4.*sw2)/sw2*C9llH0(xt,yt,lu)-D9H0(yt,lu);
	complex_t C10H_0=-C9llH0(xt,yt,lu)/sw2;

	if (C_match.empty()) C_match.resize(1);
	auto& C_LO = C_match[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));

	C_LO[static_cast<size_t>(WilsonCoefficient::C7)] = C7H_0;
	C_LO[static_cast<size_t>(WilsonCoefficient::C8)] = C8H_0;
	C_LO[static_cast<size_t>(WilsonCoefficient::C9)] = C9H_0;
	C_LO[static_cast<size_t>(WilsonCoefficient::C10)] = C10H_0;

	Logger::getInstance()->info("THDM LO Wilson Coefficient Initialized at scale " +std::to_string(scale)+" terminated successfully");

}

void THDM_NLO_Strategy::init(double scale, WilsonSet& C_match) {

	Parameters* sm = Parameters::GetInstance();
	// Parameters* susy = Parameters::GetInstance(1);
	Parameters* mod = nullptr;

    if (lu == -1 || ld == -1) {
		mod = Parameters::GetInstance(2);
		lu=(*mod)("YU", 22);
		ld=(*mod)("YD", 22);
	}
	else {
		mod = Parameters::GetInstance(1);
	}

    double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale); //mass bottom 6 (at pole)

    double sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.); //1 = param-> gp and 2 = param->g2
	double m_H = (*mod)("MASS", 37); // Charged Higgs mass (37)
    double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/m_H,2.); // param->mass_H (25)
    complex_t C4H_1=EH(yt,lu);

    complex_t C7H_1= G7H(yt,lu,ld)+Delta7H(yt,lu,ld)*log(pow(scale/m_H,2.))-4./9.*C4H_1;
	complex_t C8H_1= G8H(yt,lu,ld)+Delta8H(yt,lu,ld)*log(pow(scale/m_H,2.))-1./6.*C4H_1;
	complex_t C9H_1=(1.-4.*sw2)/sw2*C9llH1(xt,yt,lu,log(pow(scale/m_H,2.)))-D9H1(yt,lu,log(pow(scale/m_H,2.)));
	complex_t C10H_1=-C9llH1(xt,yt,lu,log(pow(scale/m_H,2.)))/sw2;
	Logger::getInstance()->info("C10H_1 : " + std::to_string(std::real(C10H_1)));
	if (C_match.size() < 2) C_match.resize(2);
	auto& C_NLO = C_match[1];
	C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));

	double alphas_mu = sm->QCDRunner.runningAlphasCalculation(scale);

	THDM_LO_Strategy::init(scale, C_match);
	auto& C_LO = C_match[0];

	auto adjustCoefficient = [&](std::complex<double>& Cx_NLO, int index) {
        double ratio = alphas_mu / (4.0 * PI);
        double absCx_NLO = std::abs(Cx_NLO) * ratio;
        if (absCx_NLO > std::abs(C_LO[index])) {
            Cx_NLO *= std::abs(C_LO[index]) / std::abs(Cx_NLO) * (1.0 / ratio);
        }
    };

	Logger::getInstance()->info("C8H_1 in THDM " + std::to_string(std::real(C8H_1)));
	// adjustCoefficient(C7H_1, 7);
	// adjustCoefficient(C8H_1, 8);
	// adjustCoefficient(C9H_1, 9);
	// adjustCoefficient(C10H_1, 10);
	Logger::getInstance()->info("C8H_1 in THDM " + std::to_string(std::real(C8H_1)));
	Logger::getInstance()->info("C7H_1 in THDM " + std::to_string(std::real(C7H_1)));
	Logger::getInstance()->info("C4H_1 in THDM " + std::to_string(std::real(C4H_1)));
	// Logger::getInstance()->info("C4Char_1 in THDM " + std::to_string(std::real(C_LO[1])));
	// Logger::getInstance()->info("C4Char_1 in THDM " + std::to_string(std::real(C_NLO[5])));
	C_NLO[static_cast<size_t>(WilsonCoefficient::C4)] = C4H_1;
	C_NLO[static_cast<size_t>(WilsonCoefficient::C7)] = C7H_1;
	C_NLO[static_cast<size_t>(WilsonCoefficient::C8)] = C8H_1;
	C_NLO[static_cast<size_t>(WilsonCoefficient::C9)] = C9H_1;
	C_NLO[static_cast<size_t>(WilsonCoefficient::C10)] = C10H_1;

	Logger::getInstance()->info("THDM NLO Wilson Coefficient Initialized at scale " +std::to_string(scale)+" terminated successfully");
}

void THDM_NNLO_Strategy::init(double scale, WilsonSet& C_match) {

	Parameters* sm = Parameters::GetInstance();
	// Parameters* susy = Parameters::GetInstance(1);
	// Parameters* thdm = Parameters::GetInstance(2);
	Parameters* mod = nullptr;

    if (lu == -1 || ld == -1) {
		mod = Parameters::GetInstance(2);
		lu=(*mod)("YU", 22);
		ld=(*mod)("YD", 22);
	}
	else {
		mod = Parameters::GetInstance(1);
	}
    double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale); //mass bottom 6 (at pole)

    double sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.); //1 = param-> gp and 2 = param->g2

    double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*mod)("MASS",37),2.); // param->mass_H (25)

    complex_t C4H_1=EH(yt,lu);

    complex_t C3H_2=G3H(yt,lu)+Delta3H(yt,lu)*log(pow(scale/(*mod)("MASS",37),2.));
	complex_t C4H_2=G4H(yt,lu)+Delta4H(yt,lu)*log(pow(scale/(*mod)("MASS",37),2.));
	complex_t C5H_2=-C3H_2/10.+2./15.*C4H_1;
	complex_t C6H_2=-3./16.*C3H_2+1./4.*C4H_1;

	complex_t C7H_2 =C7H2(yt,lu,ld,log(pow(scale/mass_top_muW, 2.)));
	complex_t C8H_2 =C8H2(yt,lu,ld,log(pow(scale/mass_top_muW, 2.)));

	


	if (C_match.empty()) C_match.resize(3);
	auto& C_NNLO = C_match[2];
	C_NNLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));

	C_NNLO[static_cast<size_t>(WilsonCoefficient::C3)] = C3H_2;
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C4)] = C4H_2;
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C5)] = C5H_2;
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C6)] = C6H_2;
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C7)] = C7H_2;
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C8)] = C8H_2;

	Logger::getInstance()->info("THDM NNLO Wilson Coefficient Initialized at scale " +std::to_string(scale)+" terminated successfully");
}

void THDM_LO_Strategy::init_scalar(double Q_match,double Q,int gen, WilsonSet& C) {
	/* Wilson coefficients CQ1 et CQ2 in 2HDM */ 
	
	Parameters* sm = Parameters::GetInstance(0);
	Parameters* thdm = Parameters::GetInstance(2);
    double ml;

	
	if(gen==1) ml=(*sm)("MASS", 11);
	else if(gen==3) ml=(*sm)("MASS", 15);
	else {gen=2; ml=(*sm)("MASS", 13);}
	
	double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),Q_match); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), Q_match); //mass bottom 6 (at pole)
	double sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.);

	double xt=pow(mass_top_muW/(*sm)("MASS",24),2.);

	double xh=pow((*thdm)("MASS",25)/(*sm)("MASS",24),2.);
	Logger::getInstance()->info("xh : " + std::to_string(xh));
	int nf=5;
	double beta0 = 11.-2./3.*nf;

	double alphas_muW=(*sm).QCDRunner.runningAlphasCalculation(Q_match);
	double alphas_mu=(*sm).QCDRunner.runningAlphasCalculation(Q);	
	double eta_mu=alphas_muW/alphas_mu;

	double alpha=(*thdm)("ALPHA", 0);
	Logger::getInstance()->info("alPHA : "+ std::to_string(alpha));
	double beta=atan((*thdm)("HMIX", 2));
	Logger::getInstance()->info("beta : "+ std::to_string(beta));
	double xH=pow((*thdm)("MASS",37)/(*sm)("MASS",24),2.);
	double xH0=pow((*thdm)("MASS",35)/(*sm)("MASS",24),2.);
	double xA=pow((*thdm)("MASS",36)/(*sm)("MASS",24),2.);
	
	Logger::getInstance()->info("xt : "+ std::to_string(xt));
	Logger::getInstance()->info("xH : "+ std::to_string(xH));
	Logger::getInstance()->info("xH0 : "+ std::to_string(xH0));
	Logger::getInstance()->info("mass_A : "+ std::to_string((*thdm)("MASS",36)));
	Logger::getInstance()->info("mass_W : "+ std::to_string((*sm)("MASS",24)));

	Logger::getInstance()->info("YD : "+ std::to_string(((*thdm)("YD",22))));
	Logger::getInstance()->info("YU : "+ std::to_string((*thdm)("YU",22)));
	Logger::getInstance()->info("YL : "+ std::to_string((*thdm)("YL",10*(gen-1)+gen-1)));


	double G1=-3./4.+(*thdm)("YD",22)*(*thdm)("YU",22)*F4SP(xt,xH)+(*thdm)("YU",22)*(*thdm)("YU",22)*F5SP(xt,xH);
	Logger::getInstance()->info("G1 : "+ std::to_string(G1));
	double G2=(*thdm)("YD",22)*((*thdm)("YD",22)*(*thdm)("YU",22)+1.)*F6SP(xt,xH)-(*thdm)("YD",22)*(*thdm)("YU",22)*(*thdm)("YU",22)*F7SP(xt,xH)
	+(*thdm)("YU",22)*(*thdm)("YU",22)*((*thdm)("YD",22)*F8SP(xt,xH)+(*thdm)("YU",22)*F9SP(xt,xH)-(*thdm)("YU",22)*F10SP(xt,xH))+(*thdm)("YU",22)*F11SP(xt,xH)-(*thdm)("YU",22)*F12SP(xt,xH);
	Logger::getInstance()->info("G2 : "+ std::to_string(G2));
	double G3=(*thdm)("YD",22)*((*thdm)("YD",22)*(*thdm)("YU",22)+1.)*F6SP(xt,xH)+(*thdm)("YD",22)*(*thdm)("YU",22)*(*thdm)("YU",22)*F7SP(xt,xH)
	+(*thdm)("YU",22)*(*thdm)("YU",22)*((*thdm)("YD",22)*F8SP(xt,xH)+(*thdm)("YU",22)*F9SP(xt,xH)+(*thdm)("YU",22)*F10SP(xt,xH))+(*thdm)("YU",22)*F11SP(xt,xH)+(*thdm)("YD",22)*F12SP(xt,xH);
	Logger::getInstance()->info("G3 : "+ std::to_string(G3));
	double CSn_2HDM=xt*(F0SP(xt)+(*thdm)("YL",10*(gen-1)+gen-1)*((*thdm)("YD",22)*F1SP(xt,xH)+(*thdm)("YU",22)*F2SP(xt,xH))+(*thdm)("YL",10*(gen-1)+gen-1)*(*thdm)("YU",22)*F3SP(xt,xH))
	+xt/2./xh*(sin(alpha-beta)+cos(alpha-beta)*(*thdm)("YL",10*(gen-1)+gen-1))*(sin(alpha-beta)*G1+cos(alpha-beta)*G2)
	+xt/2./xH0*(cos(alpha-beta)-sin(alpha-beta)*(*thdm)("YL",10*(gen-1)+gen-1))*(cos(alpha-beta)*G1-sin(alpha-beta)*G2);

	Logger::getInstance()->info("1er terme : "+ std::to_string(xt*(F0SP(xt)+(*thdm)("YL",10*(gen-1)+gen-1)*((*thdm)("YD",22)*F1SP(xt,xH)+(*thdm)("YU",22)*F2SP(xt,xH))+(*thdm)("YL",10*(gen-1)+gen-1)*(*thdm)("YU",22)*F3SP(xt,xH))));
	Logger::getInstance()->info("2eme terme : "+ std::to_string(xt/2./xh*(sin(alpha-beta)+cos(alpha-beta)*(*thdm)("YL",10*(gen-1)+gen-1))*(sin(alpha-beta)*G1+cos(alpha-beta)*G2)));
	Logger::getInstance()->info("3eme terme : "+ std::to_string(xt/2./xH0*(cos(alpha-beta)-sin(alpha-beta)*(*thdm)("YL",10*(gen-1)+gen-1))*(cos(alpha-beta)*G1-sin(alpha-beta)*G2)));

	
	Logger::getInstance()->info("F0SP(xt) : "+ std::to_string(F0SP(xt)));
	Logger::getInstance()->info("F1SP(xt,xH) : "+ std::to_string(F1SP(xt,xH)));
	Logger::getInstance()->info("F2SP(xt,xH) : "+ std::to_string(F2SP(xt,xH)));
	Logger::getInstance()->info("F3SP(xt,xH) : "+ std::to_string(F3SP(xt,xH)));
	
	Logger::getInstance()->info("sin : "+ std::to_string(sin(alpha-beta)));
	Logger::getInstance()->info("cos : "+ std::to_string(cos(alpha-beta)));

	Logger::getInstance()->info("CSn_2HDM : "+ std::to_string(CSn_2HDM));
	double CPn_2HDM=xt*(-(*thdm)("YL",10*(gen-1)+gen-1)*((*thdm)("YD",22)*F1SP(xt,xH)+(*thdm)("YU",22)*F2SP(xt,xH))+(*thdm)("YL",10*(gen-1)+gen-1)*(*thdm)("YU",22)*F3SP(xt,xH))+xt/2./xA*((*thdm)("YL",10*(gen-1)+gen-1))*G3;
	Logger::getInstance()->info("CPn_2HDM : "+ std::to_string(CPn_2HDM));
	double CQ1H_0=CSc_2HDM(xH,xt,(*thdm)("YU",22),(*thdm)("YD",22),(*thdm)("YL",10*(gen-1)+gen-1))+CSn_2HDM;
	double CQ2H_0=CPc_2HDM(xH,xt,(*thdm)("YU",22),(*thdm)("YD",22),(*thdm)("YL",10*(gen-1)+gen-1),sw2)+CPn_2HDM;
	
	CQ1H_0*=(ml*mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/sw2;
	CQ2H_0*=(ml*mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/sw2;
	
	if (C.size() < 1) C.resize(1); 
    auto& C_LO = C[0]; 
    C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

    C_LO[static_cast<size_t>(WilsonCoefficient::CQ1)] = CQ1H_0*pow(eta_mu,-4./beta0);
    C_LO[static_cast<size_t>(WilsonCoefficient::CQ2)]= CQ2H_0*pow(eta_mu,-4./beta0);

	Logger::getInstance()->info("SUSY LO Wilson Scalars Coefficient Initialized from scale " +std::to_string(Q_match)+" to scale" + std::to_string(Q) + " terminated successfully");

}