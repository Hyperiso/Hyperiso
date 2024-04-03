#include "Wilson_THDM.h"



void THDM_LO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

    double C7Heps_0,C8Heps_0,C7Heps2_0,C8Heps2_0;
	double lu,ld;

    lu=(*sm)("YUKAWA_CH_U", 33);
	ld=(*sm)("YUKAWA_CH_D", 33);

    double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale); //mass bottom 6 (at pole)

    double sw2=pow(sin(atan((*sm)("Coupling",1)/(*sm)("Coupling",2))),2.); //1 = param-> gp and 2 = param->g2

    double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

    double C7H_0=1./3.*lu*lu*F7_1(yt) - lu*ld*F7_2(yt);
	double C8H_0=1./3.*lu*lu*F8_1(yt) - lu*ld*F8_2(yt);

	double C9H_0=(1.-4.*sw2)/sw2*C9llH0(xt,yt,lu)-D9H0(yt,lu);
	double C10H_0=-C9llH0(xt,yt,lu)/sw2;
}

void THDM_NLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

    double C7Heps_0,C8Heps_0,C7Heps2_0,C8Heps2_0;
	double lu,ld;

    lu=(*sm)("YUKAWA_CH_U", 33);
	ld=(*sm)("YUKAWA_CH_D", 33);

    double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale); //mass bottom 6 (at pole)

    double sw2=pow(sin(atan((*sm)("Coupling",1)/(*sm)("Coupling",2))),2.); //1 = param-> gp and 2 = param->g2

    double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

    double C4H_1=EH(yt,lu);

    double C7H_1= G7H(yt,lu,ld)+Delta7H(yt,lu,ld)*log(pow(scale/(*sm)("MASS",25),2.))-4./9.*C4H_1;
	double C8H_1= G8H(yt,lu,ld)+Delta8H(yt,lu,ld)*log(pow(scale/(*sm)("MASS",25),2.))-1./6.*C4H_1;
	double C9H_1=(1.-4.*sw2)/sw2*C9llH1(xt,yt,lu,log(pow(scale/(*sm)("MASS",25),2.)))-D9H1(yt,lu,log(pow(scale/(*sm)("MASS",25),2.)));
	double C10H_1=-C9llH1(xt,yt,lu,log(pow(scale/(*sm)("MASS",25),2.)))/sw2;
}

void THDM_NNLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

    double C7Heps_0,C8Heps_0,C7Heps2_0,C8Heps2_0;
	double lu,ld;

    lu=(*sm)("YUKAWA_CH_U", 33);
	ld=(*sm)("YUKAWA_CH_D", 33);

    double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale); //mass bottom 6 (at pole)

    double sw2=pow(sin(atan((*sm)("Coupling",1)/(*sm)("Coupling",2))),2.); //1 = param-> gp and 2 = param->g2

    double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

    double C4H_1=EH(yt,lu);

    double C3H_2=G3H(yt,lu)+Delta3H(yt,lu)*log(pow(scale/(*sm)("MASS",25),2.));
	double C4H_2=G4H(yt,lu)+Delta4H(yt,lu)*log(pow(scale/(*sm)("MASS",25),2.));
	double C5H_2=-C3H_2/10.+2./15.*C4H_1;
	double C6H_2=-3./16.*C3H_2+1./4.*C4H_1;
}

void THDM_LO_Strategy::init_scalar(double Q_match,double Q,int gen, WilsonSet& C) {
	/* Wilson coefficients CQ1 et CQ2 in 2HDM */ 
	
	Parameters* sm = Parameters::GetInstance(0);
    double ml;

	
	if(gen==1) ml=(*sm)("MASS", 11);
	else if(gen==3) ml=(*sm)("MASS", 13);
	else {gen=2; ml=(*sm)("MASS", 15);}
	
	double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),Q_match); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), Q_match); //mass bottom 6 (at pole)
	double sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.);

	double xt=pow(mass_top_muW/(*sm)("MASS",24),2.);

	double xh=pow((*sm)("MASS",25)/(*sm)("MASS",24),2.);

	int nf=5;
	double beta0 = 11.-2./3.*nf;

	double alphas_muW=(*sm).QCDRunner.runningAlphasCalculation(Q_match);
	double alphas_mu=(*sm).QCDRunner.runningAlphasCalculation(Q);	
	double eta_mu=alphas_muW/alphas_mu;

	double alpha=(*sm)("ALPHA", 42);
	double beta=atan((*sm)("EXTPAR", 25));

	double xH=pow((*sm)("MASS",37)/(*sm)("MASS",24),2.);
	double xH0=pow((*sm)("MASS",35)/(*sm)("MASS",24),2.);
	double xA=pow((*sm)("MASS",36)/(*sm)("MASS",24),2.);
	
	double G1=-3./4.+(*sm)("YUKAWA_CH_D",33)*(*sm)("YUKAWA_CH_U",33)*F4SP(xt,xH)+(*sm)("YUKAWA_CH_U",33)*(*sm)("YUKAWA_CH_U",33)*F5SP(xt,xH);
	
	double G2=(*sm)("YUKAWA_CH_D",33)*((*sm)("YUKAWA_CH_D",33)*(*sm)("YUKAWA_CH_U",33)+1.)*F6SP(xt,xH)-(*sm)("YUKAWA_CH_D",33)*(*sm)("YUKAWA_CH_U",33)*(*sm)("YUKAWA_CH_U",33)*F7SP(xt,xH)
	+(*sm)("YUKAWA_CH_U",33)*(*sm)("YUKAWA_CH_U",33)*((*sm)("YUKAWA_CH_D",33)*F8SP(xt,xH)+(*sm)("YUKAWA_CH_U",33)*F9SP(xt,xH)-(*sm)("YUKAWA_CH_U",33)*F10SP(xt,xH))+(*sm)("YUKAWA_CH_U",33)*F11SP(xt,xH)-(*sm)("YUKAWA_CH_U",33)*F12SP(xt,xH);
	
	double G3=(*sm)("YUKAWA_CH_D",33)*((*sm)("YUKAWA_CH_D",33)*(*sm)("YUKAWA_CH_U",33)+1.)*F6SP(xt,xH)+(*sm)("YUKAWA_CH_D",33)*(*sm)("YUKAWA_CH_U",33)*(*sm)("YUKAWA_CH_U",33)*F7SP(xt,xH)
	+(*sm)("YUKAWA_CH_U",33)*(*sm)("YUKAWA_CH_U",33)*((*sm)("YUKAWA_CH_D",33)*F8SP(xt,xH)+(*sm)("YUKAWA_CH_U",33)*F9SP(xt,xH)+(*sm)("YUKAWA_CH_U",33)*F10SP(xt,xH))+(*sm)("YUKAWA_CH_U",33)*F11SP(xt,xH)+(*sm)("YUKAWA_CH_D",33)*F12SP(xt,xH);

	double CSn_2HDM=xt*(F0SP(xt)+(*sm)("YUKAWA_CH_L",gen*gen)*((*sm)("YUKAWA_CH_D",33)*F1SP(xt,xH)+(*sm)("YUKAWA_CH_U",33)*F2SP(xt,xH))+(*sm)("YUKAWA_CH_L",gen*gen)*(*sm)("YUKAWA_CH_U",33)*F3SP(xt,xH))
	+xt/2./xh*(sin(alpha-beta)+cos(alpha-beta)*(*sm)("YUKAWA_CH_L",gen*gen))*(sin(alpha-beta)*G1+cos(alpha-beta)*G2)
	+xt/2./xH0*(cos(alpha-beta)-sin(alpha-beta)*(*sm)("YUKAWA_CH_L",gen*gen))*(cos(alpha-beta)*G1-sin(alpha-beta)*G2);
	
	double CPn_2HDM=xt*(-(*sm)("YUKAWA_CH_L",gen*gen)*((*sm)("YUKAWA_CH_D",33)*F1SP(xt,xH)+(*sm)("YUKAWA_CH_U",33)*F2SP(xt,xH))+(*sm)("YUKAWA_CH_L",gen*gen)*(*sm)("YUKAWA_CH_U",33)*F3SP(xt,xH))+xt/2./xA*(*sm)("YUKAWA_CH_L",gen*gen)*G3;

	double CQ1H_0=CSc_2HDM(xH,xt,(*sm)("YUKAWA_CH_U",33),(*sm)("YUKAWA_CH_D",33),(*sm)("YUKAWA_CH_L",gen*gen))+CSn_2HDM;
	double CQ2H_0=CPc_2HDM(xH,xt,(*sm)("YUKAWA_CH_U",33),(*sm)("YUKAWA_CH_D",33),(*sm)("YUKAWA_CH_L",gen*gen),sw2)+CPn_2HDM;
	
	CQ1H_0*=(ml*mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/sw2;
	CQ2H_0*=(ml*mass_b_muW/(*sm)("MASS",24)/(*sm)("MASS",24))/sw2;
	
	if (C.size() < 1) C.resize(1); 
    auto& C_LO = C[0]; 
    C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

    C_LO[static_cast<size_t>(WilsonCoefficient::CQ1)] = CQ1H_0*pow(eta_mu,-4./beta0);
    C_LO[static_cast<size_t>(WilsonCoefficient::CQ2)]= CQ2H_0*pow(eta_mu,-4./beta0);

}