#include "Wilson_THDM.h"



void THDM_LO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

    double C7Heps_0,C8Heps_0,C7Heps2_0,C8Heps2_0;
	double lu,ld;

    lu=(*sm)("YUKAWA_CH_U", 33);
	ld=(*sm)("YUKAWA_CH_D", 33);

    double mass_top_muW=(*sm).run.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale,  (*sm)("MASS",6),(*sm)("MASS",5)); //mass top at top ?
	double mass_b_muW=(*sm).run.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale,  (*sm)("MASS",6), (*sm)("MASS",5)); //mass bottom 6 (at pole)

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

    double mass_top_muW=(*sm).run.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale,  (*sm)("MASS",6),(*sm)("MASS",5)); //mass top at top ?
	double mass_b_muW=(*sm).run.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale,  (*sm)("MASS",6), (*sm)("MASS",5)); //mass bottom 6 (at pole)

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

    double mass_top_muW=(*sm).run.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale,  (*sm)("MASS",6),(*sm)("MASS",5)); //mass top at top ?
	double mass_b_muW=(*sm).run.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale,  (*sm)("MASS",6), (*sm)("MASS",5)); //mass bottom 6 (at pole)

    double sw2=pow(sin(atan((*sm)("Coupling",1)/(*sm)("Coupling",2))),2.); //1 = param-> gp and 2 = param->g2

    double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

    double C4H_1=EH(yt,lu);

    double C3H_2=G3H(yt,lu)+Delta3H(yt,lu)*log(pow(scale/(*sm)("MASS",25),2.));
	double C4H_2=G4H(yt,lu)+Delta4H(yt,lu)*log(pow(scale/(*sm)("MASS",25),2.));
	double C5H_2=-C3H_2/10.+2./15.*C4H_1;
	double C6H_2=-3./16.*C3H_2+1./4.*C4H_1;
}