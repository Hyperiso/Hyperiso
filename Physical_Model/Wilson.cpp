#include "Wilson.h"

void WilsonInitializer::init() {
    
}


WilsonManager* WilsonManager::instance = nullptr;

void SM_LO_Strategy::init(Parameters* sm, double scale, WilsonSet& C, QCDParameters& run, Wilson_parameters& W_param) {

    double mass_top_muW=run.runningAlphasCalculation((*sm)("MASS",6)); //mass top at top ?
	double mass_b_muW=run.runningAlphasCalculation((*sm)("MASS",5)); //mass bottom 6 (at pole)
    // or the opposite

    double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
 	double sw2=pow(sin(atan((*sm)("COUPLING",1)/(*sm)("COUPLING",2))),2.); //1 = param-> gp and 2 = param->g2

	//printf("sw2=%.5e\n",sw2);

	double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

	/* LO */
	
	double C2SM_0 = 1.;
	double C7SM_0 = -0.5*A0t(xt)-23./36.;
	double C8SM_0 = -0.5*F0t(xt)-1./3.;
	double C9SM_0 = (1.-4.*sw2)/sw2*C0t(xt)-B0t(xt)/sw2-D0t(xt) +1./4./sw2+38./27.-4./9.*L;
	double C10SM_0 = (B0t(xt)-C0t(xt))/sw2-1./4./sw2;

	if (C.size() < 1) C.resize(1);  
    auto& C_LO = C[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_LO[static_cast<size_t>(WilsonCoefficient::C2)] = complex_t(C2SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C7)] = complex_t(C7SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C8)] = complex_t(C8SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C9)] = complex_t(C9SM_0, 0);
    C_LO[static_cast<size_t>(WilsonCoefficient::C10)] = complex_t(C10SM_0, 0);



}

void SM_NLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C, QCDParameters& run, Wilson_parameters& W_param) {

	double mass_top_muW=run.runningAlphasCalculation((*sm)("MASS",6)); //mass top at top ?
	double mass_b_muW=run.runningAlphasCalculation((*sm)("MASS",5)); //mass bottom 6 (at pole)
    // or the opposite

	double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

	double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
 	double sw2=pow(sin(atan((*sm)("COUPLING",1)/(*sm)("COUPLING",2))),2.); //1 = param-> gp and 2 = param->g2


    double C1SM_1 = 15.+6.*L;
	double C4SM_1 = E0t(xt)-7./9.+2./3.*L;
	double C7SM_1 = -0.5*A1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))+713./243.+4./81.*L-4./9.*C4SM_1;
	double C8SM_1 = -0.5*F1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))+91./324.-4./27.*L-C4SM_1/6.;
	double C9SM_1 = (1.-4.*sw2)/sw2*C1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))-B1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))/sw2-D1t(xt,log(scale*scale/mass_top_muW/mass_top_muW)) +1./sw2+524./729.-128./243.*pi*pi-16./3.*L-128./81.*L*L;
	double C10SM_1 = (B1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))-C1t(xt,log(scale*scale/mass_top_muW/mass_top_muW)))/sw2-1./sw2;

	if (C.size() < 1) C.resize(1);  
    auto& C_NLO = C[1];
	C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_NLO[static_cast<size_t>(WilsonCoefficient::C1)] = complex_t(C1SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C4)] = complex_t(C4SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C7)] = complex_t(C7SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C8)] = complex_t(C8SM_1, 0);
    C_NLO[static_cast<size_t>(WilsonCoefficient::C9)] = complex_t(C9SM_1, 0);
	C_NLO[static_cast<size_t>(WilsonCoefficient::C10)] = complex_t(C10SM_1, 0);
}


void SM_NNLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C, QCDParameters& run, Wilson_parameters& W_param) {

	double mass_top_muW=run.runningAlphasCalculation((*sm)("MASS",6)); //mass top at top ?
	double mass_b_muW=run.runningAlphasCalculation((*sm)("MASS",5)); //mass bottom 6 (at pole)
    // or the opposite

	double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

	double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
 	double sw2=pow(sin(atan((*sm)("COUPLING",1)/(*sm)("COUPLING",2))),2.); //1 = param-> gp and 2 = param->g2

    double C1SM_2 = -T(xt)+7987./72.+17.*pi*pi/3.+475./6.*L+17.*L*L;
	double C2SM_2 = 127./18.+4./3.*pi*pi+46./3.*L+4.*L*L;
	double C3SM_2 = G1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))-680./243.-20./81.*pi*pi-68./81.*L-20./27.*L*L;
	double C4SM_2 = E1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))+950./243.+10./81.*pi*pi+124./27.*L+10./27.*L*L;
	double C5SM_2 = -G1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))/10.+2./15.*E0t(xt)+68./243.+2./81.*pi*pi+14./81.*L+2./27.*L*L;
	double C6SM_2 = -3./16.*G1t(xt,log(scale*scale/mass_top_muW/mass_top_muW))+E0t(xt)/4.+85./162.+5./108.*pi*pi+35./108.*L+5./36.*L*L;

	double xtW=pow(run.runningAlphasCalculation((*sm)("MASS",6))/(*sm)("MASS",24),2.); // mass top at mass top
	double xtt=pow((*sm)("MASS",6)/(*sm)("MASS",24),2.); // 24 -> W
	
	double C7SM_2 = (C7t2mt(xtt)+log(scale*scale/mass_top_muW/mass_top_muW)*((-592.*pow(xt,5.)-22.*pow(xt,4.)+12814.*pow(xt,3.)-6376.*xt*xt+512.*xt)/27./pow(xt-1.,5.)*Li2(1.-1./xt)
	+(-26838.*pow(xt,5.)+25938.*pow(xt,4.)+627367.*pow(xt,3.)-331956.*xt*xt+16989.*xt-460.)/729./pow(xt-1.,6.)*log(xt)
	+(34400.*pow(xt,5.)+276644.*pow(xt,4.)-2668324.*pow(xt,3.)+1694437.*xt*xt-323354.*xt+53077.)/2187./pow(xt-1.,5.)
	+log(scale*scale/mass_top_muW/mass_top_muW)*((-63.*pow(xt,5.)+532.*pow(xt,4.)+2089.*pow(xt,3.)-1118.*xt*xt)/9./pow(xt-1.,6.)*log(xt)
	+(1186.*pow(xt,5.)-2705.*pow(xt,4.)-24791.*pow(xt,3.)-16099.*xt*xt+19229.*xt-2740.)/162./pow(xt-1.,5.))) )
	-(C7c2MW(xtW)+13763./2187.*log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24))+814./729.*pow(log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)),2.));

	double C8SM_2 = (C8t2mt(xtt)+log(scale*scale/mass_top_muW/mass_top_muW)*((-148.*pow(xt,5.)+1052.*pow(xt,4.)-4811.*pow(xt,3.)-3520.*xt*xt-61.*xt)/18./pow(xt-1.,5.)*Li2(1.-1./xt)
	+(-15984.*pow(xt,5.)+152379.*pow(xt,4.)-1358060.*pow(xt,3.)-1201653.*xt*xt-74190.*xt+9188.)/1944./pow(xt-1.,6.)*log(xt)
	+(109669.*pow(xt,5.)-1112675.*pow(xt,4.)+6239377.*pow(xt,3.)+8967623.*xt*xt+768722.*xt-42796.)/11664./pow(xt-1.,5.)
	+log(scale*scale/mass_top_muW/mass_top_muW)*((-139.*pow(xt,4.)-2938.*pow(xt,3.)-2683.*xt*xt)/12./pow(xt-1.,6.)*log(xt)
	+(1295.*pow(xt,5.)-7009.*pow(xt,4.)+29495.*pow(xt,3.)+64513.*xt*xt+17458.*xt-2072.)/216./pow(xt-1.,5.))) )
	-(C8c2MW(xtW)+16607./5832.*log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24))+397./486.*pow(log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)),2.));
	
	if (C.size() < 1) C.resize(1);  
    auto& C_NNLO = C[2];
	C_NNLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));

	C_NNLO[static_cast<size_t>(WilsonCoefficient::C1)] = complex_t(C1SM_2, 0);
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C2)] = complex_t(C2SM_2, 0);
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C3)] = complex_t(C3SM_2, 0);
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C4)] = complex_t(C4SM_2, 0);
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C5)] = complex_t(C5SM_2, 0);
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C6)] = complex_t(C6SM_2, 0);
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C7)] = complex_t(C7SM_2, 0);
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C8)] = complex_t(C8SM_2, 0);


}
