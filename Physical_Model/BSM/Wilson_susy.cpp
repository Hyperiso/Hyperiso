#include "Wilson_susy.h"
#include "susy_parameters.h"



void SUSY_LO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

	EpsilonCalculator* epsi = EpsilonCalculator::GetInstance();

	Parameters* susy = Parameters::GetInstance(1);
	susy_parameters* sus_param = susy_parameters::GetInstance(scale);


	double C7SMeps_0= ((*sus_param).epsilonb-(*sus_param).epsilonbp)/(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25))*(*susy)("EXTPAR",25)*F7_2((*sus_param).xt);
	double C8SMeps_0= ((*sus_param).epsilonb-(*sus_param).epsilonbp)/(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25))*(*susy)("EXTPAR",25)*F8_2((*sus_param).xt);


	double C7Heps_0=(-(*sus_param).epsilon0p-(*sus_param).epsilonb)/(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25))*(*susy)("EXTPAR",25)*F7_2((*sus_param).yt);
	double C8Heps_0=(-(*sus_param).epsilon0p-(*sus_param).epsilonb)/(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25))*(*susy)("EXTPAR",25)*F8_2((*sus_param).yt);

	double C7Heps2_0=0.;
	double C8Heps2_0=0.;

	if(((*sus_param).mass_A02==0.)&&((*sus_param).mass_H03==0.))
	{
		C7Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("EXTPAR",25),2.)/(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25))/(1.+(*sus_param).epsilon0*(*susy)("EXTPAR",25))*F7_2((*sus_param).yt);
		C7Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25),2.)*(1.+pow((*susy)("EXTPAR",25),2.))/(1.+(*sus_param).epsilon0*(*susy)("EXTPAR",25))/72.		*((cos((*susy)("ALPHA",0))+sin((*susy)("ALPHA",0))*(*susy)("EXTPAR",25))*(-sin((*susy)("ALPHA",0))+(*sus_param).epsilonb*cos((*susy)("ALPHA",0)))*pow((*sus_param).mass_b_muW/(*susy)("MASS",25),2.)
		+(sin((*susy)("ALPHA",0))-cos((*susy)("ALPHA",0))*(*susy)("EXTPAR",25))*(cos((*susy)("ALPHA",0))+(*sus_param).epsilonb*sin((*susy)("ALPHA",0)))*pow((*sus_param).mass_b_muW/(*susy)("MASS",35),2.)			+(-cos(atan((*susy)("EXTPAR",25)))-sin(atan((*susy)("EXTPAR",25)))*(*susy)("EXTPAR",25))*(sin(atan((*susy)("EXTPAR",25)))-(*sus_param).epsilonb*cos(atan((*susy)("EXTPAR",25))))*pow((*sus_param).mass_b_muW/(*susy)("MASS",36),2.));

		C8Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("EXTPAR",25),2.)/(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25))/(1.+(*sus_param).epsilon0*(*susy)("EXTPAR",25))*F8_2((*sus_param).yt);
		C8Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25),2.)*(1.+pow((*susy)("EXTPAR",25),2.))/(1.+(*sus_param).epsilon0*(*susy)("EXTPAR",25))/72.		*((cos((*susy)("ALPHA",0))+sin((*susy)("ALPHA",0))*(*susy)("EXTPAR",25))*(-sin((*susy)("ALPHA",0))+(*sus_param).epsilonb*cos((*susy)("ALPHA",0)))*pow((*sus_param).mass_b_muW/(*susy)("MASS",25),2.)
		+(sin((*susy)("ALPHA",0))-cos((*susy)("ALPHA",0))*(*susy)("EXTPAR",25))*(cos((*susy)("ALPHA",0))+(*sus_param).epsilonb*sin((*susy)("ALPHA",0)))*pow((*sus_param).mass_b_muW/(*susy)("MASS",35),2.)			+(-cos(atan((*susy)("EXTPAR",25)))-sin(atan((*susy)("EXTPAR",25)))*(*susy)("EXTPAR",25))*(sin(atan((*susy)("EXTPAR",25)))-(*sus_param).epsilonb*cos(atan((*susy)("EXTPAR",25))))*pow((*sus_param).mass_b_muW/(*susy)("MASS",36),2.));
	}
	else
	{		
		C7Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("EXTPAR",25),2.)/(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25))/(1.+(*sus_param).epsilon0*(*susy)("EXTPAR",25))*F7_2((*sus_param).yt);
		C7Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25),2.)*(1.+pow((*susy)("EXTPAR",25),2.))/(1.+(*sus_param).epsilon0*(*susy)("EXTPAR",25))/72.	*(((*susy)("HMIX", 11)+(*susy)("HMIX", 12)*(*susy)("EXTPAR",25))*(-(*susy)("HMIX", 12)+(*sus_param).epsilonb*(*susy)("HMIX", 11))*pow((*sus_param).mass_b_muW/(*susy)("MASS",25),2.)
		+((*susy)("HMIX", 21)+(*susy)("HMIX", 22)*(*susy)("EXTPAR",25))*(-(*susy)("HMIX", 22)+(*sus_param).epsilonb*(*susy)("HMIX", 21))*pow((*sus_param).mass_b_muW/(*susy)("MASS",35),2.)
		+((*susy)("HMIX", 31)+(*susy)("HMIX", 32)*(*susy)("EXTPAR",25))*(-(*susy)("HMIX", 32)+(*sus_param).epsilonb*(*susy)("HMIX", 31))*pow((*sus_param).mass_b_muW/(*sus_param).mass_H03,2.)

		+((*susy)("AMIX", 11)+(*susy)("AMIX", 12)*(*susy)("EXTPAR",25))*(-(*susy)("AMIX", 12)+(*sus_param).epsilonb*(*susy)("AMIX", 11))*pow((*sus_param).mass_b_muW/(*susy)("MASS",36),2.) //mass_A0 = 36 ? = HO3 ?
		+((*susy)("AMIX", 21)+(*susy)("AMIX", 22)*(*susy)("EXTPAR",25))*(-(*susy)("AMIX", 22)+(*sus_param).epsilonb*(*susy)("AMIX", 21))*pow((*sus_param).mass_b_muW/(*sus_param).mass_A02,2.));
		C8Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("EXTPAR",25),2.)/(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25))/(1.+(*sus_param).epsilon0*(*susy)("EXTPAR",25))*F8_2((*sus_param).yt);
		C8Heps2_0+=-3.*(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("EXTPAR",25),2.)*(1.+pow((*susy)("EXTPAR",25),2.))/(1.+(*sus_param).epsilon0*(*susy)("EXTPAR",25))/72.
		*(((*susy)("HMIX", 11)+(*susy)("HMIX", 12)*(*susy)("EXTPAR",25))*(-(*susy)("HMIX", 12)+(*sus_param).epsilonb*(*susy)("HMIX", 11))*pow((*sus_param).mass_b_muW/(*susy)("MASS",25),2.)
		+((*susy)("HMIX", 21)+(*susy)("HMIX", 22)*(*susy)("EXTPAR",25))*(-(*susy)("HMIX", 22)+(*sus_param).epsilonb*(*susy)("HMIX", 21))*pow((*sus_param).mass_b_muW/(*susy)("MASS",35),2.)
		+((*susy)("HMIX", 31)+(*susy)("HMIX", 32)*(*susy)("EXTPAR",25))*(-(*susy)("HMIX", 32)+(*sus_param).epsilonb*(*susy)("HMIX", 31))*pow((*sus_param).mass_b_muW/(*sus_param).mass_H03,2.)		

+((*susy)("AMIX", 11)+(*susy)("AMIX", 12)*(*susy)("EXTPAR",25))*(-(*susy)("AMIX", 12)+(*sus_param).epsilonb*(*susy)("AMIX", 11))*pow((*sus_param).mass_b_muW/(*susy)("MASS",36),2.)
		+((*susy)("AMIX", 21)+(*susy)("AMIX", 22)*(*susy)("EXTPAR",25))*(-(*susy)("AMIX", 22)+(*sus_param).epsilonb*(*susy)("AMIX", 21))*pow((*sus_param).mass_b_muW/(*sus_param).mass_A02,2.));
		}


	double C7H_0=1./3.*(*sus_param).lu*(*sus_param).lu*F7_1((*sus_param).yt) - (*sus_param).lu*(*sus_param).ld*F7_2((*sus_param).yt);
	double C8H_0=1./3.*(*sus_param).lu*(*sus_param).lu*F8_1((*sus_param).yt) - (*sus_param).lu*(*sus_param).ld*F8_2((*sus_param).yt);

	double C9H_0=(1.-4.*(*sus_param).sw2)/(*sus_param).sw2*C9llH0((*sus_param).xt,(*sus_param).yt,(*sus_param).lu)-D9H0((*sus_param).yt,(*sus_param).lu);
	double C10H_0=-C9llH0((*sus_param).xt,(*sus_param).yt,(*sus_param).lu)/(*sus_param).sw2;


	auto calculateContribution = [&](auto hFunc, const Array3D_3x7x4& X, int ie, int ae, bool isChargeps) -> double {
		double ratio = std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2);
		double msqOverMchSquared = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
		double factor = isChargeps ? (-(*sus_param).epsilonb / (1.0 + (*sus_param).epsilonb * (*susy)("EXTPAR",25)) * (*susy)("EXTPAR",25)) : 1.0;
		return ratio * (
			X[ie][ae][1] * X[ie][ae][2] * hFunc(msqOverMchSquared) +
			(*sus_param).Mch[ie] / (*sus_param).mass_b_muW * X[ie][ae][1] * X[ie][ae][2] * hFunc(msqOverMchSquared)
		) * (*sus_param).kappaFactor * factor;
	};

	double C7charg_0 = 0.0;
	double C8charg_0 = 0.0;
	double C7_chargeps_0 = 0.0;
	double C8_chargeps_0 = 0.0;

	for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) {
			C7charg_0 += calculateContribution(h10, (*sus_param).X_UL, ie, ae, false) + calculateContribution(h20, (*sus_param).X_UR, ie, ae, false);
			C8charg_0 += calculateContribution(h50, (*sus_param).X_UL, ie, ae, false) + calculateContribution(h60, (*sus_param).X_UR, ie, ae, false);
			C7_chargeps_0 += calculateContribution(h20, (*sus_param).X_UR, ie, ae, true);
			C8_chargeps_0 += calculateContribution(h60, (*sus_param).X_UR, ie, ae, true);
		}
	}

    double C9charg_0 = (1.0 - 4.0 * (*sus_param).sw2) / (*sus_param).sw2 * (*sus_param).C90c - (*sus_param).B90c / (*sus_param).sw2 - (*sus_param).D90c;
    double C10charg_0 = ((*sus_param).B100c - (*sus_param).C90c) / (*sus_param).sw2;

	double C1squark_2 = 0.0;

	if ((*sus_param).test) {
		C1squark_2 = -208.0 / 3.0;
		for (int ae = 1; ae <= 6; ++ae) { // Continue avec l'indexation à partir de 1
			double xsqaU = std::pow((*sus_param).MsqU[ae] / (*sm)("MASS", 24), 2.0);
			double xsqaD = std::pow((*sus_param).MsqD[ae] / (*sm)("MASS", 24), 2.0);

			// Ajoute les contributions de MsqU et MsqD séparément
			C1squark_2 += -2.0 * std::pow(4.0 * xsqaU - 1.0, 1.5) * Cl2(2.0 * std::asin(0.5 / std::sqrt(xsqaU))) + 8.0 * (xsqaU - 1.0 / 3.0) * std::log(xsqaU) + 16.0 * xsqaU;
			C1squark_2 += -2.0 * std::pow(4.0 * xsqaD - 1.0, 1.5) * Cl2(2.0 * std::asin(0.5 / std::sqrt(xsqaD))) + 8.0 * (xsqaD - 1.0 / 3.0) * std::log(xsqaD) + 16.0 * xsqaD;
		}
	} else {
		C1squark_2 = 0.0;
	}

	if (C_match.empty()) C_match.resize(1);
	auto& C_LO = C_match[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));

	C_LO[static_cast<size_t>(WilsonCoefficient::C2)] = std::complex<double>(0, 0); // Exemple pour C2
	C_LO[static_cast<size_t>(WilsonCoefficient::C7)] = std::complex<double>(C7SMeps_0 + C7H_0 + C7Heps_0 + C7Heps2_0 + C7charg_0 + C7_chargeps_0, 0);
	C_LO[static_cast<size_t>(WilsonCoefficient::C8)] = std::complex<double>(C8SMeps_0 + C8H_0 + C8Heps_0 + C8Heps2_0 + C8charg_0 + C8_chargeps_0, 0);
	C_LO[static_cast<size_t>(WilsonCoefficient::C9)] = std::complex<double>(C9H_0 + C9charg_0, 0);
	C_LO[static_cast<size_t>(WilsonCoefficient::C10)] = std::complex<double>(C10H_0 + C10charg_0, 0);

	// Affichage pour vérification, si nécessaire
	std::cout << "C7: " << C_LO[static_cast<size_t>(WilsonCoefficient::C7)] << std::endl;
	std::cout << "C8: " << C_LO[static_cast<size_t>(WilsonCoefficient::C8)] << std::endl;

}


// void SUSY_NLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

// 	EpsilonCalculator* epsi = EpsilonCalculator::GetInstance();

// 	Parameters* susy = Parameters::GetInstance(1);
// 	Parameters* sm = Parameters::GetInstance(0);
// 	susy_parameters* sus_param = susy_parameters::GetInstance(scale);

// 	double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale,  (*sm)("MASS",6),(*sm)("MASS",5)); //mass top at top ?
// 	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale,  (*sm)("MASS",6), (*sm)("MASS",5)); //mass bottom 6 (at pole)

// 	double C4charg_1=0.;		
// 	for(int ie=1;ie<=2;ie++) for(int ae=1;ae<=6;ae++) C4charg_1+= pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UL[ie][ae][3]*h40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));	
// 	C4charg_1*=(*sus_param).kappa; 


// 	double C7charg_1=0.;
// 	for(int ie=1;ie<=2;ie++) for(int ae=1;ae<=6;ae++) C7charg_1+=pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UL[ie][ae][3]*h11(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))) + (*sus_param).Mch[ie]/mass_b_muW*(*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UR[ie][ae][3]*h21(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))));		
// 	C7charg_1*=-0.5*(*sus_param).kappa; 
	

// 	double C8charg_1=0.;
// 	for(int ie=1;ie<=2;ie++) for(int ae=1;ae<=6;ae++) C8charg_1+=pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UL[ie][ae][3]*h51(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))) + (*sus_param).Mch[ie]/mass_b_muW*(*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UR[ie][ae][3]*h61(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))));	
// 	C8charg_1*=-0.5*(*sus_param).kappa;

// 	complex_t C91f = 0;
// 	for(int ie=1;ie<=2;ie++) for(int ae=1;ae<=6;ae++) for(int ce=1;ce<=6;ce++) for(int de=1;de<=6;de++) for(int fe=1;fe<=3;fe++) for(int ke=1;ke<=6;ke++) C91f+=(*sus_param).P_U[de][ke]*pow((*sus_param).MsqU[ke]/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ke][ce]*(1.+log(pow(scale/(*sus_param).MsqU[ke],2.)))*(*sus_param).X_UL[ie][de][2]*(*sus_param).X_UL[ie][ae][3]*
// f50(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2.))*(*sus_param).Gamma_UL[ce][fe]*(*sus_param).Gamma_UL[ae][fe];
	
// 	complex_t B1f1=0.;
// 	complex_t B1f2=0.;
// 	for(int ie=1;ie<=2;ie++) for(int je=1;je<=2;je++) for(int fe=1;fe<=3;fe++) for(int ae=1;ae<=6;ae++) for(int be=1;be<=6;be++) for(int ce=1;ce<=6;ce++)
// 	{ B1f1+=pow(((*sm)("MASS", 24))/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][be]*pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[be][ce]*(1.+log(pow(scale/(*sus_param).MsqU[be],2.)))
// 	*(*sus_param).X_UL[je][ae][2]*(*sus_param).X_UL[ie][ce][3]*(
// 	0.5*f90(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[fe]/(*sus_param).Mch[ie],2.))*(*sus_param).X_NL[ie][fe][2]*(*sus_param).X_NL[je][fe][2]);
	
// 	B1f2+=pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][be]*pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[be][ce]*(1.+log(pow(scale/(*sus_param).MsqU[be],2.)))
// 	*(*sus_param).X_UL[je][ae][2]*(*sus_param).X_UL[ie][ce][3]*(fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f100(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[fe]/(*sus_param).Mch[ie],2.))*(*sus_param).X_NR[ie][fe][2]*(*sus_param).X_NR[je][fe][2]
// 	);
// 	}

// 	complex_t B91f=(B1f1-B1f2)*2./3.*(*sus_param).kappa/(*sm)("COUPLING",2)/(*sm)("COUPLING",2);
				
// 	complex_t B101f=-(B1f1+B1f2)*2./3.*(*sus_param).kappa/(*sm)("COUPLING",2)/(*sm)("COUPLING",2);

// 	C91f*=(*sus_param).kappa/6.;

// 	complex_t D91f=0.;		
// 	for(int ie=1;ie<=2;ie++) for(int ae=1;ae<=6;ae++) for(int be=1;be<=6;be++) for(int ce=1;ce<=6;ce++) D91f+= pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][be]*pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[be][ce]*(1.+log(pow(scale/(*sus_param).MsqU[be],2.)))*(*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UL[ie][ce][3]*q51(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.));
// 	D91f*=(*sus_param).kappa; 

// 	complex_t C9four_1=(1.-4.*(*sus_param).sw2)/(*sus_param).sw2*C91f-B91f/(*sus_param).sw2-D91f;
// 	complex_t C10four_1=(B101f-C91f)/(*sus_param).sw2;

// 	int test=1;
// 	for(int ae=1;ae<=6;ae++) test=test&&(fabs((*sus_param).MsqU[ae])>(*sm)("MASS", 24)/2.)&&(fabs((*sus_param).MsqD[ae])>(*sm)("MASS", 24)/2.);
	
// 	if(test)
// 	{		
// 		double C1squark_2=-208./3.;
// 		double xsqa;
// 		for(int ae=1;ae<=6;ae++) 
// 		{
// 			xsqa=pow((*sus_param).MsqU[ae]/(*sm)("MASS", 24),2.);
// 			C1squark_2+=-2.*pow(4.*xsqa-1.,1.5)*Cl2(2.*asin(0.5/sqrt(xsqa))) +8.*(xsqa-1./3.)*log(xsqa)+16.*xsqa;
		
// 			xsqa=pow((*sus_param).MsqD[ae]/(*sm)("MASS", 24),2.);
// 			C1squark_2+=-2.*pow(4.*xsqa-1.,1.5)*Cl2(2.*asin(0.5/sqrt(xsqa))) +8.*(xsqa-1./3.)*log(xsqa)+16.*xsqa;
// 		}
// 	}
// 	else double C1squark_2=0.;

// 	double alphas_mu=(*sm).QCDRunner.runningAlphasCalculation(scale);

// 	if(fabs(C7H_1)*alphas_mu/4./Pi>fabs(C0w[7])) C7H_1*=fabs(C0w[7])/fabs(C7H_1)*4.*Pi/alphas_mu;
// 	if(fabs(C7charg_1)*alphas_mu/4./Pi>fabs(C0w[7])) C7charg_1*=fabs(C0w[7])/fabs(C7charg_1)*4.*Pi/alphas_mu;
// 	if(fabs(C7four_1)*alphas_mu/4./Pi>fabs(C0w[7])) C7four_1*=fabs(C0w[7])/fabs(C7four_1)*4.*Pi/alphas_mu;

// 	if(fabs(C8H_1)*alphas_mu/4./Pi>fabs(C0w[8])) C8H_1*=fabs(C0w[8])/fabs(C8H_1)*4.*Pi/alphas_mu;
// 	if(fabs(C8charg_1)*alphas_mu/4./Pi>fabs(C0w[8])) C8charg_1*=fabs(C0w[8])/fabs(C8charg_1)*4.*Pi/alphas_mu;
// 	if(fabs(C8four_1)*alphas_mu/4./Pi>fabs(C0w[8])) C8four_1*=fabs(C0w[8])/fabs(C8four_1)*4.*Pi/alphas_mu;   

// 	if(fabs(C9H_1)*alphas_mu/4./Pi>fabs(C0w[9])) C9H_1*=fabs(C0w[9])/fabs(C9H_1)*4.*Pi/alphas_mu;
// 	if(fabs(C9charg_1)*alphas_mu/4./Pi>fabs(C0w[9])) C9charg_1*=fabs(C0w[9])/fabs(C9charg_1)*4.*Pi/alphas_mu;
// 	if(fabs(C9four_1)*alphas_mu/4./Pi>fabs(C0w[9])) C9four_1*=fabs(C0w[9])/fabs(C9four_1)*4.*Pi/alphas_mu;

// 	if(fabs(C10H_1)*alphas_mu/4./Pi>fabs(C0w[10])) C10H_1*=fabs(C0w[10])/fabs(C10H_1)*4.*Pi/alphas_mu;
// 	if(fabs(C10charg_1)*alphas_mu/4./Pi>fabs(C0w[10])) C10charg_1*=fabs(C0w[10])/fabs(C10charg_1)*4.*Pi/alphas_mu;
// 	if(fabs(C10four_1)*alphas_mu/4./Pi>fabs(C0w[10])) C10four_1*=fabs(C0w[10])/fabs(C10four_1)*4.*Pi/alphas_mu;


// }


void SUSY_NLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

	auto* epsi = EpsilonCalculator::GetInstance();
	auto* susy = Parameters::GetInstance(1);
	auto* sm_param = Parameters::GetInstance(0);
	auto* sus_param = susy_parameters::GetInstance(scale);

	double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale,  (*sm)("MASS",6),(*sm)("MASS",5)); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale,  (*sm)("MASS",6), (*sm)("MASS",5)); //mass bottom 6 (at pole)
	complex_t C4charg_1 = 0.;
    complex_t C7charg_1 = 0.;
    complex_t C8charg_1 = 0.;
    // Initialisation des variables complexes
	complex_t C91f = 0.0;
	complex_t B1f1 = 0.0;
	complex_t B1f2 = 0.0;
	complex_t D91f = 0.0;
	complex_t B1c1=0.;
	complex_t B1c2=0.;
	complex_t D91c= 0.;
	complex_t C91c = 0.;

	complex_t C7four_1 = 0;
	complex_t C8four_1 = 0;
	// Fusion des boucles pour C91f, B1f1, B1f2
	for (int ie = 1; ie <= 2; ie++) {
		for (int ae = 1; ae <= 6; ae++) {
			double mass24_Mch_ie_squared = pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2.0);
			double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
			double log_mu_W_MsqU_ae = std::log(std::pow(scale / (*sus_param).MsqU[ae], 2.0));

			C4charg_1+= pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UL[ie][ae][3]*h40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
			C7charg_1+=pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UL[ie][ae][3]*h11(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))) + (*sus_param).Mch[ie]/mass_b_muW*(*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UR[ie][ae][3]*h21(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))));
			C8charg_1+=pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UL[ie][ae][3]*h51(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))) + (*sus_param).Mch[ie]/mass_b_muW*(*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UR[ie][ae][3]*h61(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))));
			D91c += std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UL[ie][ae][3] *
                h31(ratio_MsqU_ae_Mch_ie, log_mu_W_MsqU_ae);

			for (int je = 1; je <= 2; je++) {
				double ratio_Mch_je_ie = (*sus_param).Mch[je] / (*sus_param).Mch[ie];
				double factor_abs = 2.0 * std::fabs(ratio_Mch_je_ie);
				double factor_f31_f30 = (f31(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
										4.0 * (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
										(f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
										f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);
				double factor_f41_f40 = (f41(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
										4.0 * (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
										(f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
										f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);

				// Première partie de C91c (pour chaque ie, je, ae)
				C91c += (*sus_param).X_UL[je][ae][2] * (*sus_param).X_UL[ie][ae][3] *
						((factor_abs * factor_f31_f30 * (*susy)("UMIX", je*10+1) * (*susy)("UMIX", ie*10+1)) -
						(factor_f41_f40 * (*susy)("VMIX", je*10+1) * (*susy)("VMIX", ie*10+1)));
			}

			for (int ce = 1; ce <= 6; ce++) {
				for (int fe = 1; fe <= 3; fe++) {
					for (int de = 1; de <= 6; de++) {
						// C91f Part Specific
						for (int ke = 1; ke <= 6; ke++) {
							double MsqU_ke_Mch_ie_squared = pow((*sus_param).MsqU[ke] / (*sus_param).Mch[ie], 2.0);
							double log_scale_MsqU_ke = log(pow(scale / (*sus_param).MsqU[ke], 2.0));
							C91f += (*sus_param).P_U[de][ke] * MsqU_ke_Mch_ie_squared * (*sus_param).P_U[ke][ce] * 
									(1.0 + log_scale_MsqU_ke) * (*sus_param).X_UL[ie][de][2] * (*sus_param).X_UL[ie][ae][3] * 
									f50(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0)) * 
									(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[ae][fe];
						}

						// B1f1 and B1f2 Part Specific for ie, ce, de, fe - using je loop outside
						for (int je = 1; je <= 2; je++) {
							double factor_common = mass24_Mch_ie_squared * (*sus_param).P_U[ae][ce] * (*sus_param).X_UL[je][ae][2] * (*sus_param).X_UL[ie][ce][3];
							B1f1 += factor_common * 0.5 * f90(pow((*sus_param).Mch[je] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).Msn[fe] / (*sus_param).Mch[ie], 2.0)) *
									(*sus_param).X_NL[ie][fe][2] * (*sus_param).X_NL[je][fe][2];

							B1f2 += factor_common * fabs((*sus_param).Mch[je] / (*sus_param).Mch[ie]) * 
									f100(pow((*sus_param).Mch[je] / (*sus_param).Mch[ie], 2.0), 
											pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
											pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
											pow((*sus_param).Msn[fe] / (*sus_param).Mch[ie], 2.0)) *
									(*sus_param).X_NR[ie][fe][2] * (*sus_param).X_NR[je][fe][2];

							double ratio_Mch = (*sus_param).Mch[je] / (*sus_param).Mch[ie];
                        double ratio_MsqU = (*sus_param).MsqU[ae] / (*sus_param).Mch[ie];
                        double ratio_Msn = (*sus_param).Msn[fe] / (*sus_param).Mch[ie];
                        
                        B1c1 += (*sus_param).X_UL[je][ae][2] * (*sus_param).X_UL[ie][ae][3] / pow((*sus_param).Mch[ie], 2) * 
                                (0.5 * (*sus_param).X_NL[ie][fe][2] * (*sus_param).X_NL[je][fe][2] * 
                                (f81(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
                                4 * (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
                                (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
                                 f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
                                 log(pow(scale / (*sus_param).MsqU[ae], 2))));

                        B1c2 += (*sus_param).X_UL[je][ae][2] * (*sus_param).X_UL[ie][ae][3] / pow((*sus_param).Mch[ie], 2) * 
                                ((*sus_param).X_NR[ie][fe][2] * (*sus_param).X_NR[je][fe][2] * fabs(ratio_Mch) * 
                                (f91(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
                                4 * (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
                                (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
                                 f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
                                 log(pow(scale / (*sus_param).MsqU[ae], 2))));
						}
					}

					C91c += (*sus_param).X_UL[ie][ce][2] * (*sus_param).X_UL[ie][ae][3] *
                        (f51(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) +
                         4.0 * ((f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) + 
                         (f40(ratio_MsqU_ae_Mch_ie * 1.0001, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) - 
                          f40(ratio_MsqU_ae_Mch_ie * 0.9999, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0))) / 0.0002)) * log_mu_W_MsqU_ae) *
                        (*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[ae][fe];
				}
				// D91f Part Specific for ae, ce
				double MsqU_ce_Mch_ie_squared = pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0);
				double log_scale_MsqU_ce = 1.0 + log(pow(scale / (*sus_param).MsqU[ce], 2.0));
				D91f += mass24_Mch_ie_squared * (*sus_param).P_U[ae][ce] * MsqU_ce_Mch_ie_squared * log_scale_MsqU_ce * 
						(*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UL[ie][ce][3] * 
						q51(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), MsqU_ce_Mch_ie_squared);

				double log_mu_W_MsqU_ce = std::log(std::pow(scale / (*sus_param).MsqU[ce], 2.0));
				double MsqU_be_Mch_ie_ratio = (*sus_param).MsqU[ce] / (*sus_param).Mch[ie];
				
				for (int de = 1; de <= 6; de++) {
					// Deuxième partie des calculs pour C91c et D91c, et calculs pour C7four_1 et C8four_1
					double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
					double ratio_MsqU_de_Mch_ie = std::pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0);
					

					// Calculs pour C7four_1 et C8four_1
					C7four_1 += std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).P_U[ae][ce] * MsqU_be_Mch_ie_ratio * (*sus_param).P_U[ce][de] * (1.0 + log_mu_W_MsqU_ce) *
								((*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UL[ie][de][3] * (-q11(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) + 2.0 / 3.0 * q21(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie)) +
								(*sus_param).Mch[ie] / (*sus_param).mass_b_muW * (*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UR[ie][de][3] * (-q31(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) + 2.0 / 3.0 * q41(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie)));

					C8four_1 += std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).P_U[ae][ce] * MsqU_be_Mch_ie_ratio * (*sus_param).P_U[ce][de] * (1.0 + log_mu_W_MsqU_ce) *
								((*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UL[ie][de][3] * q21(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) +
								(*sus_param).Mch[ie] / (*sus_param).mass_b_muW * (*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UR[ie][de][3] * q41(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie));
				}
			}
		}
	}

	C4charg_1*=(*sus_param).kappa;
	C7charg_1*=-0.5*(*sus_param).kappa;
	C8charg_1*=-0.5*(*sus_param).kappa;

	C91c *= -(*sus_param).kappa / 8.0;
	D91c *= (*sus_param).kappa;
	complex_t B91c = -(B1c1 - B1c2) * (*sus_param).kappa * (*sm)("MASS", 24) * (*sm)("MASS", 24) / 2.0 / pow((*sm)("COUPLING", 2), 2);
	complex_t B101c = (B1c1 + B1c2) * (*sus_param).kappa * (*sm)("MASS", 24) * (*sm)("MASS", 24) / 2.0 / pow((*sm)("COUPLING", 2), 2);

	C91f *= (*sus_param).kappa / 6.0;
	B1f1 *= (*sus_param).kappa;
	B1f2 *= (*sus_param).kappa;
	D91f *= (*sus_param).kappa;

	// B91f and B101f final calculation
	complex_t B91f = (B1f1 - B1f2) * 2.0 / 3.0 * (*sus_param).kappa / pow((*sm)("COUPLING", 2), 2);
	complex_t B101f = -(B1f1 + B1f2) * 2.0 / 3.0 * (*sus_param).kappa / pow((*sm)("COUPLING", 2), 2);

	

    std::complex<double> C9four_1 = (1. - 4. * (*sus_param).sw2) / (*sus_param).sw2 * C91f - B91f / (*sus_param).sw2 - D91f;
    std::complex<double> C10four_1 = (B101f - C91f) / (*sus_param).sw2;	// Testez pour C1squark_2


	complex_t C9charg_1=(1.-4.*(*sus_param).sw2)/(*sus_param).sw2*C91c-B91c/(*sus_param).sw2-D91c;
	complex_t C10charg_1=(B101c-C91c)/(*sus_param).sw2;


	

	complex_t C4H_1=EH((*sus_param).yt,(*sus_param).lu);
	
 	complex_t C7H_1= G7H((*sus_param).yt,(*sus_param).lu,(*sus_param).ld)+Delta7H((*sus_param).yt,(*sus_param).lu,(*sus_param).ld)*log(pow(scale/((*susy)("MASS", 25)),2.))-4./9.*C4H_1; // param->mass_H -> 25 ?
	complex_t C8H_1= G8H((*sus_param).yt,(*sus_param).lu,(*sus_param).ld)+Delta8H((*sus_param).yt,(*sus_param).lu,(*sus_param).ld)*log(pow(scale/(*susy)("MASS", 25),2.))-1./6.*C4H_1;
	complex_t C9H_1=(1.-4.*(*sus_param).sw2)/(*sus_param).sw2*C9llH1((*sus_param).xt,(*sus_param).yt,(*sus_param).lu,log(pow(scale/(*susy)("MASS", 25),2.)))-D9H1((*sus_param).yt,(*sus_param).lu,log(pow(scale/(*susy)("MASS", 25),2.)));
	complex_t C10H_1=-C9llH1((*sus_param).xt,(*sus_param).yt,(*sus_param).lu,log(pow(scale/(*susy)("MASS", 25),2.)))/(*sus_param).sw2;

	double alphas_mu = sm_param->QCDRunner.runningAlphasCalculation(scale);
	

	SUSY_LO_Strategy::init(sm, scale, C_match);

	if (C_match.size() < 2) C_match.resize(2);
    auto& C_LO = C_match[0]; // Coefficients à l'ordre LO
    auto& C_NLO = C_match[1]; // Coefficients à l'ordre NLO
    if (C_NLO.empty()) C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));

	auto adjustCoefficient = [&](std::complex<double>& Cx_NLO, int index) {
        double ratio = alphas_mu / (4.0 * Pi);
        double absCx_NLO = std::abs(Cx_NLO) * ratio;
        if (absCx_NLO > std::abs(C_LO[index])) {
            Cx_NLO *= std::abs(C_LO[index]) / std::abs(Cx_NLO) * (1.0 / ratio);
        }
    };

	adjustCoefficient(C7H_1, 7);
	adjustCoefficient(C7charg_1, 7);
	adjustCoefficient(C7four_1, 7);

	adjustCoefficient(C8H_1, 8);
	adjustCoefficient(C8charg_1, 8);
	adjustCoefficient(C8four_1, 8);

	adjustCoefficient(C9H_1, 9);
	adjustCoefficient(C9charg_1, 9);
	adjustCoefficient(C9four_1, 9);

	adjustCoefficient(C10H_1, 10);
	adjustCoefficient(C10charg_1, 10);
	adjustCoefficient(C10four_1, 10);

	C_NLO[1] += 0.;
    C_NLO[4] += C4H_1 + C4charg_1;
    C_NLO[7] += C7H_1 + C7charg_1 + C7four_1;
    C_NLO[8] += C8H_1 + C8charg_1 + C8four_1;
    C_NLO[9] += C9H_1 + C9charg_1 + C9four_1;
    C_NLO[10] += C10H_1 + C10charg_1 + C10four_1;
}


void SUSY_NNLO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

	auto* epsi = EpsilonCalculator::GetInstance();
	auto* susy = Parameters::GetInstance(1);
	auto* sm_param = Parameters::GetInstance(0);
	auto* sus_param = susy_parameters::GetInstance(scale);

	complex_t C3charg_2 = 0.0;
	complex_t C4charg_2 = 0.0;
	complex_t C4four_2 = 0.0;
	complex_t C4charg_1 = 0.0;

	for(int ie = 1; ie <= 2; ie++) {
		for(int ae = 1; ae <= 6; ae++) {
			double ratio_mass_W_Mch = std::pow((*sm)("MASS",24)/ (*sus_param).Mch[ie], 2.0);
			double ratio_MsqU_Mch = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
			double log_mu_W_MsqU = std::log(std::pow(scale / (*sus_param).MsqU[ae], 2.0));

			// Calculs pour C3charg_2 et C4charg_2
			C3charg_2 += ratio_mass_W_Mch * (*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UL[ie][ae][3] * h71(ratio_MsqU_Mch, log_mu_W_MsqU);
			C4charg_2 += ratio_mass_W_Mch * (*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UL[ie][ae][3] * h41(ratio_MsqU_Mch, log_mu_W_MsqU);
			C4charg_1+= pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][2]*(*sus_param).X_UL[ie][ae][3]*h40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
			for(int be = 1; be <= 6; be++) {
				for(int ce = 1; ce <= 6; ce++) {
					// Calcul spécifique pour C4four_2 qui nécessite des boucles imbriquées supplémentaires
					C4four_2 += ratio_mass_W_Mch * (*sus_param).P_U[ae][be] * (*sus_param).MsqU[be] / (*sus_param).Mch[ie] * (*sus_param).P_U[be][ce] *
								(1.0 + log_mu_W_MsqU) * (*sus_param).X_UL[ie][ae][2] * (*sus_param).X_UL[ie][ce][3] *
								q61(ratio_MsqU_Mch, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0));
				}
			}
		}
	}

	C3charg_2 *= (*sus_param).kappa;
	C4charg_2 *= (*sus_param).kappa;
	C4four_2 *= (*sus_param).kappa;

	complex_t C5charg_2 = -C3charg_2 / 10.0 + 2.0 / 15.0 * C4charg_1;
    complex_t C6charg_2 = -3.0 / 16.0 * C3charg_2 + 1.0 / 4.0 * C4charg_1;

	double C1squark_2 = 0.0;
	if (std::all_of(begin((*sus_param).MsqU), end((*sus_param).MsqU), [&](double m) { return std::abs(m) > (*sm_param)("MASS", 24) / 2.0; })) {
		C1squark_2 = -208.0 / 3.0;
		for (int ae = 1; ae <= 6; ++ae) {
			double xsqa = std::pow((*sus_param).MsqU[ae] / (*sm_param)("MASS", 24), 2.0);
			if (4.0 * xsqa > 1.0) {
				double angle = 2.0 * asin(0.5 / sqrt(xsqa));
				C1squark_2 += -2.0 * std::pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
			}
			C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;

			xsqa = std::pow((*sus_param).MsqD[ae] / (*sm_param)("MASS", 24), 2.0);
			if (4.0 * xsqa > 1.0) {
				double angle = 2.0 * asin(0.5 / sqrt(xsqa));
				C1squark_2 += -2.0 * std::pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
			}
			C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;
		}
	}

	complex_t C4H_1=EH((*sus_param).yt,(*sus_param).lu);
	complex_t C3H_2=G3H((*sus_param).yt,(*sus_param).lu)+Delta3H((*sus_param).yt,(*sus_param).lu)*log(pow(scale/(*susy)("MASS",25),2.));
	complex_t C4H_2=G4H((*sus_param).yt,(*sus_param).lu)+Delta4H((*sus_param).yt,(*sus_param).lu)*log(pow(scale/(*susy)("MASS",25),2.));
	complex_t C5H_2=-C3H_2/10.+2./15.*C4H_1;
	complex_t C6H_2=-3./16.*C3H_2+1./4.*C4H_1;
	
	complex_t C7H_2=C7H2((*sus_param).yt,(*sus_param).lu,(*sus_param).ld,log(pow(scale/(*sus_param).mass_top_muW,2.)));
	complex_t C8H_2=C8H2((*sus_param).yt,(*sus_param).lu,(*sus_param).ld,log(pow(scale/(*sus_param).mass_top_muW,2.)));

	if (C_match.size() < 3) C_match.resize(3);
    auto& C_NNLO = C_match[2]; // Coefficients à l'ordre NLO
    if (C_NNLO.empty()) C_NNLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));

	C_NNLO[1] += C1squark_2;
	C_NNLO[3] += C3H_2+C3charg_2;
    C_NNLO[4] += C4H_2+C4charg_2+C4four_2;
	C_NNLO[5] += C5H_2+C5charg_2;
	C_NNLO[6] += C6H_2+C6charg_2;
    C_NNLO[7] += C7H_2;
    C_NNLO[8] += C8H_2;

}
