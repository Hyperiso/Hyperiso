#include "Wilson_susy.h"


void SUSY_LO_Strategy::init(Parameters* sm, double scale, WilsonSet& C_match) {

	EpsilonCalculator* epsi = EpsilonCalculator::GetInstance();
    // Parameters* sm = Parameters::GetInstance();
    double epsilonbp,epsilon0p,epsilon0,epsilon2,epsilon1p,epsilonb;
	
	epsilonbp=(*epsi).epsilon_bp();
	epsilon0p=(*epsi).epsilon_0p();
	epsilon0=(*epsi).epsilon_0();
	epsilon2=(*epsi).epsilon_2();
	epsilon1p=(*epsi).epsilon_1p();
	epsilonb=epsilon0+epsilon2;

	double mass_top_muW=(*sm).run.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale,  (*sm)("MASS",6),(*sm)("MASS",5));
	double mass_b_muW=(*sm).run.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale,  (*sm)("MASS",6), (*sm)("MASS",5)); //mass bottom 6 (at pole)

	double L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
 	double sw2=pow(sin(atan((*sm)("Coupling",1)/(*sm)("Coupling",2))),2.); //1 = param-> gp and 2 = param->g2

	double xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	double yt= pow(mass_top_muW/(*sm)("MASS",25),2.); // param->mass_H (25)

	double C7SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F7_2(xt);
	double C8SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F8_2(xt);


	double C7Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F7_2(yt);
	double C8Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F8_2(yt);

	double C7Heps2_0=0.;
	double C8Heps2_0=0.;

	if((param->mass_A02==0.)&&(param->mass_H03==0.))
	{
		C7Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F7_2(yt);
		C7Heps2_0+=epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.		*((cos(param->alpha)+sin(param->alpha)*param->tan_beta)*(-sin(param->alpha)+epsilonb*cos(param->alpha))*pow(mass_b_muW/param->mass_h0,2.)
		+(sin(param->alpha)-cos(param->alpha)*param->tan_beta)*(cos(param->alpha)+epsilonb*sin(param->alpha))*pow(mass_b_muW/param->mass_H0,2.)			+(-cos(atan(param->tan_beta))-sin(atan(param->tan_beta))*param->tan_beta)*(sin(atan(param->tan_beta))-epsilonb*cos(atan(param->tan_beta)))*pow(mass_b_muW/param->mass_A0,2.));

		C8Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F8_2(yt);
		C8Heps2_0+=epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.		*((cos(param->alpha)+sin(param->alpha)*param->tan_beta)*(-sin(param->alpha)+epsilonb*cos(param->alpha))*pow(mass_b_muW/param->mass_h0,2.)
		+(sin(param->alpha)-cos(param->alpha)*param->tan_beta)*(cos(param->alpha)+epsilonb*sin(param->alpha))*pow(mass_b_muW/param->mass_H0,2.)			+(-cos(atan(param->tan_beta))-sin(atan(param->tan_beta))*param->tan_beta)*(sin(atan(param->tan_beta))-epsilonb*cos(atan(param->tan_beta)))*pow(mass_b_muW/param->mass_A0,2.));
	}
	else
	{		
		C7Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F7_2(yt);
		C7Heps2_0+=epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.	*((param->H0_mix[1][1]+param->H0_mix[1][2]*param->tan_beta)*(-param->H0_mix[1][2]+epsilonb*param->H0_mix[1][1])*pow(mass_b_muW/param->mass_h0,2.)
		+(param->H0_mix[2][1]+param->H0_mix[2][2]*param->tan_beta)*(-param->H0_mix[2][2]+epsilonb*param->H0_mix[2][1])*pow(mass_b_muW/param->mass_H0,2.)
		+(param->H0_mix[3][1]+param->H0_mix[3][2]*param->tan_beta)*(-param->H0_mix[3][2]+epsilonb*param->H0_mix[3][1])*pow(mass_b_muW/param->mass_H03,2.)

		+(param->A0_mix[1][1]+param->A0_mix[1][2]*param->tan_beta)*(-param->A0_mix[1][2]+epsilonb*param->A0_mix[1][1])*pow(mass_b_muW/param->mass_A0,2.)
		+(param->A0_mix[2][1]+param->A0_mix[2][2]*param->tan_beta)*(-param->A0_mix[2][2]+epsilonb*param->A0_mix[2][1])*pow(mass_b_muW/param->mass_A02,2.));
		C8Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F8_2(yt);
		C8Heps2_0+=-3.*epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.
		*((param->H0_mix[1][1]+param->H0_mix[1][2]*param->tan_beta)*(-param->H0_mix[1][2]+epsilonb*param->H0_mix[1][1])*pow(mass_b_muW/param->mass_h0,2.)
		+(param->H0_mix[2][1]+param->H0_mix[2][2]*param->tan_beta)*(-param->H0_mix[2][2]+epsilonb*param->H0_mix[2][1])*pow(mass_b_muW/param->mass_H0,2.)
		+(param->H0_mix[3][1]+param->H0_mix[3][2]*param->tan_beta)*(-param->H0_mix[3][2]+epsilonb*param->H0_mix[3][1])*pow(mass_b_muW/param->mass_H03,2.)		

+(param->A0_mix[1][1]+param->A0_mix[1][2]*param->tan_beta)*(-param->A0_mix[1][2]+epsilonb*param->A0_mix[1][1])*pow(mass_b_muW/param->mass_A0,2.)
		+(param->A0_mix[2][1]+param->A0_mix[2][2]*param->tan_beta)*(-param->A0_mix[2][2]+epsilonb*param->A0_mix[2][1])*pow(mass_b_muW/param->mass_A02,2.));
		}

	double lu=1./param->tan_beta;
	double ld=-param->tan_beta;

	double C7H_0=1./3.*lu*lu*F7_1(yt) - lu*ld*F7_2(yt);
	double C8H_0=1./3.*lu*lu*F8_1(yt) - lu*ld*F8_2(yt);

	double C9H_0=(1.-4.*sw2)/sw2*C9llH0(xt,yt,lu)-D9H0(yt,lu);
	double C10H_0=-C9llH0(xt,yt,lu)/sw2;


	/* ...........................................................*/

	double C4charg_1,C4charg_2;
	double C3charg_2,C5charg_2,C6charg_2;
	double C7charg_0,C8charg_0,C7_chargeps_0,C8_chargeps_0,C7charg_1,C8charg_1;
	double C9charg_0,C9charg_1,C10charg_0,C10charg_1;
	double C7four_1,C8four_1,C9four_1,C10four_1,C4four_2;
	double C1squark_2;
	
	double Gamma_UL[7][4],Gamma_UR[7][4],Gamma_NL[4][4],Gamma_NR[4][4];
	double Gamma_U[7][7],I_LR[7][7],P_U[7][7];
	double X_UL[3][7][4],X_UR[3][7][4],X_NL[3][4][4],X_NR[3][4][4];
	double MU[4],MD[4],ME[4],VCKM[4][4],Mch[3],MsqU[7],MsqD[7],Msn[4],mintmp;
	double kappa,ag,aY,cosb,sinb,st,ct,alphas_mg;
	int ae,be,ce,de,ee,fe,ge,je,ke;


	alphas_mg=(*sm).run.runningAlphasCalculation(param->mass_gluino);
	ag=1.-7./12./pi*alphas_mg;
	aY=1.+alphas_mg/4./pi;
	
	kappa=1./(param->g2*param->g2*param->Vtb*param->Vts);
	
	VCKM[1][1]=param->Vud;
	VCKM[1][2]=param->Vus;
	VCKM[1][3]=-(param->Vts*param->Vtb+param->Vcs*param->Vcb)/param->Vus; /* Vub from unitarity */
	VCKM[2][1]=param->Vcd;
	VCKM[2][2]=param->Vcs;
	VCKM[2][3]=param->Vcb;
	VCKM[3][1]=param->Vtd;
	VCKM[3][2]=param->Vts;
	VCKM[3][3]=param->Vtb;
	
	sinb=sin(atan(param->tan_beta));
	cosb=cos(atan(param->tan_beta));
	ct=param->stop_mix[2][2];
	st=param->stop_mix[1][2];
	
	MU[1]=param->mass_u;
	MU[2]=param->mass_c;
	MU[3]=mass_top_muW;

	MD[1]=param->mass_u;
	MD[2]=param->mass_s;
	MD[3]=mass_b_muW;

	ME[1]=param->mass_e;
	ME[2]=param->mass_mu;
	ME[3]=param->mass_tau;

	Mch[1]=param->mass_cha1;
	Mch[2]=param->mass_cha2;
	
	MsqU[1]=param->mass_upl;
	MsqU[2]=param->mass_chl;
	MsqU[3]=param->mass_t1;
	MsqU[4]=param->mass_upr;
	MsqU[5]=param->mass_chr;
	MsqU[6]=param->mass_t2;
	
	Msn[1]=param->mass_nuel;
	Msn[2]=param->mass_numl;
	Msn[3]=param->mass_nutl;

	constexpr double Pi = 3.14159265358979323846;
	constexpr int N_UL_UR = 7, M_UL_UR = 4, N_NL_NR = 4, M_NL_NR = 4, N_Gamma_U = 7, M_Gamma_U = 7, N_X = 3, N_Mch = 3, N_MsqU = 7, N_MsqD = 7, N_Msn = 4;
	
	using Array2D_7x4 = std::array<std::array<double, M_UL_UR>, N_UL_UR>;
	using Array2D_4x4 = std::array<std::array<double, M_NL_NR>, N_NL_NR>;
	using Array2D_7x7 = std::array<std::array<double, M_Gamma_U>, N_Gamma_U>;
	using Array3D_3x7x4 = std::array<Array2D_7x4, N_X>;
	using Array1D_4 = std::array<double, 4>;
	using Array1D_3 = std::array<double, N_Mch>;
	using Array1D_7 = std::array<double, N_MsqU>;

	// Initialisation des variables (utilisation de std::array)
	Array2D_7x4 Gamma_UL, Gamma_UR;
	Array2D_4x4 Gamma_NL, Gamma_NR;
	Array2D_7x7 Gamma_U, I_LR, P_U;
	Array3D_3x7x4 X_UL, X_UR, X_NL, X_NR;
	Array1D_4 MU, MD, ME, VCKM[4]; // VCKM est un tableau de Array1D_4
	Array1D_3 Mch;
	Array1D_7 MsqU, MsqD;
	Array1D_4 Msn;
	double kappa, ag, aY, cosb, sinb, st, ct, alphas_mg;
	double C4charg_1, C4charg_2, C3charg_2, C5charg_2, C6charg_2, C7charg_0, C8charg_0, C7_chargeps_0, C8_chargeps_0, C7charg_1, C8charg_1, C9charg_0, C9charg_1, C10charg_0, C10charg_1, C7four_1, C8four_1, C9four_1, C10four_1, C4four_2, C1squark_2;

	// Exemple de calcul initial (les autres calculs doivent être adaptés de manière similaire)
	alphas_mg = sm->run.runningAlphasCalculation(param->mass_gluino);
	ag = 1.0 - 7.0 / (12.0 * Pi * alphas_mg);
	aY = 1.0 + alphas_mg / (4.0 * Pi);
	kappa = 1.0 / (param->g2 * param->g2 * param->Vtb * param->Vts);


	VCKM[0] = {param->Vud, param->Vus, -(param->Vts * param->Vtb + param->Vcs * param->Vcb) / param->Vus, 0.0}; // Ajout d'un élément pour la compatibilité de taille
	VCKM[1] = {param->Vcd, param->Vcs, param->Vcb, 0.0}; // Idem
	VCKM[2] = {param->Vtd, param->Vts, param->Vtb, 0.0};

	sinb = std::sin(std::atan(param->tan_beta));
	cosb = std::cos(std::atan(param->tan_beta));
	ct = param->stop_mix[1][1]; // Ajustement des indices pour base-0
	st = param->stop_mix[0][1]; // Ajustement des indices pour base-0

	// Initialisation des masses
	MU = {0.0, param->mass_u, param->mass_c, param->mass_top_muW}; // Ajout d'un élément vide pour compatibilité d'indice
	MD = {0.0, param->mass_d, param->mass_s, param->mass_b_muW}; // Correction pour inclure mass_d et ajustement pour base-0
	ME = {0.0, param->mass_e, param->mass_mu, param->mass_tau}; // Ajout d'un élément vide pour compatibilité d'indice

	Mch = {param->mass_cha1, param->mass_cha2, 0.0}; // Ajout d'un élément pour compatibilité de taille
	MsqU = {param->mass_upl, param->mass_chl, param->mass_t1, param->mass_upr, param->mass_chr, param->mass_t2, 0.0}; // Ajout d'un élément pour compatibilité de taille
	Msn = {param->mass_nuel, param->mass_numl, param->mass_nutl, 0.0}; // Ajout d'un élément pour compatibilité de taille

	// Vérification du mélange sU_mix et initialisation conditionnelle de Gamma_UL et Gamma_UR
	bool isNonZeroMix = std::any_of(param->sU_mix.begin(), param->sU_mix.end(), [](const auto& mix) {
		return std::accumulate(mix.begin(), mix.end(), 1.0, std::multiplies<>()) != 0.0;
	});

	if (isNonZeroMix) {
		for (int i = 0; i < 6; ++i) { // Pour N_UL_UR - 1
			for (int j = 0; j < 3; ++j) { // Pour M_UL_UR - 1
				Gamma_UL[i][j] = param->sU_mix[i][j];
				Gamma_UR[i][j] = param->sU_mix[i][j + 3];
			}
		}
	} else {
		// Réinitialisation de Gamma_UL et Gamma_UR à des valeurs par défaut
		for (auto& row : Gamma_UL) {
			row.fill(0.0); // Remplit toute la ligne avec des zéros
		}
		for (auto& row : Gamma_UR) {
			row.fill(0.0); // Remplit toute la ligne avec des zéros
		}

		// Configuration spécifique basée sur les exigences décrites
		Gamma_UL[0][0] = 1.0; // Correspond à Gamma_UL[1][1] = 1 dans l'indexation originale de base-1
		Gamma_UL[1][1] = 1.0; // Correspond à Gamma_UL[2][2] = 1
		Gamma_UL[2][2] = ct;  // ct déjà calculé précédemment
		Gamma_UL[5][2] = -st; // -st, ajustement pour base-0

		Gamma_UR[3][0] = 1.0; // Correspond à Gamma_UR[4][1] = 1 dans l'indexation originale de base-1
		Gamma_UR[4][1] = 1.0; // Correspond à Gamma_UR[5][2] = 1
		Gamma_UR[2][2] = st;  // st
		Gamma_UR[5][2] = ct;  // ct
	}



/* .............................................................................;; */
			
	if(param->sU_mix[1][1]*param->sU_mix[1][2]*param->sU_mix[1][3]*param->sU_mix[1][4]*param->sU_mix[1][5]*param->sU_mix[1][6]!=0.)
	{
	
		for(je=1;je<=5;je++)
		{
			for(ie=je+1;ie<=6;ie++) if (MsqU[je]>MsqU[ie])
			{
				mintmp=MsqU[je];
				MsqU[je]=MsqU[ie];
				MsqU[ie]=mintmp;
			}
		}
	
	
		for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UL[ae][ie]=param->sU_mix[ae][ie];
		for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UR[ae][ie]=param->sU_mix[ae][ie+3];
	}
	else
	{
		Gamma_UL[1][1]=1.;
		Gamma_UL[2][1]=0.;
		Gamma_UL[3][1]=0.;
		Gamma_UL[4][1]=0.;
		Gamma_UL[5][1]=0.;
		Gamma_UL[6][1]=0.;
		Gamma_UL[1][2]=0.;
		Gamma_UL[2][2]=1.;
		Gamma_UL[3][2]=0.;
		Gamma_UL[4][2]=0.;
		Gamma_UL[5][2]=0.;
		Gamma_UL[6][2]=0.;
		Gamma_UL[1][3]=0.;
		Gamma_UL[2][3]=0.;
		Gamma_UL[3][3]=ct;
		Gamma_UL[4][3]=0.;
		Gamma_UL[5][3]=0.;
		Gamma_UL[6][3]=-st;
	
		Gamma_UR[1][1]=0.;
		Gamma_UR[2][1]=0.;
		Gamma_UR[3][1]=0.;
		Gamma_UR[4][1]=1.;
		Gamma_UR[5][1]=0.;
		Gamma_UR[6][1]=0.;
		Gamma_UR[1][2]=0.;
		Gamma_UR[2][2]=0.;
		Gamma_UR[3][2]=0.;
		Gamma_UR[4][2]=0.;
		Gamma_UR[5][2]=1.;
		Gamma_UR[6][2]=0.;
		Gamma_UR[1][3]=0.;
		Gamma_UR[2][3]=0.;
		Gamma_UR[3][3]=st;
		Gamma_UR[4][3]=0.;
		Gamma_UR[5][3]=0.;
		Gamma_UR[6][3]=ct;
	}


	for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++)
	{
		Gamma_U[ae][ie]=Gamma_UL[ae][ie];
		Gamma_U[ae][ie+3]=Gamma_UR[ae][ie];
	}
	
	for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) I_LR[ae][be]=0.;
	for(ae=1;ae<=3;ae++) I_LR[ae][ae]=1.;
	for(ae=4;ae<=6;ae++) I_LR[ae][ae]=-1.;
	
	for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) for(de=1;de<=6;de++) P_U[ae][be]=Gamma_U[ae][ce]*I_LR[ce][de]*Gamma_U[be][de];
	
	for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++) if(ae==be) Gamma_NL[ae][be]=Gamma_NR[ae][be]=1.; else Gamma_NL[ae][be]=Gamma_NR[ae][be]=0.;
			
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
	{
		X_UL[ie][ae][be]=0.;
		for(ce=1;ce<=3;ce++) X_UL[ie][ae][be]+=-param->g2*(ag*param->charg_Vmix[ie][1]*Gamma_UL[ae][ce]
		-aY*param->charg_Vmix[ie][2]*Gamma_UR[ae][ce]*MU[ce]/(sqrt(2.)*param->mass_W*sinb))*VCKM[ce][be];
	
		X_UR[ie][ae][be]=0.;
		for(ce=1;ce<=3;ce++) X_UR[ie][ae][be]+=param->g2*aY*param->charg_Umix[ie][2]*Gamma_UL[ae][ce]*VCKM[ce][be]*MD[be]/(sqrt(2)*param->mass_W*cosb);
	}

	for(ie=1;ie<=2;ie++) for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++)
	{
		X_NL[ie][ae][be]=-param->g2*param->charg_Vmix[ie][1]*Gamma_NL[ae][be];
		X_NR[ie][ae][be]=param->g2*param->charg_Umix[ie][2]*Gamma_NL[ae][be]*ME[be]/(sqrt(2.)*param->mass_W*cosb);
	}

	/* LO */

	C7charg_0=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C7charg_0+=pow(param->mass_W/Mch[ie],2.)*(X_UL[ie][ae][2]*X_UL[ie][ae][3]*h10(pow(MsqU[ae]/Mch[ie],2.)) + Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h20(pow(MsqU[ae]/Mch[ie],2.)));		
	C7charg_0*=-0.5*kappa; 
	

	C8charg_0=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C8charg_0+=pow(param->mass_W/Mch[ie],2.)*(X_UL[ie][ae][2]*X_UL[ie][ae][3]*h50(pow(MsqU[ae]/Mch[ie],2.)) + Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h60(pow(MsqU[ae]/Mch[ie],2.)));		
	C8charg_0*=-0.5*kappa; 
	

	C7_chargeps_0=0.;		
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C7_chargeps_0+= -epsilonb/(1.+epsilonb*param->tan_beta)*param->tan_beta*pow(param->mass_W/Mch[ie],2.)*(Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h20(pow(MsqU[ae]/Mch[ie],2.)));				
	C7_chargeps_0*=-0.5*kappa; 
	

	C8_chargeps_0=0.;		
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C8_chargeps_0+= -epsilonb/(1.+epsilonb*param->tan_beta)*param->tan_beta*pow(param->mass_W/Mch[ie],2.)*(Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h60(pow(MsqU[ae]/Mch[ie],2.)));		
	C8_chargeps_0*=-0.5*kappa; 

	double B0c1=0.;
	double B0c2=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
	{ 	B0c1+=X_UL[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(0.5*X_NL[ie][be][2]*X_NL[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));	 B0c2+=X_UL[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NR[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	}
	
	double B90c=-(B0c1-B0c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
	
	double B100c=(B0c1+B0c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
	
	
	double C90c=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) C90c+=X_UL[je][ae][2]*X_UL[ie][ae][3]*(2.*fabs(Mch[je]/Mch[ie])*f30(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))*param->charg_Umix[je][1]*param->charg_Umix[ie][1] -f40(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))*param->charg_Vmix[je][1]*param->charg_Vmix[ie][1]);
	
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=3;ce++) C90c+=X_UL[ie][be][2]*X_UL[ie][ae][3]*f40(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))*Gamma_UL[be][ce]*Gamma_UL[ae][ce];
	C90c*=-kappa/8.;	
	
	
	double D90c=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) D90c+=pow(param->mass_W/Mch[ie],2.)*X_UL[ie][ae][2]*X_UL[ie][ae][3]*h30(pow(MsqU[ae]/Mch[ie],2.));
	D90c*=kappa;

	C9charg_0=(1.-4.*sw2)/sw2*C90c-B90c/sw2-D90c;
	C10charg_0=(B100c-C90c)/sw2;

	MsqD[1]=param->mass_dnl;
	MsqD[2]=param->mass_stl;
	MsqD[3]=param->mass_b1;
	MsqD[4]=param->mass_dnr;
	MsqD[5]=param->mass_str;
	MsqD[6]=param->mass_b2;
	
	int test=1;
	for(ae=1;ae<=6;ae++) test=test&&(fabs(MsqU[ae])>param->mass_W/2.)&&(fabs(MsqD[ae])>param->mass_W/2.);
	
	if(test)
	{		
		C1squark_2=-208./3.;
		double xsqa;
		for(ae=1;ae<=6;ae++) 
		{
			xsqa=pow(MsqU[ae]/param->mass_W,2.);
			C1squark_2+=-2.*pow(4.*xsqa-1.,1.5)*Cl2(2.*asin(0.5/sqrt(xsqa))) +8.*(xsqa-1./3.)*log(xsqa)+16.*xsqa;
		
			xsqa=pow(MsqD[ae]/param->mass_W,2.);
			C1squark_2+=-2.*pow(4.*xsqa-1.,1.5)*Cl2(2.*asin(0.5/sqrt(xsqa))) +8.*(xsqa-1./3.)*log(xsqa)+16.*xsqa;
		}
	}
	else C1squark_2=0.;


	double alphas_mu=(*sm).run.runningAlphasCalculation(scale);

	// if(cabs(C7H_1)*alphas_mu/4./pi>cabs(C0w[7])) C7H_1*=cabs(C0w[7])/cabs(C7H_1)*4.*pi/alphas_mu;
	// if(cabs(C7charg_1)*alphas_mu/4./pi>cabs(C0w[7])) C7charg_1*=cabs(C0w[7])/cabs(C7charg_1)*4.*pi/alphas_mu;
	// if(cabs(C7four_1)*alphas_mu/4./pi>cabs(C0w[7])) C7four_1*=cabs(C0w[7])/cabs(C7four_1)*4.*pi/alphas_mu;

	// if(cabs(C8H_1)*alphas_mu/4./pi>cabs(C0w[8])) C8H_1*=cabs(C0w[8])/cabs(C8H_1)*4.*pi/alphas_mu;
	// if(cabs(C8charg_1)*alphas_mu/4./pi>cabs(C0w[8])) C8charg_1*=cabs(C0w[8])/cabs(C8charg_1)*4.*pi/alphas_mu;
	// if(cabs(C8four_1)*alphas_mu/4./pi>cabs(C0w[8])) C8four_1*=cabs(C0w[8])/cabs(C8four_1)*4.*pi/alphas_mu;   

	// if(cabs(C9H_1)*alphas_mu/4./pi>cabs(C0w[9])) C9H_1*=cabs(C0w[9])/cabs(C9H_1)*4.*pi/alphas_mu;
	// if(cabs(C9charg_1)*alphas_mu/4./pi>cabs(C0w[9])) C9charg_1*=cabs(C0w[9])/cabs(C9charg_1)*4.*pi/alphas_mu;
	// if(cabs(C9four_1)*alphas_mu/4./pi>cabs(C0w[9])) C9four_1*=cabs(C0w[9])/cabs(C9four_1)*4.*pi/alphas_mu;

	// if(cabs(C10H_1)*alphas_mu/4./pi>cabs(C0w[10])) C10H_1*=cabs(C0w[10])/cabs(C10H_1)*4.*pi/alphas_mu;
	// if(cabs(C10charg_1)*alphas_mu/4./pi>cabs(C0w[10])) C10charg_1*=cabs(C0w[10])/cabs(C10charg_1)*4.*pi/alphas_mu;
	// if(cabs(C10four_1)*alphas_mu/4./pi>cabs(C0w[10])) C10four_1*=cabs(C0w[10])/cabs(C10four_1)*4.*pi/alphas_mu;


	C0w[2]= 0.;
	C0w[7]= C7SMeps_0+C7H_0+C7Heps_0+C7Heps2_0+C7charg_0+C7_chargeps_0;
	C0w[8]= C8SMeps_0+C8H_0+C8Heps_0+C8Heps2_0+C8charg_0+C8_chargeps_0;
	C0w[9]= C9H_0+C9charg_0;
	C0w[10]= C10H_0+C10charg_0;
	// C1w[1]= C1SM_1;
	// C1w[4]= C4SM_1+C4H_1+C4charg_1;
	// C1w[7]= C7SM_1+C7H_1+C7charg_1+C7four_1;
	// C1w[8]= C8SM_1+C8H_1+C8charg_1+C8four_1;
	// C1w[9]= C9SM_1+C9H_1+C9charg_1+C9four_1;
	// C1w[10]= C10SM_1+C10H_1+C10charg_1+C10four_1;
	// C2w[1]= C1SM_2+C1squark_2;
	// C2w[2]= C2SM_2;
	// C2w[3]= C3SM_2+C3H_2+C3charg_2;
	// C2w[4]= C4SM_2+C4H_2+C4charg_2+C4four_2;
	// C2w[5]= C5SM_2+C5H_2+C5charg_2;
	// C2w[6]= C6SM_2+C6H_2+C6charg_2;
	// C2w[7]= C7SM_2;
	// C2w[8]= C8SM_2;
}
