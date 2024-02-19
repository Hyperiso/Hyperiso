#include "Wilson_susy.h"

constexpr int N_UL_UR = 7, M_UL_UR = 4, N_NL_NR = 4, M_NL_NR = 4, N_Gamma_U = 7, M_Gamma_U = 7, N_X = 3, N_Mch = 3, N_MsqU = 7, N_MsqD = 7, N_Msn = 4;

using Array2D_7x4 = std::array<std::array<double, M_UL_UR>, N_UL_UR>;
using Array2D_4x4 = std::array<std::array<double, M_NL_NR>, N_NL_NR>;
using Array2D_7x7 = std::array<std::array<double, M_Gamma_U>, N_Gamma_U>;
using Array3D_3x7x4 = std::array<Array2D_7x4, N_X>;
using Array1D_4 = std::array<double, 4>;
using Array1D_3 = std::array<double, N_Mch>;
using Array1D_7 = std::array<double, N_MsqU>;

// double computeContributions(const Array2D_7x4& X, double (func)(double), int ie, double factor) {
//     double result = 0.0;
//     for (int ae = 0; ae < 6; ++ae) {
//         double msqOverMchSquared = std::pow(MsqU[ae] / Mch[ie], 2.0);
//         result += (X[ie][ae][1] * X[ie][ae][2] * func(msqOverMchSquared)) * factor;
//     }
//     return result;
// }

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

	double C7SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*(*sm)("EXTPAR",25))*(*sm)("EXTPAR",25)*F7_2(xt);
	double C8SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*(*sm)("EXTPAR",25))*(*sm)("EXTPAR",25)*F8_2(xt);


	double C7Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*(*sm)("EXTPAR",25))*(*sm)("EXTPAR",25)*F7_2(yt);
	double C8Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*(*sm)("EXTPAR",25))*(*sm)("EXTPAR",25)*F8_2(yt);

	double C7Heps2_0=0.;
	double C8Heps2_0=0.;

	if((param->mass_A02==0.)&&(param->mass_H03==0.))
	{
		C7Heps2_0=-epsilon2*epsilon1p*pow((*sm)("EXTPAR",25),2.)/(1.+epsilonb*(*sm)("EXTPAR",25))/(1.+epsilon0*(*sm)("EXTPAR",25))*F7_2(yt);
		C7Heps2_0+=epsilon2/pow(1.+epsilonb*(*sm)("EXTPAR",25),2.)*(1.+pow((*sm)("EXTPAR",25),2.))/(1.+epsilon0*(*sm)("EXTPAR",25))/72.		*((cos(param->alpha)+sin(param->alpha)*(*sm)("EXTPAR",25))*(-sin(param->alpha)+epsilonb*cos(param->alpha))*pow(mass_b_muW/param->mass_h0,2.)
		+(sin(param->alpha)-cos(param->alpha)*(*sm)("EXTPAR",25))*(cos(param->alpha)+epsilonb*sin(param->alpha))*pow(mass_b_muW/param->mass_H0,2.)			+(-cos(atan((*sm)("EXTPAR",25)))-sin(atan((*sm)("EXTPAR",25)))*(*sm)("EXTPAR",25))*(sin(atan((*sm)("EXTPAR",25)))-epsilonb*cos(atan((*sm)("EXTPAR",25))))*pow(mass_b_muW/param->mass_A0,2.));

		C8Heps2_0=-epsilon2*epsilon1p*pow((*sm)("EXTPAR",25),2.)/(1.+epsilonb*(*sm)("EXTPAR",25))/(1.+epsilon0*(*sm)("EXTPAR",25))*F8_2(yt);
		C8Heps2_0+=epsilon2/pow(1.+epsilonb*(*sm)("EXTPAR",25),2.)*(1.+pow((*sm)("EXTPAR",25),2.))/(1.+epsilon0*(*sm)("EXTPAR",25))/72.		*((cos(param->alpha)+sin(param->alpha)*(*sm)("EXTPAR",25))*(-sin(param->alpha)+epsilonb*cos(param->alpha))*pow(mass_b_muW/param->mass_h0,2.)
		+(sin(param->alpha)-cos(param->alpha)*(*sm)("EXTPAR",25))*(cos(param->alpha)+epsilonb*sin(param->alpha))*pow(mass_b_muW/param->mass_H0,2.)			+(-cos(atan((*sm)("EXTPAR",25)))-sin(atan((*sm)("EXTPAR",25)))*(*sm)("EXTPAR",25))*(sin(atan((*sm)("EXTPAR",25)))-epsilonb*cos(atan((*sm)("EXTPAR",25))))*pow(mass_b_muW/param->mass_A0,2.));
	}
	else
	{		
		C7Heps2_0=-epsilon2*epsilon1p*pow((*sm)("EXTPAR",25),2.)/(1.+epsilonb*(*sm)("EXTPAR",25))/(1.+epsilon0*(*sm)("EXTPAR",25))*F7_2(yt);
		C7Heps2_0+=epsilon2/pow(1.+epsilonb*(*sm)("EXTPAR",25),2.)*(1.+pow((*sm)("EXTPAR",25),2.))/(1.+epsilon0*(*sm)("EXTPAR",25))/72.	*((param->H0_mix[1][1]+param->H0_mix[1][2]*(*sm)("EXTPAR",25))*(-param->H0_mix[1][2]+epsilonb*param->H0_mix[1][1])*pow(mass_b_muW/param->mass_h0,2.)
		+(param->H0_mix[2][1]+param->H0_mix[2][2]*(*sm)("EXTPAR",25))*(-param->H0_mix[2][2]+epsilonb*param->H0_mix[2][1])*pow(mass_b_muW/param->mass_H0,2.)
		+(param->H0_mix[3][1]+param->H0_mix[3][2]*(*sm)("EXTPAR",25))*(-param->H0_mix[3][2]+epsilonb*param->H0_mix[3][1])*pow(mass_b_muW/param->mass_H03,2.)

		+(param->A0_mix[1][1]+param->A0_mix[1][2]*(*sm)("EXTPAR",25))*(-param->A0_mix[1][2]+epsilonb*param->A0_mix[1][1])*pow(mass_b_muW/param->mass_A0,2.)
		+(param->A0_mix[2][1]+param->A0_mix[2][2]*(*sm)("EXTPAR",25))*(-param->A0_mix[2][2]+epsilonb*param->A0_mix[2][1])*pow(mass_b_muW/param->mass_A02,2.));
		C8Heps2_0=-epsilon2*epsilon1p*pow((*sm)("EXTPAR",25),2.)/(1.+epsilonb*(*sm)("EXTPAR",25))/(1.+epsilon0*(*sm)("EXTPAR",25))*F8_2(yt);
		C8Heps2_0+=-3.*epsilon2/pow(1.+epsilonb*(*sm)("EXTPAR",25),2.)*(1.+pow((*sm)("EXTPAR",25),2.))/(1.+epsilon0*(*sm)("EXTPAR",25))/72.
		*((param->H0_mix[1][1]+param->H0_mix[1][2]*(*sm)("EXTPAR",25))*(-param->H0_mix[1][2]+epsilonb*param->H0_mix[1][1])*pow(mass_b_muW/param->mass_h0,2.)
		+(param->H0_mix[2][1]+param->H0_mix[2][2]*(*sm)("EXTPAR",25))*(-param->H0_mix[2][2]+epsilonb*param->H0_mix[2][1])*pow(mass_b_muW/param->mass_H0,2.)
		+(param->H0_mix[3][1]+param->H0_mix[3][2]*(*sm)("EXTPAR",25))*(-param->H0_mix[3][2]+epsilonb*param->H0_mix[3][1])*pow(mass_b_muW/param->mass_H03,2.)		

+(param->A0_mix[1][1]+param->A0_mix[1][2]*(*sm)("EXTPAR",25))*(-param->A0_mix[1][2]+epsilonb*param->A0_mix[1][1])*pow(mass_b_muW/param->mass_A0,2.)
		+(param->A0_mix[2][1]+param->A0_mix[2][2]*(*sm)("EXTPAR",25))*(-param->A0_mix[2][2]+epsilonb*param->A0_mix[2][1])*pow(mass_b_muW/param->mass_A02,2.));
		}

	double lu=1./(*sm)("EXTPAR",25);
	double ld=-(*sm)("EXTPAR",25);

	double C7H_0=1./3.*lu*lu*F7_1(yt) - lu*ld*F7_2(yt);
	double C8H_0=1./3.*lu*lu*F8_1(yt) - lu*ld*F8_2(yt);

	double C9H_0=(1.-4.*sw2)/sw2*C9llH0(xt,yt,lu)-D9H0(yt,lu);
	double C10H_0=-C9llH0(xt,yt,lu)/sw2;


	/* ...........................................................*/

	constexpr double Pi = 3.14159265358979323846;
	

	// Initialisation des variables (utilisation de std::array)
	Array2D_7x4 Gamma_UL, Gamma_UR;
	Array2D_4x4 Gamma_NL, Gamma_NR;
	Array2D_7x7 Gamma_U, I_LR, P_U;
	Array3D_3x7x4 X_UL, X_UR, X_NL, X_NR;
	Array2D_4x4  VCKM = {{{0.,0.,0.,0.},
					{param->Vud, param->Vus, -(param->Vts * param->Vtb + param->Vcs * param->Vcb) / param->Vus, 0.0},
					{param->Vcd, param->Vcs, param->Vcb, 0.0},
					{param->Vtd, param->Vts, param->Vtb, 0.0}
					}}; 
	Array1D_7 MsqD;
	constexpr size_t NumSquarks = 6;
	std::array<std::array<double, 4>, 7> sU_mix;
	double kappa, ag, aY, cosb, sinb, st, ct, alphas_mg;
	double C4charg_1, C4charg_2, C3charg_2, C5charg_2, C6charg_2, C7charg_0, C8charg_0, C7_chargeps_0, C8_chargeps_0, C7charg_1, C8charg_1, C9charg_0, C9charg_1, C10charg_0, C10charg_1, C7four_1, C8four_1, C9four_1, C10four_1, C4four_2, C1squark_2;

	// Exemple de calcul initial (les autres calculs doivent être adaptés de manière similaire)
	alphas_mg = sm->run.runningAlphasCalculation(param->mass_gluino);
	ag = 1.0 - 7.0 / (12.0 * Pi * alphas_mg);
	aY = 1.0 + alphas_mg / (4.0 * Pi);
	kappa = 1.0 / (param->g2 * param->g2 * param->Vtb * param->Vts);


	sinb = std::sin(std::atan((*sm)("EXTPAR",25)));
	cosb = std::cos(std::atan((*sm)("EXTPAR",25)));
	ct = param->stop_mix[1][1]; // Ajustement des indices pour base-0
	st = param->stop_mix[0][1]; // Ajustement des indices pour base-0

	// Initialisation des masses
	Array1D_4 MU = {0.0, param->mass_u, param->mass_c, param->mass_top_muW}; // Ajout d'un élément vide pour compatibilité d'indice
	Array1D_4 MD = {0.0, param->mass_d, param->mass_s, param->mass_b_muW}; // Correction pour inclure mass_d et ajustement pour base-0
	Array1D_4 ME = {0.0, param->mass_e, param->mass_mu, param->mass_tau}; // Ajout d'un élément vide pour compatibilité d'indice

	Array1D_3 Mch = {0.0,param->mass_cha1, param->mass_cha2 }; // Ajout d'un élément pour compatibilité de taille
	Array1D_7 MsqU = { 0.0, param->mass_upl, param->mass_chl, param->mass_t1, param->mass_upr, param->mass_chr, param->mass_t2}; // Ajout d'un élément pour compatibilité de taille
	Array1D_7 MsqD = {
    0.0, // L'indice 0 est laissé à 0 pour la compatibilité avec l'indexation 1-based utilisée dans votre exemple.
    param->mass_dnl,
    param->mass_stl,
    param->mass_b1,
    param->mass_dnr,
    param->mass_str,
    param->mass_b2
};
	Array1D_4 Msn = {param->mass_nuel, param->mass_numl, param->mass_nutl, 0.0}; // Ajout d'un élément pour compatibilité de taille

	// Vérification du mélange sU_mix et initialisation conditionnelle de Gamma_UL et Gamma_UR
	bool isNonZeroMix = true;
	for (size_t i = 0; i < NumSquarks; ++i) {
		double product = 1.0;
		for (size_t j = 0; j < 6; ++j) { // Pour chaque colonne dans la rangée i
			product *= sU_mix[i][j];
		}
		if (product == 0.0) {
			isNonZeroMix = false;
			break;
		}
	}

	if (isNonZeroMix) {
		// Tri de MsqU si la condition est vraie
		std::sort(MsqU.begin(), MsqU.end());

		// Remplissage des matrices Gamma_UL et Gamma_UR
		for (size_t ae = 0; ae < NumSquarks; ++ae) {
			for (size_t ie = 0; ie < 3; ++ie) { // Indices ajustés pour base-0
				Gamma_UL[ae][ie] = sU_mix[ae][ie];
				Gamma_UR[ae][ie] = sU_mix[ae][ie + 3];
			}
		}
}
	else {
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

	for (int ae = 0; ae < 6; ++ae) {
        for (int ie = 0; ie < 3; ++ie) {
            Gamma_U[ae][ie] = Gamma_UL[ae][ie];
            Gamma_U[ae][ie+3] = Gamma_UR[ae][ie];
        }
    }

    I_LR.fill({});
    for (int i = 0; i < 3; ++i) {
        I_LR[i][i] = 1;
        I_LR[i+3][i+3] = -1;
    }

    // Calcul de P_U
    for (int ae = 0; ae < 6; ++ae) {
        for (int be = 0; be < 6; ++be) {
            for (int ce = 0; ce < 6; ++ce) {
                for (int de = 0; de < 6; ++de) {
                    P_U[ae][be] += Gamma_U[ae][ce] * I_LR[ce][de] * Gamma_U[be][de];
                }
            }
        }
    }

    // Calculs pour X_UL et X_UR
    for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) { // Supposant que 6 est correct pour tous
			for (int be = 0; be < 3; ++be) {
				// Réinitialisation pour X_UL et X_UR
				X_UL[ie][ae][be] = 0.0;
				X_UR[ie][ae][be] = 0.0;

				// Calculs pour X_UL et X_UR
				for (int ce = 0; ce < 3; ++ce) {
					X_UL[ie][ae][be] += -param->g2 * (
						ag * param->charg_Vmix[ie][1] * Gamma_UL[ae][ce] -
						aY * param->charg_Vmix[ie][2] * Gamma_UR[ae][ce] * MU[ce] / (sqrt(2.0) * param->mass_W * sinb)
					) * VCKM[ce][be];

					X_UR[ie][ae][be] += param->g2 * aY * param->charg_Umix[ie][2] * Gamma_UL[ae][ce] * VCKM[ce][be] * MD[be] / (sqrt(2.0) * param->mass_W * cosb);
				}

				// Condition pour éviter le dépassement dans X_NL et X_NR si ae > 2
				if (ae < 3) {
					X_NL[ie][ae][be] = -param->g2 * param->charg_Vmix[ie][1] * Gamma_NL[ae][be];
					X_NR[ie][ae][be] = param->g2 * param->charg_Umix[ie][2] * Gamma_NL[ae][be] * ME[be] / (sqrt(2.0) * param->mass_W * cosb);
				}
			}
		}
	}

	auto computeContributions = [&](int ie, auto func, double additionalFactor = 1.0) {
		double result = 0.0;
		for (int ae = 0; ae < 6; ++ae) {
			double msqOverMchSquared = std::pow(MsqU[ae] / Mch[ie], 2.0);
			result += (X_UL[ie][ae][1] * X_UL[ie][ae][2] * func(msqOverMchSquared) +
					Mch[ie] / mass_b_muW * X_UL[ie][ae][1] * X_UR[ie][ae][2] * func(msqOverMchSquared)) * additionalFactor;
		}
		return result;
	};


	auto hFunc10 = [](double x) { return h10(x); };
	auto hFunc20 = [](double x) { return h20(x); };
	auto hFunc50 = [](double x) { return h50(x); };
	auto hFunc60 = [](double x) { return h60(x); };

	double kappaFactor = -0.5 * kappa;

	auto calculateContribution = [&](auto hFunc, const Array3D_3x7x4& X, int ie, int ae, bool isChargeps) -> double {
		double ratio = std::pow(param->mass_W / Mch[ie], 2);
		double msqOverMchSquared = std::pow(MsqU[ae] / Mch[ie], 2.0);
		double factor = isChargeps ? (-epsilonb / (1.0 + epsilonb * (*sm)("EXTPAR",25)) * (*sm)("EXTPAR",25)) : 1.0;
		return ratio * (
			X[ie][ae][1] * X[ie][ae][2] * hFunc(msqOverMchSquared) +
			Mch[ie] / mass_b_muW * X[ie][ae][1] * X[ie][ae][2] * hFunc(msqOverMchSquared)
		) * kappaFactor * factor;
	};

	C7charg_0 = 0.0;
	C8charg_0 = 0.0;
	C7_chargeps_0 = 0.0;
	C8_chargeps_0 = 0.0;

	for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) {
			C7charg_0 += calculateContribution(h10, X_UL, ie, ae, false) + calculateContribution(h20, X_UR, ie, ae, false);
			C8charg_0 += calculateContribution(h50, X_UL, ie, ae, false) + calculateContribution(h60, X_UR, ie, ae, false);
			C7_chargeps_0 += calculateContribution(h20, X_UR, ie, ae, true);
			C8_chargeps_0 += calculateContribution(h60, X_UR, ie, ae, true);
		}
	}



	double B0c1 = 0.0, B0c2 = 0.0, B90c = 0.0, B100c = 0.0, C90c = 0.0, D90c = 0.0;

	for (int ie = 0; ie < 2; ++ie) {
		for (int je = 0; je < 2; ++je) {
			for (int ae = 0; ae < 6; ++ae) {
				double mchRatioSquared = std::pow(Mch[je] / Mch[ie], 2.0);
				double msqOverMchSquared = std::pow(MsqU[ae] / Mch[ie], 2.0);

				// Calculations for B0c1, B0c2, and part of C90c
				for (int be = 0; be < 3; ++be) {
					double msnOverMchSquared = std::pow(Msn[be] / Mch[ie], 2.0);
					B0c1 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (Mch[ie] * Mch[ie]) * (0.5 * X_NL[ie][be][1] * X_NL[je][be][1] * f50(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));
					B0c2 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (Mch[ie] * Mch[ie]) * (X_NR[ie][be][1] * X_NR[je][be][1] * std::fabs(Mch[je] / Mch[ie]) * f60(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));
					
					if (ie == je) { // Optimization to avoid duplicate calculations
						C90c += X_UL[je][ae][1] * X_UL[ie][ae][2] * (2.0 * std::fabs(Mch[je] / Mch[ie]) * f30(mchRatioSquared, msqOverMchSquared) * param->charg_Umix[je][1] * param->charg_Umix[ie][1] - f40(mchRatioSquared, msqOverMchSquared) * param->charg_Vmix[je][1] * param->charg_Vmix[ie][1]);
					}
				}

				// Part of D90c calculation
				D90c += std::pow(param->mass_W / Mch[ie], 2.0) * X_UL[ie][ae][1] * X_UL[ie][ae][2] * h30(msqOverMchSquared);
			}
		}
	}

	// Additional loop for the rest of C90c that depends on ae and be
	for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) {
			for (int be = 0; be < 6; ++be) { // Notice that both ae and be go up to 6
				double msqOverMchSquaredAe = std::pow(MsqU[ae] / Mch[ie], 2.0);
				double msqOverMchSquaredBe = std::pow(MsqU[be] / Mch[ie], 2.0);
				for (int ce = 0; ce < 3; ++ce) {
					C90c += X_UL[ie][be][1] * X_UL[ie][ae][2] * f40(msqOverMchSquaredAe, msqOverMchSquaredBe) * Gamma_UL[be][ce] * Gamma_UL[ae][ce];
				}
			}
		}
	}

	// Final adjustments
	B90c = -(B0c1 - B0c2) * kappa * std::pow(param->mass_W, 2.0) / (2.0 * std::pow(param->g2, 2.0));
	B100c = (B0c1 + B0c2) * kappa * std::pow(param->mass_W, 2.0) / (2.0 * std::pow(param->g2, 2.0));
	C90c *= -kappa / 8.0;
	D90c *= kappa;

    C9charg_0 = (1.0 - 4.0 * sw2) / sw2 * C90c - B90c / sw2 - D90c;
    C10charg_0 = (B100c - C90c) / sw2;


	bool test = true;
	for (int ae = 1; ae <= 6; ++ae) { // Conserve l'indexation à partir de 1
		if (!(std::fabs(MsqU[ae]) > param->mass_W / 2. && std::fabs(MsqD[ae]) > param->mass_W / 2.)) {
			test = false;
			break; // Sort de la boucle dès qu'une condition n'est pas remplie
		}
	}

	double C1squark_2 = 0.0;

	if (test) {
		C1squark_2 = -208.0 / 3.0;
		for (int ae = 1; ae <= 6; ++ae) { // Continue avec l'indexation à partir de 1
			double xsqaU = std::pow(MsqU[ae] / param->mass_W, 2.0);
			double xsqaD = std::pow(MsqD[ae] / param->mass_W, 2.0);

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
