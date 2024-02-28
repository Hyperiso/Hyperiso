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
