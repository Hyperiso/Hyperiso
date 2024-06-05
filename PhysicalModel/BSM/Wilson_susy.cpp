#include "Wilson_susy.h"
#include "susy_parameters.h"
#include <sstream>
#include <iomanip>




void SUSY_LO_Strategy::init(double scale, WilsonSet& C_match) {

	Logger *logger = Logger::getInstance();

	EpsilonCalculator* epsi = EpsilonCalculator::GetInstance();

	Parameters* sm = Parameters::GetInstance();
	Parameters* susy = Parameters::GetInstance(1);
	
	susy_parameters* sus_param = susy_parameters::GetInstance(scale);


	complex_t C7SMeps_0= ((*sus_param).epsilonb-(*sus_param).epsilonbp)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))*(*susy)("HMIX",2)*F7_2((*sus_param).xt);
	complex_t C8SMeps_0= ((*sus_param).epsilonb-(*sus_param).epsilonbp)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))*(*susy)("HMIX",2)*F8_2((*sus_param).xt);

	logger->debug("epsilon b : " + std::to_string((*sus_param).epsilonb));
	logger->debug("epsilon bp : " + std::to_string((*sus_param).epsilonbp));

	complex_t C7Heps_0=(-(*sus_param).epsilon0p-(*sus_param).epsilonb)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))*(*susy)("HMIX",2)*F7_2((*sus_param).yt);
	complex_t C8Heps_0=(-(*sus_param).epsilon0p-(*sus_param).epsilonb)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))*(*susy)("HMIX",2)*F8_2((*sus_param).yt);

	complex_t C7Heps2_0=0.;
	complex_t C8Heps2_0=0.;


	if(((*sus_param).mass_A02==0.)&&((*sus_param).mass_H03==0.))
	{
		C7Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("HMIX",2),2.)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))*F7_2((*sus_param).yt);
		C7Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("HMIX",2),2.)*(1.+pow((*susy)("HMIX",2),2.))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))/72.		*((cos((*susy)("ALPHA",0))+sin((*susy)("ALPHA",0))*(*susy)("HMIX",2))*(-sin((*susy)("ALPHA",0))+(*sus_param).epsilonb*cos((*susy)("ALPHA",0)))*pow((*sus_param).mass_b_muW/(*susy)("MASS",25),2.)
		+(sin((*susy)("ALPHA",0))-cos((*susy)("ALPHA",0))*(*susy)("HMIX",2))*(cos((*susy)("ALPHA",0))+(*sus_param).epsilonb*sin((*susy)("ALPHA",0)))*pow((*sus_param).mass_b_muW/(*susy)("MASS",35),2.)			+(-cos(atan((*susy)("HMIX",2)))-sin(atan((*susy)("HMIX",2)))*(*susy)("HMIX",2))*(sin(atan((*susy)("HMIX",2)))-(*sus_param).epsilonb*cos(atan((*susy)("HMIX",2))))*pow((*sus_param).mass_b_muW/(*susy)("MASS",36),2.));

		C8Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("HMIX",2),2.)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))*F8_2((*sus_param).yt);
		C8Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("HMIX",2),2.)*(1.+pow((*susy)("HMIX",2),2.))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))/72.		*((cos((*susy)("ALPHA",0))+sin((*susy)("ALPHA",0))*(*susy)("HMIX",2))*(-sin((*susy)("ALPHA",0))+(*sus_param).epsilonb*cos((*susy)("ALPHA",0)))*pow((*sus_param).mass_b_muW/(*susy)("MASS",25),2.)
		+(sin((*susy)("ALPHA",0))-cos((*susy)("ALPHA",0))*(*susy)("HMIX",2))*(cos((*susy)("ALPHA",0))+(*sus_param).epsilonb*sin((*susy)("ALPHA",0)))*pow((*sus_param).mass_b_muW/(*susy)("MASS",35),2.)			+(-cos(atan((*susy)("HMIX",2)))-sin(atan((*susy)("HMIX",2)))*(*susy)("HMIX",2))*(sin(atan((*susy)("HMIX",2)))-(*sus_param).epsilonb*cos(atan((*susy)("HMIX",2))))*pow((*sus_param).mass_b_muW/(*susy)("MASS",36),2.));
	}
	else
	{		
		C7Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("HMIX",2),2.)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))*F7_2((*sus_param).yt);
		C7Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("HMIX",2),2.)*(1.+pow((*susy)("HMIX",2),2.))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))/72.	*(((*susy)("HMIX", 00)+(*susy)("HMIX", 01)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 01)+(*sus_param).epsilonb*(*susy)("HMIX", 00))*pow((*sus_param).mass_b_muW/(*susy)("MASS",25),2.)
		+((*susy)("HMIX", 10)+(*susy)("HMIX", 11)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 11)+(*sus_param).epsilonb*(*susy)("HMIX", 10))*pow((*sus_param).mass_b_muW/(*susy)("MASS",35),2.)
		+((*susy)("HMIX", 20)+(*susy)("HMIX", 21)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 21)+(*sus_param).epsilonb*(*susy)("HMIX", 20))*pow((*sus_param).mass_b_muW/(*sus_param).mass_H03,2.)

		+((*susy)("AMIX", 00)+(*susy)("AMIX", 01)*(*susy)("HMIX",2))*(-(*susy)("AMIX", 01)+(*sus_param).epsilonb*(*susy)("AMIX", 00))*pow((*sus_param).mass_b_muW/(*susy)("MASS",36),2.) //mass_A0 = 36 ? = HO3 ?
		+((*susy)("AMIX", 10)+(*susy)("AMIX", 11)*(*susy)("HMIX",2))*(-(*susy)("AMIX", 11)+(*sus_param).epsilonb*(*susy)("AMIX", 10))*pow((*sus_param).mass_b_muW/(*sus_param).mass_A02,2.));
		C8Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("HMIX",2),2.)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))*F8_2((*sus_param).yt);
		C8Heps2_0+=-3.*(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("HMIX",2),2.)*(1.+pow((*susy)("HMIX",2),2.))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))/72.
		*(((*susy)("HMIX", 00)+(*susy)("HMIX", 01)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 01)+(*sus_param).epsilonb*(*susy)("HMIX", 00))*pow((*sus_param).mass_b_muW/(*susy)("MASS",25),2.)
		+((*susy)("HMIX", 10)+(*susy)("HMIX", 11)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 11)+(*sus_param).epsilonb*(*susy)("HMIX", 10))*pow((*sus_param).mass_b_muW/(*susy)("MASS",35),2.)
		+((*susy)("HMIX", 31)+(*susy)("HMIX", 21)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 21)+(*sus_param).epsilonb*(*susy)("HMIX", 20))*pow((*sus_param).mass_b_muW/(*sus_param).mass_H03,2.)		

+((*susy)("AMIX", 00)+(*susy)("AMIX", 01)*(*susy)("HMIX",2))*(-(*susy)("AMIX", 01)+(*sus_param).epsilonb*(*susy)("AMIX", 00))*pow((*sus_param).mass_b_muW/(*susy)("MASS",36),2.)
		+((*susy)("AMIX", 10)+(*susy)("AMIX", 11)*(*susy)("HMIX",2))*(-(*susy)("AMIX", 11)+(*sus_param).epsilonb*(*susy)("AMIX", 10))*pow((*sus_param).mass_b_muW/(*sus_param).mass_A02,2.));
		}


	logger->info("truc : " + std::to_string(std::real(C7SMeps_0)));

	auto calculateContribution = [&](auto hFunc, const Array3D_3x7x4& X,const Array3D_3x7x4& X2, int ie, int ae, bool isChargeps) -> double {
		double ratio = std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2);
		
		double msqOverMchSquared = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
		
		double factor = isChargeps ? (-(*sus_param).epsilonb / (1.0 + (*sus_param).epsilonb * (*susy)("HMIX",2)) * (*susy)("HMIX",2)) : 1.0;
		// logger->info("truc1 : " + std::to_string((*sus_param).epsilonb));
		return ratio * (
			X[ie][ae][1] * X2[ie][ae][2] * hFunc(msqOverMchSquared)) * (*sus_param).kappaFactor * factor;
	};
	
	complex_t C7charg_0 = 0.0;
	complex_t C8charg_0 = 0.0;
	complex_t C7_chargeps_0 = 0.0;
	complex_t C8_chargeps_0 = 0.0;

	for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) {
			
			C7charg_0 += calculateContribution(h10, (*sus_param).X_UL,(*sus_param).X_UL, ie, ae, false) + (*sus_param).Mch[ie]/(*sus_param).mass_b_muW *calculateContribution(h20, (*sus_param).X_UL, (*sus_param).X_UR, ie, ae, false);
			C8charg_0 += calculateContribution(h50, (*sus_param).X_UL, (*sus_param).X_UL, ie, ae, false) + (*sus_param).Mch[ie]/(*sus_param).mass_b_muW*calculateContribution(h60,(*sus_param).X_UL,  (*sus_param).X_UR, ie, ae, false);
			C7_chargeps_0 += (*sus_param).Mch[ie]/(*sus_param).mass_b_muW *calculateContribution(h20, (*sus_param).X_UL, (*sus_param).X_UR, ie, ae, true);
			C8_chargeps_0 += (*sus_param).Mch[ie]/(*sus_param).mass_b_muW *calculateContribution(h60, (*sus_param).X_UL,  (*sus_param).X_UR, ie, ae, true);
		}
	}

    complex_t C9charg_0 = (1.0 - 4.0 * (*sus_param).sw2) / (*sus_param).sw2 * (*sus_param).C90c - (*sus_param).B90c / (*sus_param).sw2 - (*sus_param).D90c;
    complex_t C10charg_0 = ((*sus_param).B100c - (*sus_param).C90c) / (*sus_param).sw2;

	logger->info("C90c : " + std::to_string((*sus_param).C90c));
	logger->info("B100c : " + doubleToString((*sus_param).B100c, 20));
	complex_t C1squark_2 = 0.0;

	if ((*sus_param).test) {
		C1squark_2 = -208.0 / 3.0;
		for (int ae = 0; ae < 6; ++ae) { 
			double xsqaU = std::pow((*sus_param).MsqU[ae] / (*sm)("MASS", 24), 2.0);
			double xsqaD = std::pow((*sus_param).MsqD[ae] / (*sm)("MASS", 24), 2.0);

			// Ajoute les contributions de MsqU et MsqD séparément
			C1squark_2 += -2.0 * std::pow(4.0 * xsqaU - 1.0, 1.5) * Cl2(2.0 * std::asin(0.5 / std::sqrt(xsqaU))) + 8.0 * (xsqaU - 1.0 / 3.0) * std::log(xsqaU) + 16.0 * xsqaU;
			C1squark_2 += -2.0 * std::pow(4.0 * xsqaD - 1.0, 1.5) * Cl2(2.0 * std::asin(0.5 / std::sqrt(xsqaD))) + 8.0 * (xsqaD - 1.0 / 3.0) * std::log(xsqaD) + 16.0 * xsqaD;
		}
	} else {
		C1squark_2 = 0.0;
	}
	logger->debug("C7SMeps_0 : " + std::to_string(std::real(C7SMeps_0)));
	logger->debug("C7Heps_0 : " + std::to_string(std::real(C7Heps_0)));
	logger->debug("C7Heps2_0 : " + std::to_string(std::real(C7Heps2_0)));
	logger->debug("C7charg_0 : " + std::to_string(std::real(C7charg_0)));
	logger->debug("C7_chargeps_0 : " + std::to_string(std::real(C7_chargeps_0)));


	std::unique_ptr<THDM_LO_Strategy> thdm_lo = std::make_unique<THDM_LO_Strategy>();

	thdm_lo->set_lu(1/(*susy)("HMIX", 2));
	thdm_lo->set_ld(-(*susy)("HMIX", 2));
	thdm_lo->init(scale, C_match);

	if (C_match.empty()) C_match.resize(1);
	auto& C_LO = C_match[0];
	C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));
	
	logger->debug("CH7 : " + std::to_string(std::real(C_LO[static_cast<size_t>(WilsonCoefficient::C7)])));
	logger->info("CH10 : " + std::to_string(std::real(C_LO[static_cast<size_t>(WilsonCoefficient::C10)])));

	// C_LO[static_cast<size_t>(WilsonCoefficient::C2)] += std::complex<double>(0, 0);
	C_LO[static_cast<size_t>(WilsonCoefficient::C7)] += C7SMeps_0 + C7Heps_0 + C7Heps2_0 + C7charg_0 + C7_chargeps_0;
	C_LO[static_cast<size_t>(WilsonCoefficient::C8)] += C8SMeps_0 +  C8Heps_0 + C8Heps2_0 + C8charg_0 + C8_chargeps_0;
	C_LO[static_cast<size_t>(WilsonCoefficient::C9)] += C9charg_0;
	C_LO[static_cast<size_t>(WilsonCoefficient::C10)] += C10charg_0;

	logger->info("SUSY LO Wilson Coefficient Initialized at scale " +std::to_string(scale)+" terminated successfully");
}


void SUSY_NLO_Strategy::init(double scale, WilsonSet& C_match) {

	auto* epsi = EpsilonCalculator::GetInstance();
	Parameters* sm = Parameters::GetInstance();
	auto* susy = Parameters::GetInstance(1);

	auto* sus_param = susy_parameters::GetInstance(scale);
	Logger *logger = Logger::getInstance();

	double mass_top_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale); //mass top at top ?
	double mass_b_muW=(*sm).QCDRunner.running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale); //mass bottom 6 (at pole)
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
	for (int ie = 0; ie < 2; ie++) {

		for (int ae = 0; ae < 6; ae++) {
			double mass24_Mch_ie_squared = pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2.0);
			double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
			double log_mu_W_MsqU_ae = std::log(std::pow(scale / (*sus_param).MsqU[ae], 2.0));

			
			C4charg_1+= pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
			C7charg_1+=pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h11(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))) + (*sus_param).Mch[ie]/mass_b_muW*(*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h21(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))));
			C8charg_1+=pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h51(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))) + (*sus_param).Mch[ie]/mass_b_muW*(*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h61(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(scale/(*sus_param).MsqU[ae],2.))));
			D91c += std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ae][2] *
                h31(ratio_MsqU_ae_Mch_ie, log_mu_W_MsqU_ae);

			

			for (int je = 0; je < 2; je++) {
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

				C91c += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] *
						((factor_abs * factor_f31_f30 * (*susy)("UMIX", je*10+0) * (*susy)("UMIX", ie*10+0)) -
						(factor_f41_f40 * (*susy)("VMIX", je*10+0) * (*susy)("VMIX", ie*10+0)));

				for (int de=0; de<6; de++) {
					for (int ke=0; ke<6; ke++) {
						C91f+=(*sus_param).P_U[de][ke] * pow((*sus_param).MsqU[ke]/(*sus_param).Mch[ie],2.) * (*sus_param).P_U[ke][ae] * (1+log(pow(scale/(*sus_param).MsqU[ke],2.)))*(*sus_param).X_UL[je][de][1]*(*sus_param).X_UL[ie][ae][2] * (
							2.*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie]) * f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.), pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie], 2.))*(*susy)("UMIX", je*10+0)*(*susy)("UMIX", ie*10+0)
							-f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie], 2.)) *  (*susy)("VMIX", je * 10+0)*(*susy)("VMIX", ie * 10+0));
						
					}
				}

				for (int be=0; be<3; be++) {
					double ratio_Mch = (*sus_param).Mch[je] / (*sus_param).Mch[ie];
					double ratio_MsqU = (*sus_param).MsqU[ae] / (*sus_param).Mch[ie];
					double ratio_Msn = (*sus_param).Msn[be] / (*sus_param).Mch[ie];
					
					B1c1 += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] / pow((*sus_param).Mch[ie], 2) * 
							(0.5 * (*sus_param).X_NL[ie][be][1] * (*sus_param).X_NL[je][be][1] * 
							(f81(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
							4 * (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
							(f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
							f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
							log(pow(scale / (*sus_param).MsqU[ae], 2))));

					B1c2 += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] / pow((*sus_param).Mch[ie], 2) * 
							((*sus_param).X_NR[ie][be][1] * (*sus_param).X_NR[je][be][1] * fabs(ratio_Mch) * 
							(f91(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
							4 * (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
							(f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
							f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
							log(pow(scale / (*sus_param).MsqU[ae], 2))));
				}
			}



			

			for (int ce = 0; ce < 6; ce++) {
				for (int fe = 0; fe < 3; fe++) {
					for (int de = 0; de < 6; de++) {
						for (int ke = 0; ke < 6; ke++) {
							double MsqU_ke_Mch_ie_squared = pow((*sus_param).MsqU[ke] / (*sus_param).Mch[ie], 2.0);
							double log_scale_MsqU_ke = log(pow(scale / (*sus_param).MsqU[ke], 2.0));
							C91f += (*sus_param).P_U[de][ke] * MsqU_ke_Mch_ie_squared * (*sus_param).P_U[ke][ae] * 
									(1.0 + log_scale_MsqU_ke) * (*sus_param).X_UL[ie][ce][1] * (*sus_param).X_UL[ie][ae][2] * 
									f50(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) * 
									(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[de][fe];

							C91f += (*sus_param).P_U[de][ke] * MsqU_ke_Mch_ie_squared * (*sus_param).P_U[ke][ce] * 
									(1.0 + log_scale_MsqU_ke) * (*sus_param).X_UL[ie][de][1] * (*sus_param).X_UL[ie][ae][2] * 
									f50(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0)) * 
									(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[ae][fe];
							
						}

						for (int je = 0; je < 2; je++) {
							// logger->info(std::to_string((*sus_param).X_NL[ie][fe][1]));
							double factor_common = mass24_Mch_ie_squared * (*sus_param).P_U[ae][de] *pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2)*(*sus_param).P_U[de][ce] *
							(1+log(pow(scale/(*sus_param).MsqU[de],2))) *  (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ce][2];

							B1f1 += factor_common * 0.5 * f90(pow((*sus_param).Mch[je] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
										pow((*sus_param).Msn[fe] / (*sus_param).Mch[ie], 2.0)) *
									(*sus_param).X_NL[ie][fe][1] * (*sus_param).X_NL[je][fe][1];

							B1f2 += factor_common * fabs((*sus_param).Mch[je] / (*sus_param).Mch[ie]) * 
									f100(pow((*sus_param).Mch[je] / (*sus_param).Mch[ie], 2.0), 
											pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
											pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
											pow((*sus_param).Msn[fe] / (*sus_param).Mch[ie], 2.0)) *
									(*sus_param).X_NR[ie][fe][1] * (*sus_param).X_NR[je][fe][1];

							
						}
					}

					C91c += (*sus_param).X_UL[ie][ce][1] * (*sus_param).X_UL[ie][ae][2] *
						(f51(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) +
							4.0 * (f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) + 
							(f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)* 1.0001) - 
							f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)* 0.9999)) / 0.0002
							+(f40(ratio_MsqU_ae_Mch_ie*1.0001, std::pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.))
							-f40(ratio_MsqU_ae_Mch_ie*0.9999, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)))/0.0002)* log_mu_W_MsqU_ae) *
						(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[ae][fe];
				}

				double MsqU_ce_Mch_ie_squared = pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0);
				double log_scale_MsqU_ce = 1.0 + log(pow(scale / (*sus_param).MsqU[ce], 2.0));
				

				double log_mu_W_MsqU_ce = std::log(std::pow(scale / (*sus_param).MsqU[ce], 2.0));
				double MsqU_be_Mch_ie_ratio = (*sus_param).MsqU[ce] / (*sus_param).Mch[ie];
				
				for (int de = 0; de < 6; de++) {
					double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
					double ratio_MsqU_de_Mch_ie = std::pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0);
					
					D91f += mass24_Mch_ie_squared * (*sus_param).P_U[ae][ce] * MsqU_ce_Mch_ie_squared * (*sus_param).P_U[ce][de] * log_scale_MsqU_ce * 
						(*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][de][2] * 
						q51(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), ratio_MsqU_de_Mch_ie);

					C7four_1 += std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).P_U[ae][ce] * MsqU_be_Mch_ie_ratio * (*sus_param).P_U[ce][de] * (1.0 + log_mu_W_MsqU_ce) *
								((*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][de][2] * (-q11(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) + 2.0 / 3.0 * q21(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie)) +
								(*sus_param).Mch[ie] / (*sus_param).mass_b_muW * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UR[ie][de][2] * (-q31(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) + 2.0 / 3.0 * q41(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie)));

					C8four_1 += std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).P_U[ae][ce] * MsqU_be_Mch_ie_ratio * (*sus_param).P_U[ce][de] * (1.0 + log_mu_W_MsqU_ce) *
								((*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][de][2] * q21(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) +
								(*sus_param).Mch[ie] / (*sus_param).mass_b_muW * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UR[ie][de][2] * q41(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie));
				}
			}
		}
	}


	C4charg_1*=(*sus_param).kappa;
	C7charg_1*=-0.5*(*sus_param).kappa;
	C8charg_1*=-0.5*(*sus_param).kappa;

	logger->debug("kappa in NLO SUSY " + std::to_string((*sus_param).kappa));
	logger->info("C4charg_1  in NLO SUSY " + std::to_string(std::real(C4charg_1)));
	logger->debug("C7charg_1  in NLO SUSY " + std::to_string(std::real(C7charg_1)));
	logger->info("C8charg_1  in NLO SUSY " + std::to_string(std::real(C8charg_1)));

	C91c *= -(*sus_param).kappa / 8.0;
	D91c *= (*sus_param).kappa;
	complex_t B91c = -(B1c1 - B1c2) * (*sus_param).kappa * (*sm)("MASS", 24) * (*sm)("MASS", 24) / 2.0 / pow((*sm)("GAUGE", 2), 2);
	complex_t B101c = (B1c1 + B1c2) * (*sus_param).kappa * (*sm)("MASS", 24) * (*sm)("MASS", 24) / 2.0 / pow((*sm)("GAUGE", 2), 2);

	C91f *= (*sus_param).kappa / 6.0;
	// B1f1 *= (*sus_param).kappa;
	// B1f2 *= (*sus_param).kappa;
	D91f *= (*sus_param).kappa;

	complex_t B91f = (B1f1 - B1f2) * 2.0 / 3.0 * (*sus_param).kappa / pow((*sm)("GAUGE", 2), 2);
	complex_t B101f = -(B1f1 + B1f2) * 2.0 / 3.0 * (*sus_param).kappa / pow((*sm)("GAUGE", 2), 2);

	logger->info("D91f  in NLO SUSY " + doubleToString(std::real(D91f), 20));
	logger->info("C91f  in NLO SUSY " + doubleToString(std::real(C91f), 20));
	logger->info("B1f1  in NLO SUSY " + doubleToString(std::real(B1f1), 20));
	logger->info("B1f2  in NLO SUSY " + doubleToString(std::real(B1f2), 20));
	logger->info("B101f  in NLO SUSY " + doubleToString(std::real(B101f), 20));
	logger->info("B91f  in NLO SUSY " + doubleToString(std::real(B91f), 20));

	logger->info("B1c1  in NLO SUSY " + doubleToString(std::real(B1c1), 20));
	logger->info("B1c2  in NLO SUSY " + doubleToString(std::real(B1c2), 20));
	logger->info("B91c  in NLO SUSY " + doubleToString(std::real(B91c), 20));
	logger->info("B101c  in NLO SUSY " + doubleToString(std::real(B101c), 20));
	logger->info("C91c  in NLO SUSY " + doubleToString(std::real(C91c), 20));
	logger->info("D91c  in NLO SUSY " + doubleToString(std::real(D91c), 20));
    complex_t C9four_1 = (1. - 4. * (*sus_param).sw2) / (*sus_param).sw2 * C91f - B91f / (*sus_param).sw2 - D91f;
    complex_t C10four_1 = (B101f - C91f) / (*sus_param).sw2;	


	complex_t C9charg_1=(1.-4.*(*sus_param).sw2)/(*sus_param).sw2*C91c-B91c/(*sus_param).sw2-D91c;
	complex_t C10charg_1=(B101c-C91c)/(*sus_param).sw2;

	Logger::getInstance()->info("C10Charg_1 : " + std::to_string(std::real(C10charg_1)));
	Logger::getInstance()->info("C10four_1 : " + std::to_string(std::real(C10four_1)));

	Logger::getInstance()->info("C9Charg_1 : " + std::to_string(std::real(C9charg_1)));
	Logger::getInstance()->info("C9four_1 : " + std::to_string(std::real(C9four_1)));
	double alphas_mu = sm->QCDRunner.runningAlphasCalculation(scale);
	

	SUSY_LO_Strategy::init(scale, C_match);

	std::unique_ptr<THDM_NLO_Strategy> thdm_nlo = std::make_unique<THDM_NLO_Strategy>();

	thdm_nlo->set_lu(1/(*susy)("HMIX", 2)); 
	thdm_nlo->set_ld(-(*susy)("HMIX", 2));
	thdm_nlo->init(scale, C_match);

	if (C_match.size() < 2) C_match.resize(2);
    auto& C_LO = C_match[0];
    auto& C_NLO = C_match[1];
	
    if (C_NLO.empty()) C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));

	auto adjustCoefficient = [&](std::complex<double>& Cx_NLO, int index) {
        double ratio = alphas_mu / (4.0 * Pi);
        double absCx_NLO = std::abs(Cx_NLO) * ratio;
        if (absCx_NLO > std::abs(C_LO[index])) {
            Cx_NLO *= std::abs(C_LO[index]) / std::abs(Cx_NLO) * (1.0 / ratio);
        }
    };

	// adjustCoefficient(C7charg_1, 7);
	// adjustCoefficient(C7four_1, 7);

	// adjustCoefficient(C8charg_1, 8);
	// adjustCoefficient(C8four_1, 8);

	// adjustCoefficient(C9charg_1, 9);
	// adjustCoefficient(C9four_1, 9);

	// adjustCoefficient(C10charg_1, 10);
	// adjustCoefficient(C10four_1, 10);

	logger->info("C8H_1 " + std::to_string(std::real(C_NLO[7])));
	logger->info("C4H_1 " + std::to_string(std::real(C_NLO[3])));

	logger->info("C8Char_1 " + std::to_string(std::real(C7charg_1)));
	logger->info("C8four_1 " + std::to_string(std::real(C7four_1)));
	logger->info("C8H_1 " + std::to_string(std::real(C_NLO[7])));

	C_NLO[static_cast<size_t>(WilsonCoefficient::C1)] += 0.;
    C_NLO[static_cast<size_t>(WilsonCoefficient::C4)] += C4charg_1;
    C_NLO[static_cast<size_t>(WilsonCoefficient::C7)] += C7charg_1 + C7four_1;
    C_NLO[static_cast<size_t>(WilsonCoefficient::C8)] += C8charg_1 + C8four_1;
    C_NLO[static_cast<size_t>(WilsonCoefficient::C9)] += C9charg_1 + C9four_1;
    C_NLO[static_cast<size_t>(WilsonCoefficient::C10)] += C10charg_1 + C10four_1;

	logger->info("SUSY NLO Wilson Coefficient Initialized at scale " +std::to_string(scale)+" terminated successfully");
}


void SUSY_NNLO_Strategy::init(double scale, WilsonSet& C_match) {

	auto* epsi = EpsilonCalculator::GetInstance();
	auto* susy = Parameters::GetInstance(1);
	auto* sm = Parameters::GetInstance(0);
	auto* sus_param = susy_parameters::GetInstance(scale);

	complex_t C3charg_2 = 0.0;
	complex_t C4charg_2 = 0.0;
	complex_t C4four_2 = 0.0;
	complex_t C4charg_1 = 0.0;

	for(int ie = 0; ie < 2; ie++) {
		for(int ae = 0; ae < 6; ae++) {
			double ratio_mass_W_Mch = std::pow((*sm)("MASS",24)/ (*sus_param).Mch[ie], 2.0);
			double ratio_MsqU_Mch = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
			double log_mu_W_MsqU = std::log(std::pow(scale / (*sus_param).MsqU[ae], 2.0));

			C3charg_2 += ratio_mass_W_Mch * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ae][2] * h71(ratio_MsqU_Mch, log_mu_W_MsqU);
			C4charg_2 += ratio_mass_W_Mch * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ae][2] * h41(ratio_MsqU_Mch, log_mu_W_MsqU);
			C4charg_1+= pow((*sm)("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
			for(int be = 0; be < 6; be++) {
				for(int ce = 0; ce < 6; ce++) {
					C4four_2 += ratio_mass_W_Mch * (*sus_param).P_U[ae][be] * (*sus_param).MsqU[be] / (*sus_param).Mch[ie] * (*sus_param).P_U[be][ce] *
								(1.0 + log_mu_W_MsqU) * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ce][3] *
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
	if (std::all_of(begin((*sus_param).MsqU), end((*sus_param).MsqU), [&](double m) { return std::abs(m) > (*sm)("MASS", 24) / 2.0; })) {
		C1squark_2 = -208.0 / 3.0;
		for (int ae = 0; ae < 6; ++ae) {
			double xsqa = std::pow((*sus_param).MsqU[ae] / (*sm)("MASS", 24), 2.0);
			if (4.0 * xsqa > 1.0) {
				double angle = 2.0 * asin(0.5 / sqrt(xsqa));
				C1squark_2 += -2.0 * std::pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
			}
			C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;

			xsqa = std::pow((*sus_param).MsqD[ae] / (*sm)("MASS", 24), 2.0);
			if (4.0 * xsqa > 1.0) {
				double angle = 2.0 * asin(0.5 / sqrt(xsqa));
				C1squark_2 += -2.0 * std::pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
			}
			C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;
		}
	}
	
	
	std::unique_ptr<THDM_NNLO_Strategy> thdm_nnlo = std::make_unique<THDM_NNLO_Strategy>();

	thdm_nnlo->set_lu(1/(*susy)("HMIX", 2));
	Logger::getInstance()->info(std::to_string((*susy)("HMIX", 2)));
	thdm_nnlo->set_ld(-(*susy)("HMIX", 2));
	Logger::getInstance()->info("TOUT VA BIEEEEN0");
	thdm_nnlo->init(scale, C_match);

	if (C_match.size() < 3) C_match.resize(3);
    auto& C_NNLO = C_match[2];
    if (C_NNLO.empty()) C_NNLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, std::complex<double>(0, 0));

	C_NNLO[static_cast<size_t>(WilsonCoefficient::C1)] += C1squark_2;
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C3)] += C3charg_2;
    C_NNLO[static_cast<size_t>(WilsonCoefficient::C4)] += C4charg_2+C4four_2;
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C5)] += C5charg_2;
	C_NNLO[static_cast<size_t>(WilsonCoefficient::C6)] += C6charg_2;

	Logger::getInstance()->info("SUSY NNLO Wilson Coefficient Initialized at scale " +std::to_string(scale)+" terminated successfully");
}


void SUSY_LO_Strategy::init_prime(double Q_match,double Q,int gen, WilsonSet& C) {
	
	Parameters* sm = Parameters::GetInstance();
	Parameters* susy = Parameters::GetInstance(1);
	EpsilonCalculator* epsi = EpsilonCalculator::GetInstance();
	susy_parameters* sus_param = susy_parameters::GetInstance(Q_match);

	double ml;

	
	if(gen==1) ml=(*sm)("MASS", 11);
	else if(gen==3) ml=(*sm)("MASS", 13);
	else {gen=2; ml=(*sm)("MASS", 15);}

	double alphas_muW = (*sm).alpha_s(Q_match);
    double alphas_mu = (*sm).alpha_s(Q);
    double eta_mu = alphas_muW / alphas_mu;

    double mass_c_muW = (*sm).running_mass((*sm)("MASS", 4), (*sm)("MASS", 4), Q_match, "running");

	double epsfac=pow((1.+(*epsi).epsilon_b()*(*susy)("HMIX",2)),2.);

	double C7pH=(*sm)("MASS",3)*(*sus_param).mass_b_muW/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW*1./3.*(*sus_param).ld*(*sus_param).ld*F7_1((*sus_param).yt);
	double C8pH=(*sm)("MASS",3)*(*sus_param).mass_b_muW/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW*1./3.*(*sus_param).ld*(*sus_param).ld*F8_1((*sus_param).yt);

	double C7pcharg=0.;
	double C8pcharg=0.;

	double B10pc=0.;
	double C9pc=0.;
	double B9pc=0.;
	double D9pc=0.; 

	double BQ1pc1=0.;
	double BQ1pc2=0.;

	double Dp,Dm;
	double a0a, a0b, a0c, a0Q1, a0Q2, a1, NQ1pc, NQ2pc;

	for(int ie=0;ie<2;ie++) {
		for(int ae=1;ae<=6;ae++) {
			C7pcharg+=pow((*sm)("MASS",24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h10(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)) + (*sus_param).Mch[ie]/(*sus_param).mass_b_muW*(*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h20(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
			C8pcharg+=pow((*sm)("MASS",24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h50(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)) + (*sus_param).Mch[ie]/(*sus_param).mass_b_muW*(*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h60(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));

			for(int je=0;je<2;je++) {
				C9pc+=(*sus_param).X_UR[je][ae][2]*(*sus_param).X_UR[ie][ae][2]*(2.*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f30(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*(*susy)("VMIX", je*10+1)*(*susy)("VMIX", ie*10+1) -f40(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*(*susy)("UMIX", je*10+1)*(*susy)("UMIX", ie*10+1));
				B10pc*=(*sus_param).kappa* ((*sm)("MASS",24))*(*sm)("MASS",24)/2./(*sm)("MASS",2)/(*sm)("MASS",2);

				for(int be=0;be<3;be++) {
					B10pc+=-(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*(0.5*(*sus_param).X_NR[ie][be][2]*(*sus_param).X_NR[je][be][2]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)) +(*sus_param).X_NL[ie][be][2]*(*sus_param).X_NL[je][be][2]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
					B9pc+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*(0.5*(*sus_param).X_NR[ie][be][2]*(*sus_param).X_NR[je][be][2]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)) -(*sus_param).X_NL[ie][be][2]*(*sus_param).X_NL[je][be][2]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
					B9pc*=(*sus_param).kappa*((*sm)("MASS",24))*(*sm)("MASS",24)/2./(*sm)("GAUGE", 2)/(*sm)("GAUGE", 2);
					BQ1pc1+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UL[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][2]*(*sus_param).X_NR[je][be][2]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))); 	
					BQ1pc2+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UL[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][2]*(*sus_param).X_NL[je][be][2]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
					for(int me=0;me<3;me++) {
						for(int ne=0;ne<3;ne++) {
							Dp=0.;
							Dm=0.;
							for(int fe=0;fe<3;fe++) { 
								Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[be][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[be][fe]);
								Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[be][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[be][fe]);
							}
							a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,be);

							a0b=-(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", je*10+1)*(*susy)("VMIX", ie*10+0))*kron(ae,be);
							a0c=1./(*sm)("MASS",24)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.))*kron(ie,je);
							a0Q1=a0a+a0b+Dp*a0c;
							a0Q2=-a0a+a0b+Dm*a0c;
							a1=(*sus_param).Mch[ie]/sqrt(2.)/(*sm)("MASS",24)*f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,be);
							
							NQ1pc+=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[be][me]*(*susy)("UMIX", je*10+1)*(a0Q1+a1*(*susy)("HMIX",2));
							NQ2pc+=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[be][me]*(*susy)("UMIX", je*10+1)*(a0Q2+a1*(*susy)("HMIX",2));
						}
					}
				}

			}
			for(int be=0;be<6;be++) {
				for(int ce=0;ce<3;ce++) {
					C9pc+=-(*sus_param).X_UR[ie][be][1]*(*sus_param).X_UR[ie][ae][2]*f40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.))*(*sus_param).Gamma_UR[be][ce]*(*sus_param).Gamma_UR[ae][ce];
					}
			}
			D9pc+=pow((*sm)("MASS",24)/(*sus_param).Mch[ie],2.)*(*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.));
		}
	} 		
	C7pcharg*=-0.5*(*sus_param).kappa; 
		
	
	C8pcharg*=-0.5*(*sus_param).kappa; 
	C9pc*=-(*sus_param).kappa/8.;

	if (C.size() < 1) C.resize(1); 
    auto& C_LO = C[0]; 
    C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));


	C_LO[static_cast<size_t>(WilsonCoefficient::CP7)]+=pow(eta_mu,16./23.)*(C7pH+C7pcharg);
	C_LO[static_cast<size_t>(WilsonCoefficient::CP8)]+=pow(eta_mu,14./23.)*(C8pH+C8pcharg);
	
	/* Wilson coefficients C9 and C10 prime */ 	
	double C10pH = -(*sus_param).mass_b_muW*((*sm)("MASS",3))*((*susy)("HMIX",2)*(*susy)("HMIX",2)/8./(*sm)("MASS",24)/(*sm)("MASS",24)
	+pow(ml*(*susy)("HMIX",2)*(*susy)("HMIX",2)/4./(*sm)("MASS",24)/(*susy)("MASS",37),2.))*f20((*sus_param).yt)/(*sus_param).sw2;
	
	double C9pH =(4.*(*sus_param).sw2-1.)*C10pH - (*sm)("MASS",3)*(*sus_param).mass_b_muW/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW*D9H0((*sus_param).yt,(*sus_param).ld);
	
	double C10pcharg=(B10pc-C9pc)/(*sus_param).sw2;
	
	
	
	
	
	D9pc*=(*sus_param).kappa;

	double C9pcharg=(1.-4.*(*sus_param).sw2)/(*sus_param).sw2*C9pc-B9pc/(*sus_param).sw2-D9pc;
	
	C_LO[static_cast<size_t>(WilsonCoefficient::CP9)]+=C9pH+C9pcharg;
	C_LO[static_cast<size_t>(WilsonCoefficient::CP10)]+=C10pH+C10pcharg;	

	/* Wilson coefficients CQ1 and CQ2 prime */ 
	double NQ1pH=-ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./(*sm)("MASS",24)/(*sm)("MASS",24)*(*sus_param).xt*f30((*sus_param).xt,(*sus_param).z);
	
	double BQ1pH=ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./(*sm)("MASS",24)/(*sm)("MASS",24)*f70((*sus_param).xt,(*sus_param).z);
	
	complex_t CQ1pH=(NQ1pH+BQ1pH)*(*sm)("MASS",3)/(*sus_param).sw2;
	
	complex_t CQ2pH=CQ1pH;
	
	
	
	
	double BQ1pc=(BQ1pc1+BQ1pc2)*(*sus_param).kappa*((*sm)("MASS",24))*(*sm)("MASS",24)/2./(*sm)("GAUGE", 2)/(*sm)("GAUGE", 2)/(*sus_param).sw2;
	double BQ2pc=(BQ1pc1-BQ1pc2)*(*sus_param).kappa*((*sm)("MASS",24))*(*sm)("MASS",24)/2./(*sm)("GAUGE", 2)/(*sm)("GAUGE", 2)/(*sus_param).sw2;

	NQ1pc*=ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/(*sm)("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-(*sm)("MASS",24)*(*sm)("MASS",24))*(*sus_param).aY*((*sm)("MASS",3))/(*sus_param).sw2;
	NQ2pc*=ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/(*sm)("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-(*sm)("MASS",24)*(*sm)("MASS",24))*(*sus_param).aY*((*sm)("MASS",3))/(*sus_param).sw2;
	
	complex_t CQ1pcharg=NQ1pc+BQ1pc;
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]=CQ1pH+CQ1pcharg;
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]/=epsfac;
	

	complex_t CQ2pcharg=NQ2pc+BQ2pc;
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]=CQ2pH+CQ2pcharg;
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]/=epsfac;

	int nf=5;
	double beta0 = 11.-2./3.*nf;
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]*=pow(eta_mu,-4./beta0);
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]*=pow(eta_mu,-4./beta0);

	Logger::getInstance()->info("SUSY LO Wilson Primes Coefficient Initialized from scale " +std::to_string(Q_match)+" to scale" + std::to_string(Q) + " terminated successfully");
}

void SUSY_LO_Strategy::init_scalar(double Q_match,double Q,int gen, WilsonSet& C) {

	Parameters* sm = Parameters::GetInstance();
	Parameters* susy = Parameters::GetInstance(1);
	EpsilonCalculator* epsi = EpsilonCalculator::GetInstance();
	susy_parameters* sus_param = susy_parameters::GetInstance(Q_match);

	double ml;
	if(gen==1) ml=(*sm)("MASS", 11);
	else if(gen==3) ml=(*sm)("MASS", 13);
	else {gen=2; ml=(*sm)("MASS", 15);}

	complex_t BQ10c1=0.;
	complex_t BQ10c2=0.;

	double Dp, Dm;
	double a0a, a0b, a0c, a0Q1, a0Q2;
	double a1;

	complex_t NQ10c=0.;
	complex_t NQ20c=0.;

	for(int ie=0;ie<2;ie++) {
		for(int je=0;je<2;je++) {
			for(int ae=0;ae<6;ae++) {
				for(int be=0;be<3;be++) { 
					BQ10c1+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][2]*(*sus_param).X_NL[je][be][2]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
					BQ10c2+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][2]*(*sus_param).X_NR[je][be][2]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));

					for(int me=0;me<6;me++) {
						for(int ne=0;ne<3;ne++) {
							Dp=0.;
							Dm=0.;
							for(int fe=0;fe<3;fe++) 
							{
								Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
								Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
							}
							a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,me);
							a0b=-(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", je*10+1)*(*susy)("VMIX", ie*10+0))*kron(ae,me);
							a0c=1./(*sm)("MASS",24)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))*kron(ie,je);
							a0Q1=a0a+a0b+Dp*a0c;
							a0Q2=-a0a+a0b+Dm*a0c;
							
							a1=(*sus_param).Mch[ie]/sqrt(2.)/(*sm)("MASS",24)*f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me);
							NQ10c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q1+a1*(*susy)("HMIX",2));
							NQ20c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q2+a1*(*susy)("HMIX",2));

						}
					
					}
				}
			}
		}
	}
	complex_t BQ10c=(BQ10c1+BQ10c2)*(*sus_param).kappa*(*sm)("MASS",24)*(*sm)("MASS",24)/2./(*sm)("GAUGE",2)/(*sm)("GAUGE",2)/(*sus_param).sw2;
	complex_t BQ20c=-(BQ10c1-BQ10c2)*(*sus_param).kappa*(*sm)("MASS",24)*(*sm)("MASS",24)/2./(*sm)("GAUGE",2)/(*sm)("GAUGE",2)/(*sus_param).sw2;
	
	
	

	NQ10c*=ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/(*sm)("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-(*sm)("MASS",24)*(*sm)("MASS",24))*(*sus_param).aY*(*sus_param).mass_b_muW/(*sus_param).sw2;
	NQ20c*=-ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/(*sm)("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-(*sm)("MASS",24)*(*sm)("MASS",24))*(*sus_param).aY*(*sus_param).mass_b_muW/(*sus_param).sw2;


	complex_t CQ1charg_0=NQ10c+BQ10c;
	complex_t CQ2charg_0=NQ20c+BQ20c;
	double epsfac=pow((1.+(*epsi).epsilon_b()*(*susy)("HMIX",2)),2.);

	if (C.size() < 1) C.resize(1); 
    auto& C_LO = C[0]; 
    C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));


	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]=CQ1charg_0;
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]/= epsfac;
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]=CQ2charg_0;
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]/=epsfac;

	
	Logger::getInstance()->info("SUSY LO Wilson Scalar Coefficient Initialized from scale " +std::to_string(Q_match)+" to scale" + std::to_string(Q) + " terminated successfully");
}


void SUSY_NLO_Strategy::init_scalar(double Q_match,double Q,int gen, WilsonSet& C) {

	SUSY_LO_Strategy::init_scalar(Q_match, Q, gen, C);

	Parameters* sm = Parameters::GetInstance();
	Parameters* susy = Parameters::GetInstance(1);
	EpsilonCalculator* epsi = EpsilonCalculator::GetInstance();
	susy_parameters* sus_param = susy_parameters::GetInstance(Q_match);

	double alphas_mu = sm->QCDRunner.runningAlphasCalculation(Q);
	double eta_mu=(*sus_param).alphas_muW/alphas_mu;
	double ml;
	if(gen==1) ml=(*sm)("MASS", 11);
	else if(gen==3) ml=(*sm)("MASS", 13);
	else {gen=2; ml=(*sm)("MASS", 15);}

	int nf=5;
	double beta0 = 11.-2./3.*nf;
	
	

	/* NLO - Charged Higgs */
	
	complex_t NQ11H=-ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./(*sm)("MASS",24)/(*sm)("MASS",24)*(f141((*sus_param).xt,(*sus_param).z)+8.*(*sus_param).xt*(f30((*sus_param).xt,(*sus_param).z)+(*sus_param).xt*(f30((*sus_param).xt*1.0001,(*sus_param).z)-f30((*sus_param).xt*0.9999,(*sus_param).z))/0.0002)*log(Q_match*Q_match/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW));
	complex_t BQ11H=ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./(*sm)("MASS",24)/(*sm)("MASS",24)*(f111((*sus_param).xt,(*sus_param).z)+8.*(f70((*sus_param).xt*1.0001,(*sus_param).z)-f70((*sus_param).xt*0.9999,(*sus_param).z))/0.0002*log(Q_match*Q_match/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW));
	complex_t CQ1H_1=(NQ11H+BQ11H)*(*sus_param).mass_b_muW/(*sus_param).sw2;
	complex_t CQ2H_1=-CQ1H_1;

	/* NLO - charginos */
	
	complex_t BQ11c1=0.;
	complex_t BQ11c2=0.;
	complex_t NQ11c=0.;
	complex_t NQ21c=0.;
	complex_t BQ11f1=0.;
	complex_t BQ11f2=0.;
	complex_t NQ11f=0.;
	complex_t NQ21f=0.;

	double Dp, Dm, temp, temp2;
	double a0a, a0b, a0c, a0Q1, a0Q2, a0p, a1,a2p;

	for(int ie=0;ie<2;ie++) {
		for(int ae=0;ae<6;ae++){
			for(int je=0;je<2;je++)  {
				for(int be=0;be<3;be++) {
					BQ11c1+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][1]*(*sus_param).X_NL[je][be][1]*(f121(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))+4.*(f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))-f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.))));
					BQ11c2+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][1]*(*sus_param).X_NR[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*(f131(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))+4.*(f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))-f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.))));

					for(int me=1;me<=6;me++){ 
						for(int ne=1;ne<=3;ne++) {
							Dp=0.;
							Dm=0.;
							for(int fe=1;fe<=3;fe++) { 	
								Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
								Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
							}
							a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*(f181(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))+4.*(f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*1.0001)-f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,me);
							a0b=-((f191(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))+4.*(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*1.0001)-f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*(*susy)("UMIX", je*10+1)*(*susy)("UMIX", ie*10+0))*kron(ae,me);
							a0c=1./(*sm)("MASS",24)*(f171(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))+4.*(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))+(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))-f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)))/0.0002+(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)*1.0001)-f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)*0.9999))/0.0002)*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*kron(ie,je);
						
							a0Q1=a0a+a0b+Dp*a0c;
							a0Q2=-a0a+a0b+Dm*a0c;
							a0p=4.*(*sus_param).G_aimn[ae][ie][be][ne]/(*sm)("MASS",24)/((*sm)("VCKM", be*10+2)*(*sm)("VCKM", ne*10+1)/(*sm)("VCKM", 22)/(*sm)("VCKM", 21))/(*susy)("UMIX", je*10+1)*f151(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
							a1=(*sus_param).Mch[ie]/sqrt(2.)/(*sm)("MASS",24)*(f161(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))+4.*(f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001)-f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*kron(ie,je)*kron(ae,me);
							a2p=(*sus_param).Gamma_UL[me][be]*((*sm)("VCKM", be*10+2)*(*sm)("VCKM", ne*10+1)/(*sm)("VCKM", 22)/(*sm)("VCKM", 21))*(*susy)("UMIX", je*10+1)/2./(*sm)("MASS",24)*f151(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
							
							NQ11c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q1+a1*(*susy)("HMIX",2))
							+(*sus_param).G_aimn[ae][ie][be][ne]*(*susy)("UMIX", je*10+1)*a0p
							+(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*a2p*pow((*sm)("MASS",3)*(*susy)("HMIX",2),2.);	
							NQ21c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q2+a1*(*susy)("HMIX",2))
							+(*sus_param).G_aimn[ae][ie][be][ne]*(*susy)("UMIX", je*10+1)*a0p
							+(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*a2p*pow((*sm)("MASS",3)*(*susy)("HMIX",2),2.);
						}
					}

				}
				for(int be=1;be<=6;be++) {
					for(int ce=1;ce<=6;ce++) {
						for(int fe=1;fe<=3;fe++) {
							BQ11f1+=-(*sus_param).X_UL[je][be][1]*(*sus_param).X_UR[ie][ae][2]*pow((*sm)("MASS",24)/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][ce]*(*sus_param).MsqU[ce]/(*sus_param).Mch[ie]*(*sus_param).P_U[ce][be]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ce],2.)))	*(f90(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[fe]/(*sus_param).Mch[ie],2.))*(*sus_param).X_NR[ie][fe][2]*(*sus_param).X_NL[je][fe][2]);
							BQ11f2+=-(*sus_param).X_UL[je][be][1]*(*sus_param).X_UR[ie][ae][2]*pow((*sm)("MASS",24)/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][ce]*(*sus_param).MsqU[ce]/(*sus_param).Mch[ie]*(*sus_param).P_U[ce][be]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ce],2.)))	*(fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f100(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[fe]/(*sus_param).Mch[ie],2.))*(*sus_param).X_NL[ie][fe][1]*(*sus_param).X_NR[je][fe][1]);

						}
					}
				}
			}

			for(int me=1;me<=3;me++) {
				for(int ne=1;ne<=3;ne++) {
					for(int de=1;de<=6;de++) {
						for(int ke=1;ke<=6;ke++) {
							
							temp2 = (*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[de][me]*(*susy)("UMIX",ie*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[ie]*(*sus_param).P_U[ke][de]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*(*susy)("HMIX",2)*(*sus_param).Mch[ie]/sqrt(2.)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2.));
							NQ11f+=temp2;
							NQ21f+=temp2;
							for(int ce=1;ce<=6;ce++) {
								Dp=0.;
								Dm=0.;
								for(int fe=1;fe<=3;fe++) 
								{		
									Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[de][fe]*(*sus_param).Gamma_UL[ce][fe]+(*sus_param).Gamma_UL[de][fe]*(*sus_param).Gamma_UR[ce][fe]); 
									Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[de][fe]*(*sus_param).Gamma_UL[ce][fe]-(*sus_param).Gamma_UL[de][fe]*(*sus_param).Gamma_UR[ce][fe]); 
								}
								temp=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[ce][me]*(*susy)("UMIX",ie*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[ie]*(*sus_param).P_U[ke][de]*
								(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*f60(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.));	
								NQ11f+=Dp*temp;
								NQ21f+=Dm*temp;
							}
							for(int je=1;je<=2;je++) {
								temp=-(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[de][me]*(*susy)("UMIX",je*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[je]*(*sus_param).P_U[ke][de]*
								(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*(*sm)("MASS",24)*(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f60(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[je],2.))*(*susy)("UMIX",ie*10+2)*(*susy)("VMIX",je*10+1)+
								f50(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[je],2.))*(*susy)("UMIX",je*10+2)*(*susy)("VMIX",ie*10+1)); 
								NQ11f+=temp;
								NQ21f+=-temp;
							}

						}
					}
				}
					
			}
		}
	}
	//Warning, triple calculation in superiso, not done here
	
	complex_t BQ11c=(BQ11c1+BQ11c2)*(*sus_param).kappa*(*sm)("MASS",24)*(*sm)("MASS",24)/2./(*sm)("GAUGE",2)/(*sm)("GAUGE",2)/(*sus_param).sw2;
	complex_t BQ21c=-(BQ11c1-BQ11c2)*(*sus_param).kappa*(*sm)("MASS",24)*(*sm)("MASS",24)/2./(*sm)("GAUGE",2)/(*sm)("GAUGE",2)/(*sus_param).sw2;
	
	
	NQ11c*=ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/(*sm)("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-(*sm)("MASS",24)*(*sm)("MASS",24))*(*sus_param).aY*(*sus_param).mass_b_muW/(*sus_param).sw2;
	NQ21c*=-ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/(*sm)("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-(*sm)("MASS",24)*(*sm)("MASS",24))*(*sus_param).aY*(*sus_param).mass_b_muW/(*sus_param).sw2;
	
	complex_t CQ1charg_1=NQ11c+BQ11c;
	
	if (C.size() < 1) C.resize(2); 
    auto& C_LO = C[0];
	auto& C_NLO = C[1]; 
    C_LO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));
	C_NLO.resize(static_cast<size_t>(WilsonCoefficient::CPQ2) + 1, complex_t(0, 0));


	if(fabs(CQ1charg_1)*alphas_mu/4./Pi>fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)])) CQ1charg_1*=fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)])/fabs(CQ1charg_1)*4.*Pi/alphas_mu;
	if(fabs(CQ1H_1)*alphas_mu/4./Pi>fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)])) CQ1H_1*=fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)])/fabs(CQ1H_1)*4.*Pi/alphas_mu;
	C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ1)]=CQ1H_1+CQ1charg_1;
	
	C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ1)]/=(*sus_param).epsfac;
	
	complex_t CQ2charg_1=NQ21c+BQ21c;
	
	if(fabs(CQ2charg_1)*alphas_mu/4./Pi>fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)])) CQ2charg_1*=fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)])/fabs(CQ2charg_1)*4.*Pi/alphas_mu;
	if(fabs(CQ2H_1)*alphas_mu/4./Pi>fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)])) CQ2H_1*=fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)])/fabs(CQ2H_1)*4.*Pi/alphas_mu;
	C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ2)]=CQ2H_1+CQ2charg_1;
		
	C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ2)]/=(*sus_param).epsfac;
	
	
	/* Wilson coefficient CQ1 */ 
	/* NLO  - four points */
		
	complex_t BQ11f=(BQ11f1+BQ11f2)*2./3.*(*sus_param).kappa/(*sm)("GAUGE",2)/(*sm)("GAUGE",2)/(*sus_param).sw2;
	complex_t BQ21f=-(BQ11f1-BQ11f2)*2./3.*(*sus_param).kappa/(*sm)("GAUGE",2)/(*sm)("GAUGE",2)/(*sus_param).sw2;
	
	
	NQ11f*=-4./3.*ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/(*sm)("MASS",24)/(*sm)("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-(*sm)("MASS",24)*(*sm)("MASS",24))*(*sus_param).aY*(*sus_param).mass_b_muW/(*sus_param).sw2;

	NQ21f*=4./3.*ml*((*susy)("HMIX",2))*(*susy)("HMIX",2)/(*sm)("MASS",24)/(*sm)("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-(*sm)("MASS",24)*(*sm)("MASS",24))*(*sus_param).aY*(*sus_param).mass_b_muW/(*sus_param).sw2;


	complex_t CQ1four_1=NQ11f+BQ11f;

	if(fabs(CQ1four_1)*alphas_mu/4./Pi>fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)])) CQ1four_1*=fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)])/fabs(CQ1four_1)*4.*Pi/alphas_mu;
		
	C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ1)]+=CQ1four_1;

	complex_t CQ2four_1=NQ21f+BQ21f;
	
	if(fabs(CQ2four_1)*alphas_mu/4./Pi>fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)])) CQ2four_1*=fabs(C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)])/fabs(CQ2four_1)*4.*Pi/alphas_mu;
	C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ2)]+=CQ2four_1;

	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]*=pow(eta_mu,-4./beta0);
	C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]*=pow(eta_mu,-4./beta0);
	C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ1)]*=pow(eta_mu,-4./beta0)*eta_mu;
	C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ2)]*=pow(eta_mu,-4./beta0)*eta_mu;


	/* NMSSM */

	double lambdaNMSSM = 1;
	double lambdaSNMSSM = 1;
	double AlambdaNSSM = 1;
	double kappaNMSSM = 1;
	double m_Bs = 1;
	double mass_nutl = 1;

	if((*susy)("MASS",36)!=0.)
	{
		double s=lambdaSNMSSM/lambdaNMSSM;
		double v=sqrt(1./sqrt(2.)/(*sm)("SMINPUTS", 2));
		
		double v_deltam_s=v/s*(sqrt(2.)*AlambdaNSSM-2.*kappaNMSSM*s)/(sqrt(2.)*AlambdaNSSM+kappaNMSSM*s);
		
		C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]=0.;
		C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]=0.;
		C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ1)]=0.;
		C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ2)]=0.;
		
		double mH0[4],mA0[3],mstop[3];
		
		mstop[0]=(*susy)("MASS", 2000013); //mass upr, is that right ?
		mstop[1]=(*susy)("MASS", 1000006);
		mstop[2]=(*susy)("MASS", 2000006);
		
		mH0[1]=(*susy)("MASS", 25);
		mH0[2]=(*susy)("MASS", 35);
		mH0[3]=(*susy)("MASS", 36);
		mA0[1]=(*susy)("MASS",36);
		mA0[2]=(*susy)("MASS",36);
		
		double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
		double T1[3][4][4],T2[4][4][4];
		std::array<std::array<double,4>,4> TU;
	
		TU[1][1]=1.;
		for(int ie=1;ie<=2;ie++){
			 for(int je=1;je<=2;je++) {
				TU[ie+1][je+1]=(*susy)("STOPMIX", ie*10+je);
			}
		}

		
		double vu=sqrt(pow(sin(atan((*susy)("HMIX",2))),2.)/sqrt(2.)/(*sm)("SMINPUTS", 2));
		double vd=vu/(*susy)("HMIX",2);

		for(int je=1;je<=2;je++) {
			for(int le=1;le<=2;le++) {
				 for(int ae=1;ae<=3;ae++) {
					if (ae <3 ){
						Ralj[ae][le][je]=-(*sm)("GAUGE",2)/sqrt(2.)*((*susy)("A0MIX",ae*10+1)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je)+(*susy)("A0MIX",ae*10+2)*(*susy)("UMIX",10+le)*(*susy)("VMIX",20+je))-lambdaNMSSM/sqrt(2.)*(*susy)("A0MIX",ae*10+3)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je);
					}
					Qalj[ae][le][je]=(*sm)("GAUGE",2)/sqrt(2.)*((*susy)("H0MIX",ae*10+1)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je)+(*susy)("H0MIX",ae*10+2)*(*susy)("UMIX",10+le)*(*susy)("VMIX",20+je))-lambdaNMSSM/sqrt(2.)*(*susy)("H0MIX",ae*10+3)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je);
					for(int ke=1;ke<=3;ke++) {
						G1[ae][ke][je][le]=(TU[ae][2]*TU[ke][2]-kron(ae,1)*kron(ke,1))*(*susy)("VMIX",10+le)*(*susy)("UMIX",20+je)-(*sus_param).mass_top_muW/sqrt(2.)/sin(atan((*susy)("HMIX",2)))/(*sm)("MASS",24)*TU[ae][3]*TU[ke][2]*(*susy)("VMIX",20+le)*(*susy)("UMIX",20+je);
					}
				}
			}
			for(int ie=1;ie<=3;ie++) {
				for(int ke=1;ke<=3;ke++) {
					T1[je][ie][ke]=(TU[ie][3]*TU[ke][2]-TU[ie][2]*TU[ke][3])*((lambdaNMSSM/sqrt(2.)*(vd*(*susy)("A0MIX",je*10+3)+s*(*susy)("A0MIX",je*10+1)))-(*susy)("AU", 11)*(*susy)("A0MIX",je*10+2));
				}
			}
		}

		complex_t CQ1H=0.;
		complex_t CQ2H=0.;
		complex_t CQ1c=0.;
		complex_t CQ2c=0.;
		complex_t CAc=0.;

		for(int ae=1;ae<=3;ae++) {
			for(int ie=1;ie<=3;ie++) {
				for(int ke=1;ke<=3;ke++){
					T2[ae][ie][ke]=-(*sus_param).mass_top_muW/2./(*sm)("MASS",24)*(2.*(*sus_param).mass_top_muW*(*susy)("H0MIX",ae*10+2)*(TU[ie][2]*TU[ke][2]+TU[ie][3]*TU[ke][3])	+((lambdaNMSSM/sqrt(2.)*(vd*(*susy)("H0MIX",ae*10+3)+s*(*susy)("H0MIX",ae*10+1)))+(*susy)("AU", 11)*(*susy)("H0MIX",ae*10+2))*(TU[ie][3]*TU[ke][2]+TU[ie][2]*TU[ke][3]))	+(*sm)("MASS",23)/2./sqrt(1.-(*sus_param).sw2)*(1.-4./3.*(*sus_param).sw2)*(*susy)("H0MIX",ae*10+2)*(TU[ie][1]*TU[ke][1]+TU[ie][2]*TU[ke][2])+2./3.*(*sm)("MASS",24)*(*sus_param).sw2/(1.-(*sus_param).sw2)*(*susy)("H0MIX",ae*10+2)*TU[ie][3]*TU[ke][3];

					for(int je=1;je<=2;je++) {
						for(int le=1;le<=2;le++) {
							CQ1c+=G1[ie][ke][je][le]/mH0[ae]/mH0[ae]*( sqrt(2.)*(*susy)("H0MIX",ae*10+1)*(*susy)("H0MIX",ae*10+1)*(*sus_param).Mch[je]/(*sm)("MASS",24)/cos(atan((*susy)("HMIX",2)))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/(*sus_param).Mch[je],2.))
							-2.*sqrt(2.)*(*susy)("H0MIX",ae*10+1)/(*sm)("GAUGE",2)*kron(ie,ke)*(Qalj[ae][le][je]*f40(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))+(*sus_param).Mch[je]/(*sus_param).Mch[le]*Qalj[ae][je][le]*f30(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))		+2.*sqrt(2.)*(*susy)("H0MIX",ae*10+1)*T2[ae][ie][ke]*(*sus_param).Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow((*sus_param).Mch[je]/mstop[ke-1],2.))
							+mH0[ae]*mH0[ae]/(*sus_param).Mch[je]/(*sus_param).Mch[je]*kron(ie,ke)*((*susy)("UMIX",20+je)*(*susy)("VMIX",10+le)*f50(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))
							-(*sus_param).Mch[le]/(*sus_param).Mch[je]*(*susy)("UMIX",20+le)*(*susy)("VMIX",10+je)* f60(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))));

							CQ2c+=G1[ie][ke][je][le]/mA0[ae]/mA0[ae]*(sqrt(2.)*(*susy)("A0MIX",ae*10+1)*(*susy)("A0MIX",ae*10+1)*(*sus_param).Mch[je]/(*sm)("MASS",24)/cos(atan((*susy)("HMIX",2)))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/(*sus_param).Mch[je],2.))
							-2.*sqrt(2.)*(*susy)("A0MIX",ae*10+1)/(*sm)("GAUGE",2)*kron(ie,ke)*(-Ralj[ae][le][je]*f40(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))+(*sus_param).Mch[je]/(*sus_param).Mch[le]*Ralj[ae][je][le]*f30(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))			-sqrt(2.)*(*susy)("A0MIX",ae*10+1)*T1[ae][ie][ke]*(*sus_param).mass_top_muW*(*sus_param).Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow((*sus_param).Mch[je]/mstop[ke-1],2.))
							+mA0[ae]*mA0[ae]/(*sus_param).Mch[je]/(*sus_param).Mch[je]*kron(ie,ke)*((*susy)("UMIX",20+je)*(*susy)("VMIX",10+le)*f50(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))
							-(*sus_param).Mch[le]/(*sus_param).Mch[je]*(*susy)("UMIX",20+le)*(*susy)("VMIX",10+je)*f60(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))));
						}		
					}

				}
			}
			CQ1H+=((*susy)("MASS",37)*(*susy)("MASS",37)/(*sm)("MASS",24)/(*sm)("MASS",24)*(*susy)("H0MIX",ae*10+1)*(*susy)("H0MIX",ae*10+1)*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,(*sm)("MASS",24)*(*sm)("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)	+(*sus_param).mass_top_muW*(*sus_param).mass_top_muW*mH0[ae]*mH0[ae]/(*sm)("MASS",24)/(*sm)("MASS",24)/(*susy)("MASS",37)/(*susy)("MASS",37)*f30((*sus_param).mass_top_muW*(*sus_param).mass_top_muW/(*susy)("MASS",37)/(*susy)("MASS",37),(*sus_param).mass_top_muW*(*sus_param).mass_top_muW/(*sm)("MASS",24)/(*sm)("MASS",24)))/mH0[ae]/mH0[ae];
			if (ae < 3) {
				CQ2H+=(((*susy)("MASS",37)*(*susy)("MASS",37)/(*sm)("MASS",24)/(*sm)("MASS",24)*(*susy)("A0MIX",ae*10+1)*(*susy)("A0MIX",ae*10+1)+kron(ae,2)*(*susy)("A0MIX",ae*10+1))*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,(*sm)("MASS",24)*(*sm)("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)	+(*sus_param).mass_top_muW*(*sus_param).mass_top_muW*mA0[ae]*mA0[ae]/(*sm)("MASS",24)/(*sm)("MASS",24)/(*susy)("MASS",37)/(*susy)("MASS",37)*f30((*sus_param).mass_top_muW*(*sus_param).mass_top_muW/(*susy)("MASS",37)/(*susy)("MASS",37),(*sus_param).mass_top_muW*(*sus_param).mass_top_muW/(*sm)("MASS",24)/(*sm)("MASS",24)))/mA0[ae]/mA0[ae];
			}
			for(int je=1;je<=2;je++) {
				for(int le=1;le<=2;le++) {
					CAc = std::complex<double>(CAc.real(), CAc.imag()+((*susy)("HMIX",2))/sqrt(2.)*G1[ae][ae][je][le]*(v_deltam_s*kron(le,je)*fabs((*sus_param).Mch[je]/(*sm)("MASS",24))*f80(pow(mstop[ae-1]/(*sus_param).Mch[je],2.))-(Ralj[1][je][le]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[le])*f30(pow(mstop[ae-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))-Ralj[1][le][je]*f40(pow(mstop[ae-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))));
					}
				}
		}
	
		CQ1H*=-ml/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
		CQ2H*=ml/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
		
		complex_t CAH={0,-lambdaNMSSM*AlambdaNSSM/(*sm)("GAUGE",2)/(*sm)("MASS",24)*(*susy)("HMIX",2)*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,(*sm)("MASS",24)*(*sm)("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)};
	
			
		CQ1c*=ml/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
		CQ2c*=-ml/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);		
	
		C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]=(CQ1H+CQ1c)*(*sus_param).mass_b_muW/(*sus_param).sw2/(*sus_param).epsfac;
		C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]=(CQ2H+CQ2c)*(*sus_param).mass_b_muW/(*sus_param).sw2/(*sus_param).epsfac;

		complex_t CA=CAH+CAc;

		if((*susy)("MASS",36)>Q_match) C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]+=-v_deltam_s/2.*(*sus_param).mass_b_muW/(*sus_param).sw2*ml*CA/(*susy)("MASS",36)/(*susy)("MASS",36);

		C_LO[static_cast<size_t>(WilsonCoefficient::CPQ1)]*=pow(eta_mu,-4./beta0);
		C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]*=pow(eta_mu,-4./beta0);
		C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ1)]*=pow(eta_mu,-4./beta0)*eta_mu;
		C_NLO[static_cast<size_t>(WilsonCoefficient::CPQ2)]*=pow(eta_mu,-4./beta0)*eta_mu;
		
		if(((*susy)("MASS",36)>(*sm).QCDRunner.get_mb_pole())&&((*susy)("MASS",36)<Q_match ))
		{	
			double alphas_Ma1=(*sm).QCDRunner.runningAlphasCalculation((*susy)("MASS",36));	
			double eta_a1=alphas_Ma1/alphas_mu;
			double mass_b_ma1=(*sm).running_mass((*sm)("MASS", 5),(*sm)("MASS", 5),(*susy)("MASS",36));
			C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]+=-v_deltam_s/2.*mass_b_ma1/(*sus_param).sw2*ml*CA/(*susy)("MASS",36)/(*susy)("MASS",36)*pow(eta_a1,-4./beta0);
		}
		
		if((*susy)("MASS",36)<(*sm).QCDRunner.get_mb_pole())
		{	
			double width_A0=1.e-6;
			C_LO[static_cast<size_t>(WilsonCoefficient::CPQ2)]+=std::complex<double>{v_deltam_s/2.*(*sm)("MASS", 5)/(*sus_param).sw2*ml*CA/(m_Bs*m_Bs-(*susy)("MASS",36)*(*susy)("MASS",36),(*susy)("MASS",36)*width_A0)};
		}

	}
	Logger::getInstance()->info("SUSY NLO Wilson Scalar Coefficient Initialized from scale " +std::to_string(Q_match)+" to scale" + std::to_string(Q) + " terminated successfully");
}