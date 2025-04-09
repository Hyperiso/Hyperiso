#include "susy_parameters.h"

// CHOOSE REAL PART OF CKM WHY?????
// susy_parameters::susy_parameters(double scale) {

//     mass_H03 = 0.;
// 	mass_A02 = 0.; // for testing
	
// 	LOG_DEBUG("c11 " + std::to_string(std::real(c11)));
// 	LOG_DEBUG("c12 " + std::to_string(std::real(c12)));
// 	LOG_DEBUG("c13 " + std::to_string(std::real(c13)));
// 	LOG_DEBUG("c21 " + std::to_string(std::real(c21)));
// 	LOG_DEBUG("c22 " + std::to_string(std::real(c22)));
// 	LOG_DEBUG("c23 " + std::to_string(std::real(c23)));
// 	LOG_DEBUG("c31 " + std::to_string(std::real(c31)));
// 	LOG_DEBUG("c32 " + std::to_string(std::real(c32)));
// 	LOG_DEBUG("c33 " + std::to_string(std::real(c33)));

// 	LOG_DEBUG("imag c11 " + doubleToString(std::imag(c11), 20));
// 	LOG_DEBUG("imag c12 " + doubleToString(std::imag(c12), 20));
// 	LOG_DEBUG("imag c13 " + doubleToString(std::imag(c13), 20));
// 	LOG_DEBUG("imag c21 " + doubleToString(std::imag(c21), 20));
// 	LOG_DEBUG("imag c22 " + doubleToString(std::imag(c22), 20));
// 	LOG_DEBUG("imag c23 " + doubleToString(std::imag(c23), 20));
// 	LOG_DEBUG("imag c31 " + doubleToString(std::imag(c31), 20));
// 	LOG_DEBUG("imag c32 " + doubleToString(std::imag(c32), 20));
// 	LOG_DEBUG("imag c33 " + doubleToString(std::imag(c33), 20));
	
// 	// epsilonbp=(*epsi).epsilon_bp();
// 	// epsilon0p=(*epsi).epsilon_0p();
// 	// epsilon0=(*epsi).epsilon_0();
// 	// epsilon2=(*epsi).epsilon_2();
// 	// epsilon1p=(*epsi).epsilon_1p();
// 	// epsilonb=epsilon0+epsilon2;

// 	// LOG_DEBUG("epsilon0 " + std::to_string(std::real(epsilon0)));
// 	// LOG_DEBUG("epsilon2 " + std::to_string(std::real(epsilon2)));

//     mass_top_muW = QCDHelper::msbar_mass(6, scale, MassType::MSBAR, MassType::POLE);
// 	mass_b_muW = QCDHelper::msbar_mass(5, scale, MassType::MSBAR, MassType::POLE); //mass bottom 6 (at pole)

// 	L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
//  	sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.); //1 = param-> gp and 2 = (*sm)("GAUGE",2)

// 	LOG_DEBUG("sw2 : " + std::to_string(sw2));

// 	xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
// 	yt= pow(mass_top_muW/(*susy)("MASS",37),2.); // param->mass_H (37)

//     lu=1./(*susy)("HMIX",2);
// 	ld=-(*susy)("HMIX",2);
	
// 	alphas_muW = QCDHelper::alpha_s(scale);

//     alphas_mg = QCDHelper::alpha_s((*susy)("MASS",1000021));
// 	ag = 1.0 - 7.0 / (12.0 * Pi) * alphas_mg;
// 	aY = 1.0 + alphas_mg / (4.0 * Pi);
// 	kappa = 1.0 / ((*sm)("GAUGE",2) * (*sm)("GAUGE",2) * std::real(VCKM[2][2]*VCKM[2][1])); //VCKM 33 et 32

// 	LOG_DEBUG("prod : " + doubleToString((*sm)("GAUGE",2) * (*sm)("GAUGE",2) * std::real(VCKM[2][2]*VCKM[2][1]), 20));
// 	LOG_INFO("vckm 22 : " + doubleToString(std::real(VCKM[2][2]), 20));
// 	LOG_INFO("vckm 21 : " + doubleToString(std::real(VCKM[2][1]), 20));
// 	LOG_INFO("g2 : " + doubleToString((*sm)("GAUGE",2),20 ));
// 	z=pow((*susy)("MASS",37)/(*sm)("MASS",24),2.);
// 	sinb = std::sin(std::atan((*susy)("HMIX",2)));
// 	cosb = std::cos(std::atan((*susy)("HMIX",2)));
// 	ct = (*susy)("STOPMIX",11); // Ajustement des indices pour base-0
// 	st = (*susy)("STOPMIX",01); // Ajustement des indices pour base-0
// 	LOG_DEBUG("ST " + std::to_string((*susy)("STOPMIX",01)));

//     // Initialisation des masses
// 	MU = {(*sm)("MASS",2), (*sm)("MASS",4), mass_top_muW}; 
// 	MD = {(*sm)("MASS",1), (*sm)("MASS",3), mass_b_muW}; 
// 	ME = {(*sm)("MASS",11), (*sm)("MASS",13), (*sm)("MASS",15)}; 

// 	Mch = {(*susy)("MASS",1000024), (*susy)("MASS",1000037) };
// 	// Array1D_7 MsqU = { 0.0, param->mass_upl, param->mass_chl, param->mass_t1, param->mass_upr, param->mass_chr, param->mass_t2}; // Ajout d'un élément pour compatibilité de taille
// 	MsqU = {(*susy)("MASS",1000002), (*susy)("MASS",1000004), (*susy)("MASS",1000006), (*susy)("MASS",2000002), (*susy)("MASS",2000004), (*susy)("MASS",2000006)}; // Ajout d'un élément pour compatibilité de taille
	
// 	MsqD = {
//     (*susy)("MASS",1000001),
//     (*susy)("MASS",1000003),
//     (*susy)("MASS",1000005),
//     (*susy)("MASS",2000001),
//     (*susy)("MASS",2000003),
//     (*susy)("MASS",2000005)};
// 	Msn = {(*susy)("MASS",1000012), (*susy)("MASS",1000014), (*susy)("MASS",1000016)};
// 	LOG_INFO("Msn : " + std::to_string(Msn[0]) + " " + std::to_string(Msn[1]) + " " + std::to_string(Msn[2]));
// 	// Vérification du mélange sU_mix et initialisation conditionnelle de Gamma_UL et Gamma_UR
// 	bool isNonZeroMix = true;
// 	for (size_t i = 0; i < NumSquarks; ++i) {
// 		double product = 1.0;
// 		for (size_t j = 0; j < 6; ++j) { // Pour chaque colonne dans la rangée i
// 			product *= sU_mix[i][j];
// 		}
// 		if (product == 0.0) {
// 			isNonZeroMix = false;
// 			break;
// 		}
// 	}
// 	LOG_INFO("STILL FINE");
// 	if (isNonZeroMix) {
// 		// Tri de MsqU si la condition est vraie
// 		std::sort(MsqU.begin(), MsqU.end());

// 		// Remplissage des matrices Gamma_UL et Gamma_UR
// 		for (size_t ae = 0; ae < NumSquarks; ++ae) {
// 			for (size_t ie = 0; ie < 3; ++ie) { // Indices ajustés pour base-0
// 				Gamma_UL[ae][ie] = sU_mix[ae][ie];
// 				Gamma_UR[ae][ie] = sU_mix[ae][ie + 3];
// 			}
// 		}
// }
// 	else {
// 		// Configuration spécifique basée sur les exigences décrites
// 		Gamma_UL[0][0] = 1.0; // Correspond à Gamma_UL[1][1] = 1 dans l'indexation originale de base-1
// 		Gamma_UL[1][1] = 1.0; // Correspond à Gamma_UL[2][2] = 1
// 		Gamma_UL[2][2] = ct;  // ct déjà calculé précédemment
// 		Gamma_UL[5][2] = -st; // -st, ajustement pour base-0

// 		Gamma_UR[3][0] = 1.0; // Correspond à Gamma_UR[4][1] = 1 dans l'indexation originale de base-1
// 		Gamma_UR[4][1] = 1.0; // Correspond à Gamma_UR[5][2] = 1
// 		Gamma_UR[2][2] = st;  // st
// 		Gamma_UR[5][2] = ct;  // ct
// 	}

// 	for (int ae = 0; ae < 6; ++ae) {
//         for (int ie = 0; ie < 3; ++ie) {
//             Gamma_U[ae][ie] = Gamma_UL[ae][ie];
//             Gamma_U[ae][ie+3] = Gamma_UR[ae][ie];
// 			if (ae <3 && ae==ie) {
// 				Gamma_NL[ae][ie] = 1.;
// 			}
//         }
//     }

//     I_LR.fill({});
//     for (int i = 0; i < 3; ++i) {
//         I_LR[i][i] = 1.;
//         I_LR[i+3][i+3] = -1.;
//     }
//     // Calcul de P_U
//     for (int ae = 0; ae < 6; ++ae) {
//         for (int be = 0; be < 6; ++be) {
//             for (int ce = 0; ce < 6; ++ce) {
//                 for (int de = 0; de < 6; ++de) {
//                     P_U[ae][be] = Gamma_U[ae][ce] * I_LR[ce][de] * Gamma_U[be][de];
// 					LOG_DEBUG("P_U[" + std::to_string(ae) + "][" + std::to_string(be) + "] = " + std::to_string(P_U[ae][be]));
//                 }
//             }
//         }
//     }

	
// 	LOG_DEBUG("CT " + std::to_string(isNonZeroMix));
// 	LOG_DEBUG("ST " + std::to_string(st));

//     // Calculs pour X_UL et X_UR
//     for (int ie = 0; ie < 2; ++ie) {
// 		for (int ae = 0; ae < 6; ++ae) { // Supposant que 6 est correct pour tous
// 			for (int be = 0; be < 3; ++be) {
// 				// Réinitialisation pour X_UL et X_UR
// 				X_UL[ie][ae][be] = 0.0;
// 				X_UR[ie][ae][be] = 0.0;

// 				// Calculs pour X_UL et X_UR
// 				for (int ce = 0; ce < 3; ++ce) {
// 					X_UL[ie][ae][be] += -(*sm)("GAUGE",2) * (
// 						ag * (*susy)("VMIX", ie*10+0) * Gamma_UL[ae][ce] -
// 						aY * (*susy)("VMIX", ie*10+1) * Gamma_UR[ae][ce] * MU[ce] / (sqrt(2.0) * (*sm)("MASS", 24) * sinb)
// 					) * std::real(VCKM[ce][be]);
// 					// LOG_INFO("GAMMA_UL " + std::to_string(ae) + " " + std::to_string(ce)  + " " + std::to_string(Gamma_UL[ae][ce]));
// 					X_UR[ie][ae][be] += (*sm)("GAUGE",2) * aY * (*susy)(std::string("UMIX"), ie*10+1) * Gamma_UL[ae][ce] * std::real(VCKM[ce][be]) * MD[be] / (sqrt(2.0) * (*sm)("MASS", 24) * cosb);

// 					G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*(*sm)("MASS",24)*(*susy)("VMIX", ie*10+0)*Gamma_UL[ae][ce]*ag-MU[ce]*(*susy)("VMIX", ie*10+1)*Gamma_UR[ae][ce]*aY)*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])*std::real(VCKM[2][1]));
// 					// LOG_INFO(std::to_string(std::real(VCKM[ce][be])) + " WAOUW");
// 				}
// 				LOG_DEBUG("X_UL[" + std::to_string(ie) + "][" +std::to_string(ae)+"]["+std::to_string(be)+"] = "+std::to_string(X_UL[ie][ae][be]));
// 				LOG_DEBUG("X_UR[" + std::to_string(ie) + "][" +std::to_string(ae)+"]["+std::to_string(be)+"] = "+std::to_string(X_UR[ie][ae][be]));
// 				// Condition pour éviter le dépassement dans X_NL et X_NR si ae > 2
// 				if (ae < 3) {
// 					X_NL[ie][ae][be] = -(*sm)("GAUGE",2) * (*susy)("VMIX", ie*10+0) * Gamma_NL[ae][be];
// 					X_NR[ie][ae][be] = (*sm)("GAUGE",2) * (*susy)("UMIX", ie*10+1) * Gamma_NL[ae][be] * ME[be] / (sqrt(2.0) * (*sm)("MASS", 24) * cosb);
// 				}
// 			}
// 		}
// 	}

// 	LOG_INFO("AG " + std::to_string(ag));
// 	LOG_INFO("MU ", MU[0], MU[1], MU[2]);
// 	LOG_INFO("AY " + std::to_string(aY));
// 	LOG_DEBUG("vmix00 " + std::to_string((*susy)("VMIX", 0)));
// 	LOG_DEBUG("vmix01 " + std::to_string((*susy)("VMIX", 01)));
// 	LOG_DEBUG("vmix10 " + std::to_string((*susy)("VMIX", 10)));
// 	LOG_DEBUG("vmix11 " + std::to_string((*susy)("VMIX", 11)));
// 	//X_UL and X_UR change from here, -1 from superiso

// 	auto computeContributions = [&](int ie, auto func, double additionalFactor = 1.0) {
// 		double result = 0.0;
// 		for (int ae = 0; ae < 6; ++ae) {
// 			double msqOverMchSquared = std::pow(MsqU[ae] / Mch[ie], 2.0);
// 			result += (X_UL[ie][ae][0] * X_UL[ie][ae][1] * func(msqOverMchSquared) +
// 					Mch[ie] / mass_b_muW * X_UL[ie][ae][0] * X_UR[ie][ae][1] * func(msqOverMchSquared)) * additionalFactor;
// 		}
// 		return result;
// 	};

// 	LOG_DEBUG("KAPPA IS " + std::to_string(kappa));

// 	auto hFunc10 = [](double x) { return h10(x); };
// 	auto hFunc20 = [](double x) { return h20(x); };
// 	auto hFunc50 = [](double x) { return h50(x); };
// 	auto hFunc60 = [](double x) { return h60(x); };

// 	kappaFactor = -0.5 * kappa;


// 	for (int ie = 0; ie < 2; ++ie) {
// 		for (int je = 0; je < 2; ++je) {
// 			for (int ae = 0; ae < 6; ++ae) {
// 				double mchRatioSquared = std::pow(Mch[je] / Mch[ie], 2.0);
// 				double msqOverMchSquared = std::pow(MsqU[ae] / Mch[ie], 2.0);

// 				for (int be = 0; be < 3; ++be) {
// 					double msnOverMchSquared = std::pow(Msn[be] / Mch[ie], 2.0);
// 					B0c1 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (Mch[ie] * Mch[ie]) * (0.5 * X_NL[ie][be][1] * X_NL[je][be][1] * f50(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));
// 					B0c2 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (Mch[ie] * Mch[ie]) * (X_NR[ie][be][1] * X_NR[je][be][1] * std::fabs(Mch[je] / Mch[ie]) * f60(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));	
// 				}

// 				C90c += X_UL[je][ae][1] * X_UL[ie][ae][2] * (2.0 * std::fabs(Mch[je] / Mch[ie]) * f30(mchRatioSquared, msqOverMchSquared) * (*susy)("UMIX", je*10+0) * (*susy)("UMIX", ie*10+0) - f40(mchRatioSquared, msqOverMchSquared) * (*susy)("VMIX", je*10+0) * (*susy)("VMIX", ie*10+0));

// 				if (ie == je)	{
// 					D90c += std::pow((*sm)("MASS", 24) / Mch[ie], 2.0) * X_UL[ie][ae][1] * X_UL[ie][ae][2] * h30(msqOverMchSquared);
// 				}
// 			}
// 		}
// 	}
// 	LOG_DEBUG("B0c1 : " + std::to_string(B0c1*1e10));
// 	LOG_DEBUG("B0c2 : " + doubleToString(B0c2*1e16, 20));
// 	for (int ie = 0; ie < 2; ++ie) {
// 		for (int ae = 0; ae < 6; ++ae) {
// 			for (int be = 0; be < 6; ++be) {
// 				double msqOverMchSquaredAe = std::pow(MsqU[ae] / Mch[ie], 2.0);
// 				double msqOverMchSquaredBe = std::pow(MsqU[be] / Mch[ie], 2.0);
// 				for (int ce = 0; ce < 3; ++ce) {
// 					C90c += X_UL[ie][be][1] * X_UL[ie][ae][2] * f40(msqOverMchSquaredAe, msqOverMchSquaredBe) * Gamma_UL[be][ce] * Gamma_UL[ae][ce];
// 				}
// 			}
// 		}
// 	}
// 	m_H = (*susy)("MASS", 37);
// 	// Final adjustments
// 	B90c = -(B0c1 - B0c2) * kappa * std::pow((*sm)("MASS", 24), 2.0) / (2.0 * std::pow((*sm)("GAUGE",2), 2.0));
// 	B100c = (B0c1 + B0c2) * kappa * std::pow((*sm)("MASS", 24), 2.0) / (2.0 * std::pow((*sm)("GAUGE",2), 2.0));
// 	C90c *= -kappa / 8.0;
// 	D90c *= kappa;

//     test = true;
// 	for (int ae = 0; ae < 6; ++ae) { // Conserve l'indexation à partir de 1
// 		if (!(std::fabs(MsqU[ae]) > (*sm)("MASS", 24) / 2. && std::fabs(MsqD[ae]) > (*sm)("MASS", 24) / 2.)) {
// 			test = false;
// 			break; // Sort de la boucle dès qu'une condition n'est pas remplie
// 		}
// 	}
	
// }

//TODO TODO TODO
// void susy_parameters::reset_PrimeCQG(double Q_match) {
// 	if (is_PrimeCQG) {return;}
// 	is_PrimeCQG = true;

// 	std::shared_ptr<Parameters> sm = Parameters::GetInstance();
// 	double mass_c_muW = QCDHelper::msbar_mass(4, Q_match);

// 	LOG_INFO("mass_c_muW", mass_c_muW);
// 	MU = {(*sm)("MASS",2), mass_c_muW, mass_top_muW};
// 	LOG_INFO("MU_0", MU[0]);
// 	LOG_INFO("MU_1", MU[1]);
// 	LOG_INFO("MU_2", MU[2]);
// 	for (int ie = 0; ie < 2; ++ie) {
// 		for (int ae = 0; ae < 6; ++ae) { 
// 			for (int be = 0; be < 3; ++be) {
// 				X_UL[ie][ae][be] = 0.0;

// 				// Calculs pour X_UL et X_UR
// 				for (int ce = 0; ce < 3; ++ce) {
// 					X_UL[ie][ae][be] += -(*sm)("GAUGE",2) * (
// 						ag * (*susy)("VMIX", ie*10+0) * Gamma_UL[ae][ce] -
// 						aY * (*susy)("VMIX", ie*10+1) * Gamma_UR[ae][ce] * MU[ce] / (sqrt(2.0) * (*sm)("MASS", 24) * sinb)
// 					) * std::real(VCKM[ce][be]);

// 					G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*(*sm)("MASS",24)*(*susy)("VMIX", ie*10+0)*Gamma_UL[ae][ce]*ag-MU[ce]*(*susy)("VMIX", ie*10+1)*Gamma_UR[ae][ce]*aY)*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])/std::real(VCKM[2][1]));
// 				}
// 			}
// 		}
// 	}
// 	LOG_INFO("U1",(*susy)("VMIX", 1));
// 	LOG_INFO("U2",(*susy)("VMIX", 11));
// 	LOG_INFO("VCKM[2][2]", VCKM[2][2]);
// 	LOG_INFO("VCKM[2][2]", VCKM[2][1]);

// 	LOG_INFO("VCKM[0][0]", VCKM[0][0]);
// 	LOG_INFO("VCKM[0][1]", VCKM[0][1]);
// 	LOG_INFO("VCKM[0][2]", VCKM[0][2]);

// 	LOG_INFO("VCKM[1][0]", VCKM[1][0]);
// 	LOG_INFO("VCKM[1][1]", VCKM[1][1]);
// 	LOG_INFO("VCKM[1]2]", VCKM[1][2]);

// 	LOG_INFO("VCKM[2][0]", VCKM[2][0]);
// 	LOG_INFO("VCKM[2][1]", VCKM[2][1]);
// 	LOG_INFO("VCKM[2][2]", VCKM[2][2]);

// 	LOG_INFO("VCKM[ce][1]*VCKM[ce][1]", VCKM[1][1]*VCKM[0][2]);
// 	LOG_INFO("VCKM[ce][1]*VCKM[ce][1]", VCKM[1][1]*VCKM[0][2]);
// 	LOG_INFO("VCKM[be][2]", VCKM[1][2]);
// }

// void susy_parameters::reset_G() {
// 	if (!is_PrimeCQG) {return;}
// 	is_PrimeCQG=false;
// 	MU = {(*sm)("MASS",2), (*sm)("MASS",4), mass_top_muW}; 
// 	for (int ie = 0; ie < 2; ++ie) {
// 		for (int ae = 0; ae < 6; ++ae) { 
// 			for (int be = 0; be < 3; ++be) {
// 				X_UL[ie][ae][be] = 0.0;

// 				// Calculs pour X_UL et X_UR
// 				for (int ce = 0; ce < 3; ++ce) {
// 					X_UL[ie][ae][be] += -(*sm)("GAUGE",2) * (
// 						ag * (*susy)("VMIX", ie*10+0) * Gamma_UL[ae][ce] -
// 						aY * (*susy)("VMIX", ie*10+1) * Gamma_UR[ae][ce] * MU[ce] / (sqrt(2.0) * (*sm)("MASS", 24) * sinb)
// 					) * std::real(VCKM[ce][be]);

// 					G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*(*sm)("MASS",24)*(*susy)("VMIX", ie*10+0)*Gamma_UL[ae][ce]*ag-MU[ce]*(*susy)("VMIX", ie*10+1)*Gamma_UR[ae][ce]*aY)*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])*std::real(VCKM[2][1]));
// 				}

// 			}
// 		}
// 	}
// }


void susy_parameters::init() {

	if (susy_parameters::initialized) {
		return;
	}

	init_scale_independant_block();
	init_matching_block();
}

void susy_parameters::init_scale_independant_block() {

	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS"}}};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double mW = src.at("MASS")->retrieve(24)->get_val();
		double alphas_mg = QCDHelper::alpha_s(src.at("MASS")->retrieve(1000021)->get_val());
		double ag = 1.0 - 7.0 / (12.0 * Pi) * alphas_mg;
		double aY = 1.0 + alphas_mg / (4.0 * Pi);

		// TODO : Ask Nazila
		double kappa = 1.0 / (pow(src.at("GAUGE")->retrieve(2)->get_val(), 2.) * 
						std::real((src.at("RECKM")->retrieve(22)->get_val() + I * src.at("IMCKM")->retrieve(22)->get_val())*(src.at("RECKM")->retrieve(21)->get_val() + I * src.at("IMCKM")->retrieve(21)->get_val()))); //VCKM 33 et 32

		double kappaFactor = -0.5 * kappa;
		double tanb = src.at("HMIX")->retrieve(2)->get_val();
		double z = pow(src.at("MASS")->retrieve(37)->get_val() / mW, 2.);
		double sinb = std::sin(std::atan(tanb));
		double cosb = std::cos(std::atan(tanb));
		double ct = src.at("STOPMIX")->retrieve(11)->get_val();
		double st = src.at("STOPMIX")->retrieve(01)->get_val();


        double lu = 1./tanb;
        double ld = -tanb;

		//TODO : IN progress
		Array1D_4 ME = {src.at("MASS")->retrieve(11)->get_val(), src.at("MASS")->retrieve(13)->get_val(), src.at("MASS")->retrieve(15)->get_val()}; 

		Array1D_3 Mch = {src.at("MASS")->retrieve(1000024)->get_val(), src.at("MASS")->retrieve(1000037)->get_val()};

		Array1D_7 MsqU = {src.at("MASS")->retrieve(1000002)->get_val(), src.at("MASS")->retrieve(1000004)->get_val(), src.at("MASS")->retrieve(1000006)->get_val(), 
				src.at("MASS")->retrieve(2000002)->get_val(), src.at("MASS")->retrieve(2000004)->get_val(), src.at("MASS")->retrieve(2000006)->get_val()};

		Array1D_7 MsqD = {src.at("MASS")->retrieve(1000001)->get_val(), src.at("MASS")->retrieve(1000003)->get_val(), src.at("MASS")->retrieve(1000005)->get_val(), 
				src.at("MASS")->retrieve(2000001)->get_val(), src.at("MASS")->retrieve(2000003)->get_val(), src.at("MASS")->retrieve(2000005)->get_val()};

		Array1D_4 Msn = {src.at("MASS")->retrieve(1000012)->get_val(), src.at("MASS")->retrieve(1000014)->get_val(), src.at("MASS")->retrieve(1000016)->get_val()};
		const size_t NumSquarks = 6;
		Array2D_7x7 sU_mix; //ERROR
		bool isNonZeroMix = true;
		for (size_t i = 0; i < NumSquarks; ++i) {
			double product = 1.0;
			for (size_t j = 0; j < 6; ++j) {
				product *= sU_mix[i][j]; //TODO, ISSUE
			}
			if (product == 0.0) {
				isNonZeroMix = false;
				break;
			}
		}

		if (isNonZeroMix) {
			std::sort(MsqU.begin(), MsqU.end());
		}

        int id {1};
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, z, 0., 0.)); //1
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, cosb, 0., 0.)); //2
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, sinb, 0., 0.)); //3
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, ct, 0., 0.)); //4
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, st, 0., 0.)); //5
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, kappa, 0., 0.)); //6
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, lu, 0., 0.)); //7
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, ld, 0., 0.)); //8
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, alphas_mg, 0., 0.)); //9
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, ag, 0., 0.)); //10
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, aY, 0., 0.)); //11
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, ME[0], 0., 0.)); //12
		dep_block->store_or_assign({id, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 1}}, ME[1], 0., 0.)); //12
		dep_block->store_or_assign({id++, 2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 2}}, ME[2], 0., 0.)); //12
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, Mch[0], 0., 0.)); //13
		dep_block->store_or_assign({id++, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 1}}, Mch[1], 0., 0.)); //13
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqU[0], 0., 0.)); //14
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqU[1], 0., 0.)); //14
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqU[2], 0., 0.)); //14
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqU[3], 0., 0.)); //14
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqU[4], 0., 0.)); //14
		dep_block->store_or_assign({id++, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 1}}, MsqU[5], 0., 0.)); //14
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqD[0], 0., 0.)); //15
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqD[1], 0., 0.)); //15
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqD[2], 0., 0.)); //15
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqD[3], 0., 0.)); //15
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqD[4], 0., 0.)); //15
		dep_block->store_or_assign({id++, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 1}}, MsqD[5], 0., 0.)); //15
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, Msn[0], 0., 0.)); //16
		dep_block->store_or_assign({id, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 1}}, Msn[1], 0., 0.)); //16
		dep_block->store_or_assign({id++, 2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 2}}, Msn[2], 0., 0.)); //16
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, (double)isNonZeroMix, 0., 0.)); //17
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, kappaFactor, 0., 0.)); //18
    };

    susy_parameters::composer.compose_block("WPARAM_SI_BSM", src, func);



}

void susy_parameters::init_matching_block() {

	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS"}}};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double yt= pow(src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val()/src.at("MASS")->retrieve(37)->get_val(),2.); // param->mass_H (25)
		Array1D_4 MU = {src.at("MASS")->retrieve(2)->get_val(), src.at("MASS")->retrieve(4)->get_val(), src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val()}; //TODO : size 3 not 4
		Array1D_4 MD = {src.at("MASS")->retrieve(1)->get_val(), src.at("MASS")->retrieve(3)->get_val(), src.at("WPARAM_MATCH_SM")->retrieve({5,1})->get_val()}; //TODO : size 3 not 4

        dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, yt, 0., 0.));
		dep_block->store_or_assign({2,0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,0}}, MU[0], 0., 0.));
		dep_block->store_or_assign({2,1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,1}}, MU[1], 0., 0.));
		dep_block->store_or_assign({2,2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,2}}, MU[2], 0., 0.));
		dep_block->store_or_assign({3,0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {3,0}}, MD[0], 0., 0.));
		dep_block->store_or_assign({3,1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {3,1}}, MD[1], 0., 0.));
		dep_block->store_or_assign({3,2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {3,2}}, MD[2], 0., 0.));

    };

    susy_parameters::composer.compose_block("WPARAM_MATCH_BSM", src, func);


	std::unordered_map<ParameterType, std::vector<std::string>> src_matrix = {{ParameterType::SM, {"MASS"}}};

    auto func_matrix = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {

		Array2D_7x4 Gamma_UL {};
		Array2D_7x4 Gamma_UR {};
		Array2D_4x4 Gamma_NL {};
		Array2D_4x4 Gamma_NR {};
		Array2D_7x7 Gamma_U{};
		Array2D_7x7 I_LR{};
		Array2D_7x7 P_U{};
		Array3D_3x7x4 X_UL{};
		Array3D_3x7x4 X_UR{};
		Array3D_3x7x4 X_NL{};
		Array3D_3x7x4 X_NR{};
		std::array<std::array<std::array<std::array<double, 4>, 4>, 3>, 7> G_aimn;

		
		complex_t c11 = src.at("RECKM")->retrieve(00)->get_val() + src.at("IMCKM")->retrieve(00)->get_val() * complex_t(0, 1);
		complex_t c12 = src.at("RECKM")->retrieve(01)->get_val() + src.at("IMCKM")->retrieve(01)->get_val() * complex_t(0, 1);
		complex_t c13 = src.at("RECKM")->retrieve(02)->get_val() + src.at("IMCKM")->retrieve(02)->get_val() * complex_t(0, 1);
		complex_t c21 = src.at("RECKM")->retrieve(10)->get_val() + src.at("IMCKM")->retrieve(10)->get_val() * complex_t(0, 1);
		complex_t c22 = src.at("RECKM")->retrieve(11)->get_val() + src.at("IMCKM")->retrieve(11)->get_val() * complex_t(0, 1);
		complex_t c23 = src.at("RECKM")->retrieve(12)->get_val() + src.at("IMCKM")->retrieve(12)->get_val() * complex_t(0, 1);
		complex_t c31 = src.at("RECKM")->retrieve(20)->get_val() + src.at("IMCKM")->retrieve(20)->get_val() * complex_t(0, 1);
		complex_t c32 = src.at("RECKM")->retrieve(21)->get_val() + src.at("IMCKM")->retrieve(21)->get_val() * complex_t(0, 1);
		complex_t c33 = src.at("RECKM")->retrieve(22)->get_val() + src.at("IMCKM")->retrieve(22)->get_val() * complex_t(0, 1);

		
		complex_t complexTerm = -(c32 * c33 + c22 * c23) / c12;

		Array2D_4x4_I VCKM = {{
			{
				c11,
				c12,
				c13
			},
			{
				c21,
				c22,
				c23
			},
			{
				c31,
				c32,
				c33
			}
		}};

		double B0c1 = 0.0, B0c2 = 0.0, B90c = 0.0, B100c = 0.0, C90c = 0.0, D90c = 0.0;
    	bool test;

        double mW = src.at("MASS")->retrieve(24)->get_val();
		double g2 = src.at("GAUGE")->retrieve(2)->get_val();

		if (src.at("WPARAM_SI_BSM")->retrieve(17)->get_val()) {
			Array2D_7x7 sU_mix;
			const size_t NumSquarks = 6;
			for (size_t ae = 0; ae < NumSquarks; ++ae) {
				for (size_t ie = 0; ie < 3; ++ie) {
					Gamma_UL[ae][ie] = sU_mix[ae][ie];
					Gamma_UR[ae][ie] = sU_mix[ae][ie + 3];
				}
			}
	}
		else {
			Gamma_UL[0][0] = 1.0; 
			Gamma_UL[1][1] = 1.0;
			Gamma_UL[2][2] = src.at("WPARAM_SI_BSM")->retrieve(4)->get_val();
			Gamma_UL[5][2] = -src.at("WPARAM_SI_BSM")->retrieve(5)->get_val();

			Gamma_UR[3][0] = 1.0;
			Gamma_UR[4][1] = 1.0;
			Gamma_UR[2][2] = src.at("WPARAM_SI_BSM")->retrieve(5)->get_val();
			Gamma_UR[5][2] = src.at("WPARAM_SI_BSM")->retrieve(4)->get_val();
		}

		for (int ae = 0; ae < 6; ++ae) {
			for (int ie = 0; ie < 3; ++ie) {
				Gamma_U[ae][ie] = Gamma_UL[ae][ie];
				Gamma_U[ae][ie+3] = Gamma_UR[ae][ie];
				if (ae <3 && ae==ie) {
					Gamma_NL[ae][ie] = 1.;
				}
			}
		}

		I_LR.fill({});
		for (int i = 0; i < 3; ++i) {
			I_LR[i][i] = 1.;
			I_LR[i+3][i+3] = -1.;
		}

		for (int ae = 0; ae < 6; ++ae) {
			for (int be = 0; be < 6; ++be) {
				for (int ce = 0; ce < 6; ++ce) {
					for (int de = 0; de < 6; ++de) {
						P_U[ae][be] = Gamma_U[ae][ce] * I_LR[ce][de] * Gamma_U[be][de];
					}
				}
			}
		}
		
		
		for (int ie = 0; ie < 2; ++ie) {
			for (int ae = 0; ae < 6; ++ae) {
				for (int be = 0; be < 3; ++be) {
					X_UL[ie][ae][be] = 0.0;
					X_UR[ie][ae][be] = 0.0;

					for (int ce = 0; ce < 3; ++ce) {
						X_UL[ie][ae][be] += -g2 * (
							src.at("WPARAM_SI_BSM")->retrieve(10)->get_val() * src.at("VMIX")->retrieve(ie*10+0)->get_val() * Gamma_UL[ae][ce] -
							src.at("WPARAM_SI_BSM")->retrieve(11)->get_val() * src.at("VMIX")->retrieve(ie*10+1)->get_val() * Gamma_UR[ae][ce] * src.at("WPARAM_MATCH_BSM")->retrieve({2,ce})->get_val() / (sqrt(2.0) * mW * src.at("WPARAM_SI_BSM")->retrieve(3)->get_val())
						) * std::real(VCKM[ce][be]);
						X_UR[ie][ae][be] += g2 * src.at("WPARAM_SI_BSM")->retrieve(11)->get_val() * src.at("UMIX")->retrieve(ie*10+1)->get_val() * Gamma_UL[ae][ce] * std::real(VCKM[ce][be]) * src.at("WPARAM_MATCH_BSM")->retrieve({2,be})->get_val() / (sqrt(2.0) * mW * src.at("WPARAM_SI_BSM")->retrieve(2)->get_val());

						G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*mW*src.at("VMIX")->retrieve(ie*10+0)->get_val()*Gamma_UL[ae][ce]*src.at("WPARAM_SI_BSM")->retrieve(10)->get_val()-src.at("WPARAM_MATCH_BSM")->retrieve({2,ce})->get_val()*src.at("VMIX")->retrieve(ie*10+1)->get_val()*Gamma_UR[ae][ce]*src.at("WPARAM_SI_BSM")->retrieve(11)->get_val())*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])*std::real(VCKM[2][1]));
					}

					if (ae < 3) {
						X_NL[ie][ae][be] = -g2 * src.at("VMIX")->retrieve(ie*10+0)->get_val() * Gamma_NL[ae][be]; 

						X_NR[ie][ae][be] = g2 * src.at("UMIX")->retrieve(ie*10+1)->get_val() * Gamma_NL[ae][be] * src.at("WPARAM_SI_BSM")->retrieve({12, be})->get_val() / (sqrt(2.0) * mW * src.at("WPARAM_SI_BSM")->retrieve(2)->get_val()); //12 -> ME
					}
				}
			}
		}
		
		auto computeContributions = [&](int ie, auto func, double additionalFactor = 1.0) {
			double result = 0.0;
			for (int ae = 0; ae < 6; ++ae) {
				double msqOverMchSquared = std::pow(src.at("WPARAM_SI_BSM")->retrieve({14, ae})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);
				result += (X_UL[ie][ae][0] * X_UL[ie][ae][1] * func(msqOverMchSquared) + 
				src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val() / src.at("WPARAM_MATCH_SM")->retrieve({5,1})->get_val() * X_UL[ie][ae][0] * X_UR[ie][ae][1] * func(msqOverMchSquared)) * additionalFactor;
			}
			return result;
		};


		auto hFunc10 = [](double x) { return h10(x); };
		auto hFunc20 = [](double x) { return h20(x); };
		auto hFunc50 = [](double x) { return h50(x); };
		auto hFunc60 = [](double x) { return h60(x); };

		double kappaFactor = -0.5 * src.at("WPARAM_SI_BSM")->retrieve(6)->get_val();

		
		for (int ie = 0; ie < 2; ++ie) {
			for (int je = 0; je < 2; ++je) {
				for (int ae = 0; ae < 6; ++ae) {
					double mchRatioSquared = pow(src.at("WPARAM_SI_BSM")->retrieve({13, je})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0); //Mch : WPARAM_SI_BSM 13
					double msqOverMchSquared = pow(src.at("WPARAM_SI_BSM")->retrieve({14, ae})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);

					for (int be = 0; be < 3; ++be) {
						double msnOverMchSquared = pow(src.at("WPARAM_SI_BSM")->retrieve({16, be})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);
						B0c1 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val() * src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val()) * (0.5 * X_NL[ie][be][1] * X_NL[je][be][1] * f50(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));
						B0c2 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val() * src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val()) * (X_NR[ie][be][1] * X_NR[je][be][1] * std::fabs(src.at("WPARAM_SI_BSM")->retrieve({13, je})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val()) * f60(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));	
					}

					C90c += X_UL[je][ae][1] * X_UL[ie][ae][2] * (2.0 * std::fabs(src.at("WPARAM_SI_BSM")->retrieve({13, je})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val()) * f30(mchRatioSquared, msqOverMchSquared) * src.at("UMIX")->retrieve(je*10+0)->get_val() * src.at("UMIX")->retrieve(ie*10+0)->get_val() - f40(mchRatioSquared, msqOverMchSquared) * src.at("VMIX")->retrieve(je*10+0)->get_val() * src.at("VMIX")->retrieve(ie*10+0)->get_val());

					if (ie == je)	{
						D90c += pow(mW / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0) * X_UL[ie][ae][1] * X_UL[ie][ae][2] * h30(msqOverMchSquared);
					}
				}
			}
		}

		for (int ie = 0; ie < 2; ++ie) {
			for (int ae = 0; ae < 6; ++ae) {
				for (int be = 0; be < 6; ++be) {
					double msqOverMchSquaredAe = pow(src.at("WPARAM_SI_BSM")->retrieve({14, ae})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);
					double msqOverMchSquaredBe = pow(src.at("WPARAM_SI_BSM")->retrieve({14, be})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);
					for (int ce = 0; ce < 3; ++ce) {
						C90c += X_UL[ie][be][1] * X_UL[ie][ae][2] * f40(msqOverMchSquaredAe, msqOverMchSquaredBe) * Gamma_UL[be][ce] * Gamma_UL[ae][ce];
					}
				}
			}
		}
		B90c = -(B0c1 - B0c2) * src.at("WPARAM_SI_BSM")->retrieve(6)->get_val() * std::pow(mW, 2.0) / (2.0 * std::pow(g2, 2.0));
		B100c = (B0c1 + B0c2) * src.at("WPARAM_SI_BSM")->retrieve(6)->get_val() * std::pow(mW, 2.0) / (2.0 * std::pow(g2, 2.0));
		C90c *= -src.at("WPARAM_SI_BSM")->retrieve(6)->get_val() / 8.0;
		D90c *= src.at("WPARAM_SI_BSM")->retrieve(6)->get_val();

		test = true;
		for (int ae = 0; ae < 6; ++ae) {
			if (!(std::fabs(src.at("WPARAM_SI_BSM")->retrieve({14, ae})->get_val()) > mW / 2. && std::fabs(src.at("WPARAM_SI_BSM")->retrieve({15, ae})->get_val()) > mW / 2.)) {
				test = false;
				break;
			}
		}


		for (size_t i = 0; i < Gamma_UL.size(); ++i) {
			for (size_t j = 0; j < Gamma_UL[i].size(); ++j) {
				dep_block->store_or_assign({1,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {1,i,j}}, Gamma_UL[i][j], 0., 0.)); //1
				dep_block->store_or_assign({2,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {2,i,j}}, Gamma_UR[i][j], 0., 0.)); //2
			}
		}

		for (size_t i = 0; i < X_UL.size(); ++i) {
			for (size_t j = 0; j < X_UL[i].size(); ++j) {
				for (size_t k = 0; k < X_UL[i][j].size(); ++k) {
					dep_block->store_or_assign({3,i,j,k}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {3,i,j,k}}, X_UL[i][j][k], 0., 0.)); //3
					dep_block->store_or_assign({4,i,j,k}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {4,i,j,k}}, X_UR[i][j][k], 0., 0.)); //4
					dep_block->store_or_assign({5,i,j,k}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {5,i,j,k}}, X_NL[i][j][k], 0., 0.)); //5
					dep_block->store_or_assign({6,i,j,k}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {6,i,j,k}}, X_NR[i][j][k], 0., 0.)); //6
				}
			}
		}

		for (size_t i = 0; i < Gamma_U.size(); ++i) {
			for (size_t j = 0; j < Gamma_U[i].size(); ++j) {
				dep_block->store_or_assign({7,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {7,i,j}}, Gamma_U[i][j], 0., 0.)); //7
				dep_block->store_or_assign({8,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {8,i,j}}, I_LR[i][j], 0., 0.)); //8
				dep_block->store_or_assign({9,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {9,i,j}}, P_U[i][j], 0., 0.)); //9
			}
		}

		for (size_t i = 0; i < Gamma_NL.size(); ++i) {
			for (size_t j = 0; j < Gamma_NL[i].size(); ++j) {
				dep_block->store_or_assign({10,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {10,i,j}}, Gamma_NL[i][j], 0., 0.)); //10
				dep_block->store_or_assign({11,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {11,i,j}}, Gamma_NR[i][j], 0., 0.)); //11
			}
		}

		for (size_t i = 0; i < G_aimn.size(); ++i) {
			for (size_t j = 0; j < G_aimn[i].size(); ++j) {
				for (size_t k = 0; k < G_aimn[i][j].size(); ++k) {
					for (size_t l = 0; l < G_aimn[i][j][k].size(); ++l) {
						dep_block->store_or_assign({12, i, j, k, l}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {12, i, j, k, l}}, G_aimn[i][j][k][l], 0., 0.)); //12
					}
				}
			}
		}

        dep_block->store_or_assign(13, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", 13}, B90c, 0., 0.)); //13
		dep_block->store_or_assign(14, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", 14}, C90c, 0., 0.)); //14
		dep_block->store_or_assign(15, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", 15}, D90c, 0., 0.)); //15
		dep_block->store_or_assign(16, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", 16}, B100c, 0., 0.)); //16
    };

    susy_parameters::composer.compose_block("MATRIX_BSM", src_matrix, func_matrix);

}


void susy_parameters::update() {

// 	ParameterProxy wilson_p {ParameterType::WILSON};
// 	ParameterProxy susy {ParameterType::SUSY};
// 	ParameterProxy sm {ParameterType::SM};
// 	// mass_H03 = 0.;
// 	// mass_A02 = 0.; // for testing
	
// 	double mW = sm("MASS", 24);
// 	double g2 = sm("GAUGE", 2);

// 	if (wilson_p("WPARAM_SI_BSM", 17)) {

// 		const size_t NumSquarks = 6;
// 		for (size_t ae = 0; ae < NumSquarks; ++ae) {
// 			for (size_t ie = 0; ie < 3; ++ie) {
// 				w_susy.Gamma_UL[ae][ie] = sU_mix[ae][ie];
// 				w_susy.Gamma_UR[ae][ie] = sU_mix[ae][ie + 3];
// 			}
// 		}
// }
// 	else {
// 		w_susy.Gamma_UL[0][0] = 1.0; 
// 		w_susy.Gamma_UL[1][1] = 1.0;
// 		w_susy.Gamma_UL[2][2] = wilson_p("WPARAM_SI_BSM", 4);
// 		w_susy.Gamma_UL[5][2] = -wilson_p("WPARAM_SI_BSM", 5);

// 		w_susy.Gamma_UR[3][0] = 1.0;
// 		w_susy.Gamma_UR[4][1] = 1.0;
// 		w_susy.Gamma_UR[2][2] = wilson_p("WPARAM_SI_BSM", 5);
// 		w_susy.Gamma_UR[5][2] = wilson_p("WPARAM_SI_BSM", 4);
// 	}

// 	for (int ae = 0; ae < 6; ++ae) {
//         for (int ie = 0; ie < 3; ++ie) {
//             w_susy.Gamma_U[ae][ie] = w_susy.Gamma_UL[ae][ie];
//             w_susy.Gamma_U[ae][ie+3] = w_susy.Gamma_UR[ae][ie];
// 			if (ae <3 && ae==ie) {
// 				w_susy.Gamma_NL[ae][ie] = 1.;
// 			}
//         }
//     }

//     w_susy.I_LR.fill({});
//     for (int i = 0; i < 3; ++i) {
//         w_susy.I_LR[i][i] = 1.;
//         w_susy.I_LR[i+3][i+3] = -1.;
//     }

// 	for (int ae = 0; ae < 6; ++ae) {
//         for (int be = 0; be < 6; ++be) {
//             for (int ce = 0; ce < 6; ++ce) {
//                 for (int de = 0; de < 6; ++de) {
//                     w_susy.P_U[ae][be] = w_susy.Gamma_U[ae][ce] * w_susy.I_LR[ce][de] * w_susy.Gamma_U[be][de];
//                 }
//             }
//         }
//     }

	
//     for (int ie = 0; ie < 2; ++ie) {
// 		for (int ae = 0; ae < 6; ++ae) {
// 			for (int be = 0; be < 3; ++be) {
// 				w_susy.X_UL[ie][ae][be] = 0.0;
// 				w_susy.X_UR[ie][ae][be] = 0.0;

// 				for (int ce = 0; ce < 3; ++ce) {
// 					w_susy.X_UL[ie][ae][be] += -g2 * (
// 						wilson_p("WPARAM_SI_BSM", 10) * susy("VMIX", ie*10+0) * w_susy.Gamma_UL[ae][ce] -
// 						wilson_p("WPARAM_SI_BSM", 11) * susy("VMIX", ie*10+1) * w_susy.Gamma_UR[ae][ce] * wilson_p("WPARAM_MATCH_BSM", {2,ce}) / (sqrt(2.0) * mW * wilson_p("WPARAM_SI_BSM", 3))
// 					) * std::real(VCKM[ce][be]);
// 					w_susy.X_UR[ie][ae][be] += g2 * wilson_p("WPARAM_SI_BSM", 11) * susy(std::string("UMIX"), ie*10+1) * w_susy.Gamma_UL[ae][ce] * std::real(VCKM[ce][be]) * wilson_p("WPARAM_MATCH_BSM", {2,be}) / (sqrt(2.0) * mW * wilson_p("WPARAM_SI_BSM", 2));

// 					w_susy.G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*mW*susy("VMIX", ie*10+0)*w_susy.Gamma_UL[ae][ce]*wilson_p("WPARAM_SI_BSM", 10)-wilson_p("WPARAM_MATCH_BSM", {2,ce})*susy("VMIX", ie*10+1)*w_susy.Gamma_UR[ae][ce]*wilson_p("WPARAM_SI_BSM", 11))*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])*std::real(VCKM[2][1]));
// 				}

// 				if (ae < 3) {
// 					w_susy.X_NL[ie][ae][be] = -g2 * susy( "VMIX", ie*10+0) * w_susy.Gamma_NL[ae][be]; 

// 					w_susy.X_NR[ie][ae][be] = g2 * susy("UMIX", ie*10+1) * w_susy.Gamma_NL[ae][be] * wilson_p("WPARAM_SI_BSM", {12, be}) / (sqrt(2.0) * mW * wilson_p("WPARAM_SI_BSM", 2)); //12 -> ME
// 				}
// 			}
// 		}
// 	}
	
// 	auto computeContributions = [&](int ie, auto func, double additionalFactor = 1.0) {
// 		double result = 0.0;
// 		for (int ae = 0; ae < 6; ++ae) {
// 			double msqOverMchSquared = std::pow(wilson_p("WPARAM_SI_BSM", {14, ae}) / wilson_p("WPARAM_SI_BSM", {13, ie}), 2.0);
// 			result += (w_susy.X_UL[ie][ae][0] * w_susy.X_UL[ie][ae][1] * func(msqOverMchSquared) +
// 			wilson_p("WPARAM_SI_BSM", {13, ie}) / wilson_p("WPARAM_MATCH_SM", {5,1}) * w_susy.X_UL[ie][ae][0] * w_susy.X_UR[ie][ae][1] * func(msqOverMchSquared)) * additionalFactor;
// 		}
// 		return result;
// 	};


// 	auto hFunc10 = [](double x) { return h10(x); };
// 	auto hFunc20 = [](double x) { return h20(x); };
// 	auto hFunc50 = [](double x) { return h50(x); };
// 	auto hFunc60 = [](double x) { return h60(x); };

// 	kappaFactor = -0.5 * wilson_p("WPARAM_SI_BSM", 6);

	
// 	for (int ie = 0; ie < 2; ++ie) {
// 		for (int je = 0; je < 2; ++je) {
// 			for (int ae = 0; ae < 6; ++ae) {
// 				double mchRatioSquared = std::pow(wilson_p("WPARAM_SI_BSM", {13, je}) / wilson_p("WPARAM_SI_BSM", {13, ie}), 2.0); //Mch : WPARAM_SI_BSM 13
// 				double msqOverMchSquared = std::pow(wilson_p("WPARAM_SI_BSM", {14, ae}) / wilson_p("WPARAM_SI_BSM", {13, ie}), 2.0);

// 				for (int be = 0; be < 3; ++be) {
// 					double msnOverMchSquared = std::pow(wilson_p("WPARAM_SI_BSM", {16, be}) / wilson_p("WPARAM_SI_BSM", {13, ie}), 2.0);
// 					B0c1 += w_susy.X_UL[je][ae][1] * w_susy.X_UL[ie][ae][2] / (wilson_p("WPARAM_SI_BSM", {13, ie}) * wilson_p("WPARAM_SI_BSM", {13, ie})) * (0.5 * w_susy.X_NL[ie][be][1] * w_susy.X_NL[je][be][1] * f50(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));
// 					B0c2 += w_susy.X_UL[je][ae][1] * w_susy.X_UL[ie][ae][2] / (wilson_p("WPARAM_SI_BSM", {13, ie}) * wilson_p("WPARAM_SI_BSM", {13, ie})) * (w_susy.X_NR[ie][be][1] * w_susy.X_NR[je][be][1] * std::fabs(wilson_p("WPARAM_SI_BSM", {13, je}) / wilson_p("WPARAM_SI_BSM", {13, ie})) * f60(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));	
// 				}

// 				C90c += w_susy.X_UL[je][ae][1] * w_susy.X_UL[ie][ae][2] * (2.0 * std::fabs(wilson_p("WPARAM_SI_BSM", {13, je}) / wilson_p("WPARAM_SI_BSM", {13, ie})) * f30(mchRatioSquared, msqOverMchSquared) * susy("UMIX", je*10+0) * susy("UMIX", ie*10+0) - f40(mchRatioSquared, msqOverMchSquared) * susy("VMIX", je*10+0) * susy("VMIX", ie*10+0));

// 				if (ie == je)	{
// 					D90c += std::pow(mW / wilson_p("WPARAM_SI_BSM", {13, ie}), 2.0) * w_susy.X_UL[ie][ae][1] * w_susy.X_UL[ie][ae][2] * h30(msqOverMchSquared);
// 				}
// 			}
// 		}
// 	}

// 	for (int ie = 0; ie < 2; ++ie) {
// 		for (int ae = 0; ae < 6; ++ae) {
// 			for (int be = 0; be < 6; ++be) {
// 				double msqOverMchSquaredAe = std::pow(wilson_p("WPARAM_SI_BSM", {14, ae}) / wilson_p("WPARAM_SI_BSM", {13, ie}), 2.0);
// 				double msqOverMchSquaredBe = std::pow(wilson_p("WPARAM_SI_BSM", {14, be}) / wilson_p("WPARAM_SI_BSM", {13, ie}), 2.0);
// 				for (int ce = 0; ce < 3; ++ce) {
// 					C90c += w_susy.X_UL[ie][be][1] * w_susy.X_UL[ie][ae][2] * f40(msqOverMchSquaredAe, msqOverMchSquaredBe) * w_susy.Gamma_UL[be][ce] * w_susy.Gamma_UL[ae][ce];
// 				}
// 			}
// 		}
// 	}
// 	B90c = -(B0c1 - B0c2) * wilson_p("WPARAM_SI_BSM", 6) * std::pow(mW, 2.0) / (2.0 * std::pow(g2, 2.0));
// 	B100c = (B0c1 + B0c2) * wilson_p("WPARAM_SI_BSM", 6) * std::pow(mW, 2.0) / (2.0 * std::pow(g2, 2.0));
// 	C90c *= -wilson_p("WPARAM_SI_BSM", 6) / 8.0;
// 	D90c *= wilson_p("WPARAM_SI_BSM", 6);

//     test = true;
// 	for (int ae = 0; ae < 6; ++ae) {
// 		if (!(std::fabs(wilson_p("WPARAM_SI_BSM", {14, ae})) > mW / 2. && std::fabs(wilson_p("WPARAM_SI_BSM", {15, ae})) > mW / 2.)) {
// 			test = false;
// 			break;
// 		}
// 	}
}