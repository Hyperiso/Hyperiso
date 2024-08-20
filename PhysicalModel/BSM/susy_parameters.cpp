#include "susy_parameters.h"

// CHOOSE REAL PART OF CKM WHY?????
susy_parameters::susy_parameters(double scale) {

    mass_H03 = 0.;
	mass_A02 = 0.; // for testing
	
	LOG_DEBUG("c11 " + std::to_string(std::real(c11)));
	LOG_DEBUG("c12 " + std::to_string(std::real(c12)));
	LOG_DEBUG("c13 " + std::to_string(std::real(c13)));
	LOG_DEBUG("c21 " + std::to_string(std::real(c21)));
	LOG_DEBUG("c22 " + std::to_string(std::real(c22)));
	LOG_DEBUG("c23 " + std::to_string(std::real(c23)));
	LOG_DEBUG("c31 " + std::to_string(std::real(c31)));
	LOG_DEBUG("c32 " + std::to_string(std::real(c32)));
	LOG_DEBUG("c33 " + std::to_string(std::real(c33)));

	LOG_DEBUG("imag c11 " + doubleToString(std::imag(c11), 20));
	LOG_DEBUG("imag c12 " + doubleToString(std::imag(c12), 20));
	LOG_DEBUG("imag c13 " + doubleToString(std::imag(c13), 20));
	LOG_DEBUG("imag c21 " + doubleToString(std::imag(c21), 20));
	LOG_DEBUG("imag c22 " + doubleToString(std::imag(c22), 20));
	LOG_DEBUG("imag c23 " + doubleToString(std::imag(c23), 20));
	LOG_DEBUG("imag c31 " + doubleToString(std::imag(c31), 20));
	LOG_DEBUG("imag c32 " + doubleToString(std::imag(c32), 20));
	LOG_DEBUG("imag c33 " + doubleToString(std::imag(c33), 20));
	
	epsilonbp=(*epsi).epsilon_bp();
	epsilon0p=(*epsi).epsilon_0p();
	epsilon0=(*epsi).epsilon_0();
	epsilon2=(*epsi).epsilon_2();
	epsilon1p=(*epsi).epsilon_1p();
	epsilonb=epsilon0+epsilon2;

	LOG_DEBUG("epsilon0 " + std::to_string(std::real(epsilon0)));
	LOG_DEBUG("epsilon2 " + std::to_string(std::real(epsilon2)));

    mass_top_muW=(*sm).running_mass((*sm)("MASS",6), (*sm)("MASS",6),scale, "running", "pole");
	mass_b_muW=(*sm).running_mass((*sm)("MASS",5), (*sm)("MASS",5), scale, "running", "pole"); //mass bottom 6 (at pole)

	L=log(scale*scale/(*sm)("MASS",24)/(*sm)("MASS",24)); // scale -> mu_W
 	sw2=pow(sin(atan((*sm)("GAUGE",1)/(*sm)("GAUGE",2))),2.); //1 = param-> gp and 2 = (*sm)("GAUGE",2)

	LOG_DEBUG("sw2 : " + std::to_string(sw2));

	xt= pow(mass_top_muW/(*sm)("MASS",24),2.); // W boson mass (24)
	yt= pow(mass_top_muW/(*susy)("MASS",37),2.); // param->mass_H (37)

    lu=1./(*susy)("HMIX",2);
	ld=-(*susy)("HMIX",2);
	
	alphas_muW=(*sm).alpha_s(scale);

    alphas_mg = sm->alpha_s((*susy)("MASS",1000021));
	ag = 1.0 - 7.0 / (12.0 * Pi) * alphas_mg;
	aY = 1.0 + alphas_mg / (4.0 * Pi);
	kappa = 1.0 / ((*sm)("GAUGE",2) * (*sm)("GAUGE",2) * std::real(VCKM[2][2]*VCKM[2][1])); //VCKM 33 et 32

	LOG_DEBUG("prod : " + doubleToString((*sm)("GAUGE",2) * (*sm)("GAUGE",2) * std::real(VCKM[2][2]*VCKM[2][1]), 20));
	LOG_INFO("vckm 22 : " + doubleToString(std::real(VCKM[2][2]), 20));
	LOG_INFO("vckm 21 : " + doubleToString(std::real(VCKM[2][1]), 20));
	LOG_INFO("g2 : " + doubleToString((*sm)("GAUGE",2),20 ));
	z=pow((*susy)("MASS",37)/(*sm)("MASS",24),2.);
	sinb = std::sin(std::atan((*susy)("HMIX",2)));
	cosb = std::cos(std::atan((*susy)("HMIX",2)));
	ct = (*susy)("STOPMIX",11); // Ajustement des indices pour base-0
	st = (*susy)("STOPMIX",01); // Ajustement des indices pour base-0
	LOG_DEBUG("ST " + std::to_string((*susy)("STOPMIX",01)));

    // Initialisation des masses
	MU = {(*sm)("MASS",2), (*sm)("MASS",4), mass_top_muW}; 
	MD = {(*sm)("MASS",1), (*sm)("MASS",3), mass_b_muW}; 
	ME = {(*sm)("MASS",11), (*sm)("MASS",13), (*sm)("MASS",15)}; 

	Mch = {(*susy)("MASS",1000024), (*susy)("MASS",1000037) };
	// Array1D_7 MsqU = { 0.0, param->mass_upl, param->mass_chl, param->mass_t1, param->mass_upr, param->mass_chr, param->mass_t2}; // Ajout d'un élément pour compatibilité de taille
	MsqU = {(*susy)("MASS",1000002), (*susy)("MASS",1000004), (*susy)("MASS",1000006), (*susy)("MASS",2000002), (*susy)("MASS",2000004), (*susy)("MASS",2000006)}; // Ajout d'un élément pour compatibilité de taille
	
	MsqD = {
    (*susy)("MASS",1000001),
    (*susy)("MASS",1000003),
    (*susy)("MASS",1000005),
    (*susy)("MASS",2000001),
    (*susy)("MASS",2000003),
    (*susy)("MASS",2000005)};
	Msn = {(*susy)("MASS",1000012), (*susy)("MASS",1000014), (*susy)("MASS",1000016)};
	LOG_INFO("Msn : " + std::to_string(Msn[0]) + " " + std::to_string(Msn[1]) + " " + std::to_string(Msn[2]));
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
	LOG_INFO("STILL FINE");
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
    // Calcul de P_U
    for (int ae = 0; ae < 6; ++ae) {
        for (int be = 0; be < 6; ++be) {
            for (int ce = 0; ce < 6; ++ce) {
                for (int de = 0; de < 6; ++de) {
                    P_U[ae][be] = Gamma_U[ae][ce] * I_LR[ce][de] * Gamma_U[be][de];
					LOG_DEBUG("P_U[" + std::to_string(ae) + "][" + std::to_string(be) + "] = " + std::to_string(P_U[ae][be]));
                }
            }
        }
    }

	
	LOG_DEBUG("CT " + std::to_string(isNonZeroMix));
	LOG_DEBUG("ST " + std::to_string(st));

    // Calculs pour X_UL et X_UR
    for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) { // Supposant que 6 est correct pour tous
			for (int be = 0; be < 3; ++be) {
				// Réinitialisation pour X_UL et X_UR
				X_UL[ie][ae][be] = 0.0;
				X_UR[ie][ae][be] = 0.0;

				// Calculs pour X_UL et X_UR
				for (int ce = 0; ce < 3; ++ce) {
					X_UL[ie][ae][be] += -(*sm)("GAUGE",2) * (
						ag * (*susy)("VMIX", ie*10+0) * Gamma_UL[ae][ce] -
						aY * (*susy)("VMIX", ie*10+1) * Gamma_UR[ae][ce] * MU[ce] / (sqrt(2.0) * (*sm)("MASS", 24) * sinb)
					) * std::real(VCKM[ce][be]);
					// LOG_INFO("GAMMA_UL " + std::to_string(ae) + " " + std::to_string(ce)  + " " + std::to_string(Gamma_UL[ae][ce]));
					X_UR[ie][ae][be] += (*sm)("GAUGE",2) * aY * (*susy)(std::string("UMIX"), ie*10+1) * Gamma_UL[ae][ce] * std::real(VCKM[ce][be]) * MD[be] / (sqrt(2.0) * (*sm)("MASS", 24) * cosb);

					G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*(*sm)("MASS",24)*(*susy)("VMIX", ie*10+0)*Gamma_UL[ae][ce]*ag-MU[ce]*(*susy)("VMIX", ie*10+1)*Gamma_UR[ae][ce]*aY)*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])*std::real(VCKM[2][1]));
					// LOG_INFO(std::to_string(std::real(VCKM[ce][be])) + " WAOUW");
				}
				LOG_DEBUG("X_UL[" + std::to_string(ie) + "][" +std::to_string(ae)+"]["+std::to_string(be)+"] = "+std::to_string(X_UL[ie][ae][be]));
				LOG_DEBUG("X_UR[" + std::to_string(ie) + "][" +std::to_string(ae)+"]["+std::to_string(be)+"] = "+std::to_string(X_UR[ie][ae][be]));
				// Condition pour éviter le dépassement dans X_NL et X_NR si ae > 2
				if (ae < 3) {
					X_NL[ie][ae][be] = -(*sm)("GAUGE",2) * (*susy)("VMIX", ie*10+0) * Gamma_NL[ae][be];
					X_NR[ie][ae][be] = (*sm)("GAUGE",2) * (*susy)("UMIX", ie*10+1) * Gamma_NL[ae][be] * ME[be] / (sqrt(2.0) * (*sm)("MASS", 24) * cosb);
				}
			}
		}
	}

	LOG_INFO("AG " + std::to_string(ag));
	LOG_INFO("MU ", MU[0], MU[1], MU[2]);
	LOG_INFO("AY " + std::to_string(aY));
	LOG_DEBUG("vmix00 " + std::to_string((*susy)("VMIX", 0)));
	LOG_DEBUG("vmix01 " + std::to_string((*susy)("VMIX", 01)));
	LOG_DEBUG("vmix10 " + std::to_string((*susy)("VMIX", 10)));
	LOG_DEBUG("vmix11 " + std::to_string((*susy)("VMIX", 11)));
	//X_UL and X_UR change from here, -1 from superiso

	auto computeContributions = [&](int ie, auto func, double additionalFactor = 1.0) {
		double result = 0.0;
		for (int ae = 0; ae < 6; ++ae) {
			double msqOverMchSquared = std::pow(MsqU[ae] / Mch[ie], 2.0);
			result += (X_UL[ie][ae][0] * X_UL[ie][ae][1] * func(msqOverMchSquared) +
					Mch[ie] / mass_b_muW * X_UL[ie][ae][0] * X_UR[ie][ae][1] * func(msqOverMchSquared)) * additionalFactor;
		}
		return result;
	};

	LOG_DEBUG("KAPPA IS " + std::to_string(kappa));

	auto hFunc10 = [](double x) { return h10(x); };
	auto hFunc20 = [](double x) { return h20(x); };
	auto hFunc50 = [](double x) { return h50(x); };
	auto hFunc60 = [](double x) { return h60(x); };

	kappaFactor = -0.5 * kappa;


	for (int ie = 0; ie < 2; ++ie) {
		for (int je = 0; je < 2; ++je) {
			for (int ae = 0; ae < 6; ++ae) {
				double mchRatioSquared = std::pow(Mch[je] / Mch[ie], 2.0);
				double msqOverMchSquared = std::pow(MsqU[ae] / Mch[ie], 2.0);

				for (int be = 0; be < 3; ++be) {
					double msnOverMchSquared = std::pow(Msn[be] / Mch[ie], 2.0);
					B0c1 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (Mch[ie] * Mch[ie]) * (0.5 * X_NL[ie][be][1] * X_NL[je][be][1] * f50(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));
					B0c2 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (Mch[ie] * Mch[ie]) * (X_NR[ie][be][1] * X_NR[je][be][1] * std::fabs(Mch[je] / Mch[ie]) * f60(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));	
				}

				C90c += X_UL[je][ae][1] * X_UL[ie][ae][2] * (2.0 * std::fabs(Mch[je] / Mch[ie]) * f30(mchRatioSquared, msqOverMchSquared) * (*susy)("UMIX", je*10+0) * (*susy)("UMIX", ie*10+0) - f40(mchRatioSquared, msqOverMchSquared) * (*susy)("VMIX", je*10+0) * (*susy)("VMIX", ie*10+0));

				if (ie == je)	{
					D90c += std::pow((*sm)("MASS", 24) / Mch[ie], 2.0) * X_UL[ie][ae][1] * X_UL[ie][ae][2] * h30(msqOverMchSquared);
				}
			}
		}
	}
	LOG_DEBUG("B0c1 : " + std::to_string(B0c1*1e10));
	LOG_DEBUG("B0c2 : " + doubleToString(B0c2*1e16, 20));
	for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) {
			for (int be = 0; be < 6; ++be) {
				double msqOverMchSquaredAe = std::pow(MsqU[ae] / Mch[ie], 2.0);
				double msqOverMchSquaredBe = std::pow(MsqU[be] / Mch[ie], 2.0);
				for (int ce = 0; ce < 3; ++ce) {
					C90c += X_UL[ie][be][1] * X_UL[ie][ae][2] * f40(msqOverMchSquaredAe, msqOverMchSquaredBe) * Gamma_UL[be][ce] * Gamma_UL[ae][ce];
				}
			}
		}
	}

	// Final adjustments
	B90c = -(B0c1 - B0c2) * kappa * std::pow((*sm)("MASS", 24), 2.0) / (2.0 * std::pow((*sm)("GAUGE",2), 2.0));
	B100c = (B0c1 + B0c2) * kappa * std::pow((*sm)("MASS", 24), 2.0) / (2.0 * std::pow((*sm)("GAUGE",2), 2.0));
	C90c *= -kappa / 8.0;
	D90c *= kappa;

    test = true;
	for (int ae = 0; ae < 6; ++ae) { // Conserve l'indexation à partir de 1
		if (!(std::fabs(MsqU[ae]) > (*sm)("MASS", 24) / 2. && std::fabs(MsqD[ae]) > (*sm)("MASS", 24) / 2.)) {
			test = false;
			break; // Sort de la boucle dès qu'une condition n'est pas remplie
		}
	}
	
}

void susy_parameters::reset_PrimeCQG(double Q_match) {
	if (is_PrimeCQG) {return;}
	is_PrimeCQG = true;
	Parameters* sm = Parameters::GetInstance();
	double mass_c_muW = (*sm).running_mass((*sm)("MASS", 4), (*sm)("MASS", 4), Q_match);
	LOG_INFO("mass_c_muW", mass_c_muW);
	MU = {(*sm)("MASS",2), mass_c_muW, mass_top_muW};
	LOG_INFO("MU_0", MU[0]);
	LOG_INFO("MU_1", MU[1]);
	LOG_INFO("MU_2", MU[2]);
	for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) { 
			for (int be = 0; be < 3; ++be) {
				X_UL[ie][ae][be] = 0.0;

				// Calculs pour X_UL et X_UR
				for (int ce = 0; ce < 3; ++ce) {
					X_UL[ie][ae][be] += -(*sm)("GAUGE",2) * (
						ag * (*susy)("VMIX", ie*10+0) * Gamma_UL[ae][ce] -
						aY * (*susy)("VMIX", ie*10+1) * Gamma_UR[ae][ce] * MU[ce] / (sqrt(2.0) * (*sm)("MASS", 24) * sinb)
					) * std::real(VCKM[ce][be]);

					G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*(*sm)("MASS",24)*(*susy)("VMIX", ie*10+0)*Gamma_UL[ae][ce]*ag-MU[ce]*(*susy)("VMIX", ie*10+1)*Gamma_UR[ae][ce]*aY)*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])/std::real(VCKM[2][1]));
				}

			}
		}
	}
	LOG_INFO("U1",(*susy)("VMIX", 1));
	LOG_INFO("U2",(*susy)("VMIX", 11));
	LOG_INFO("VCKM[2][2]", VCKM[2][2]);
	LOG_INFO("VCKM[2][2]", VCKM[2][1]);

	LOG_INFO("VCKM[0][0]", VCKM[0][0]);
	LOG_INFO("VCKM[0][1]", VCKM[0][1]);
	LOG_INFO("VCKM[0][2]", VCKM[0][2]);

	LOG_INFO("VCKM[1][0]", VCKM[1][0]);
	LOG_INFO("VCKM[1][1]", VCKM[1][1]);
	LOG_INFO("VCKM[1]2]", VCKM[1][2]);

	LOG_INFO("VCKM[2][0]", VCKM[2][0]);
	LOG_INFO("VCKM[2][1]", VCKM[2][1]);
	LOG_INFO("VCKM[2][2]", VCKM[2][2]);

	LOG_INFO("VCKM[ce][1]*VCKM[ce][1]", VCKM[1][1]*VCKM[0][2]);
	LOG_INFO("VCKM[ce][1]*VCKM[ce][1]", VCKM[1][1]*VCKM[0][2]);
	LOG_INFO("VCKM[be][2]", VCKM[1][2]);
}

void susy_parameters::reset_G() {
	if (!is_PrimeCQG) {return;}
	is_PrimeCQG=false;
	MU = {(*sm)("MASS",2), (*sm)("MASS",4), mass_top_muW}; 
	for (int ie = 0; ie < 2; ++ie) {
		for (int ae = 0; ae < 6; ++ae) { 
			for (int be = 0; be < 3; ++be) {
				X_UL[ie][ae][be] = 0.0;

				// Calculs pour X_UL et X_UR
				for (int ce = 0; ce < 3; ++ce) {
					X_UL[ie][ae][be] += -(*sm)("GAUGE",2) * (
						ag * (*susy)("VMIX", ie*10+0) * Gamma_UL[ae][ce] -
						aY * (*susy)("VMIX", ie*10+1) * Gamma_UR[ae][ce] * MU[ce] / (sqrt(2.0) * (*sm)("MASS", 24) * sinb)
					) * std::real(VCKM[ce][be]);

					G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*(*sm)("MASS",24)*(*susy)("VMIX", ie*10+0)*Gamma_UL[ae][ce]*ag-MU[ce]*(*susy)("VMIX", ie*10+1)*Gamma_UR[ae][ce]*aY)*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])*std::real(VCKM[2][1]));
				}

			}
		}
	}
}

susy_parameters* susy_parameters::instance = nullptr;