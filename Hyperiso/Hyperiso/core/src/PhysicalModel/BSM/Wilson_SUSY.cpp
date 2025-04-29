#include "Wilson_SUSY.h"

// TODO ALL FUNCTION AND RETURN

void C1_susy::NNLO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::SM, "MASS", 24},
		{ParameterType::BSM, "MASS", 1000002},
		{ParameterType::BSM, "MASS", 1000004},
		{ParameterType::BSM, "MASS", 1000006},
		{ParameterType::BSM, "MASS", 2000002},
		{ParameterType::BSM, "MASS", 2000004},
		{ParameterType::BSM, "MASS", 2000006},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {15, 0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {15, 1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {15, 2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {15, 3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {15, 4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {15, 5}}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        complex_t C1squark_2 = 0.0;
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		Array1D_7 MsqU = {src.at({ParameterType::BSM, "MASS", 1000002})->get_val(), src.at({ParameterType::BSM, "MASS", 1000004})->get_val(), src.at({ParameterType::BSM, "MASS", 1000006})->get_val(), 
		src.at({ParameterType::BSM, "MASS", 2000002})->get_val(), src.at({ParameterType::BSM, "MASS", 2000004})->get_val(), src.at({ParameterType::BSM, "MASS", 2000006})->get_val()};

		if (std::all_of(begin(MsqU), end(MsqU), [&](double m) { return std::abs(m) > mW / 2.0; })) {
			
			C1squark_2 = -208.0 / 3.0;
			for (int ae = 0; ae < 6; ++ae) {
				double xsqa = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::SM, "MASS", 24})->get_val(), 2.0);
				if (4.0 * xsqa > 1.0) {
					double angle = 2.0 * asin(0.5 / sqrt(xsqa));
					C1squark_2 += -2.0 * std::pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
				}
				C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;

				xsqa = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {15, ae}})->get_val() / src.at({ParameterType::SM, "MASS", 24})->get_val(), 2.0);
				if (4.0 * xsqa > 1.0) {
					double angle = 2.0 * asin(0.5 / sqrt(xsqa));
					C1squark_2 += -2.0 * std::pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
				}
				C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;
			}
		}
        dep_param->set_expected(C1squark_2);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3040405, 6161, 2, 1)}, sources, func);

    // complex_t C1squark_2 = 0.0;
	// if (std::all_of(begin((*sus_param).MsqU), end((*sus_param).MsqU), [&](double m) { return std::abs(m) > sm("MASS", 24) / 2.0; })) {
		
	// 	C1squark_2 = -208.0 / 3.0;
	// 	for (int ae = 0; ae < 6; ++ae) {
	// 		double xsqa = std::pow((*sus_param).MsqU[ae] / sm("MASS", 24), 2.0);
	// 		if (4.0 * xsqa > 1.0) {
	// 			double angle = 2.0 * asin(0.5 / sqrt(xsqa));
	// 			C1squark_2 += -2.0 * std::pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
	// 		}
	// 		C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;

	// 		xsqa = std::pow((*sus_param).MsqD[ae] / sm("MASS", 24), 2.0);
	// 		if (4.0 * xsqa > 1.0) {
	// 			double angle = 2.0 * asin(0.5 / sqrt(xsqa));
	// 			C1squark_2 += -2.0 * std::pow(4.0 * xsqa - 1.0, 1.5) * Cl2(angle);
	// 		}
	// 		C1squark_2 += 8.0 * (xsqa - 1.0 / 3.0) * log(xsqa) + 16.0 * xsqa;
	// 	}
	// }
    // this->set_WilsonCoeffMatching("NNLO", C1squark_2);
    // return C1squark_2;
}
void C3_susy::NNLO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "EW_SCALE", 1},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13, 0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13, 1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 0, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 1, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 2, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 3, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 4, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 0, 5, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 0, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 1, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 2, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 3, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 4, 2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3, 1, 5, 2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},
		{ParameterType::BSM, "MASS", 37}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		complex_t C3charg_2 = 0.0;
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		std::shared_ptr<susy_parameters> sus_param;
		for(int ie = 0; ie < 2; ie++) {
			for(int ae = 0; ae < 6; ae++) {
				double ratio_mass_W_Mch = pow(src.at({ParameterType::SM, "MASS", 24})->get_val()/ src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double ratio_MsqU_Mch = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double log_mu_W_MsqU = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2.0));
				
				C3charg_2 += ratio_mass_W_Mch * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * h71(ratio_MsqU_Mch, log_mu_W_MsqU);
			}
		}

		C3charg_2 *= src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		double C3H_2 = G3H(src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val(),src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val())+Delta3H(src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val(),src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val()) // 7 -> lu
						*log(pow(Q_match/src.at({ParameterType::BSM, "MASS", 37})->get_val(),2.));

        dep_param->set_expected(C3charg_2 + C3H_2);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 4133, 2, 1)}, sources, func);


    // complex_t C3charg_2 = 0.0;

	// for(int ie = 0; ie < 2; ie++) {
	// 	for(int ae = 0; ae < 6; ae++) {
	// 		double ratio_mass_W_Mch = std::pow(sm("MASS",24)/ (*sus_param).Mch[ie], 2.0);
	// 		double ratio_MsqU_Mch = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
	// 		double log_mu_W_MsqU = std::log(std::pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2.0));
	// 		C3charg_2 += ratio_mass_W_Mch * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ae][2] * h71(ratio_MsqU_Mch, log_mu_W_MsqU);
	// 	}
	// }

	// C3charg_2 *= (*sus_param).kappa;
	// double C3H_2 = G3H(sus_param->yt,sus_param->lu)+Delta3H(sus_param->yt,sus_param->lu)*log(pow(this->get_Q_match()/(*susy)("MASS",37),2.));
    // this->set_WilsonCoeffMatching("NNLO", C3charg_2 + C3H_2);
    // // return C3charg_2 + C3H_2;
}

void C4_susy::NLO_calculation() {

	std::unordered_set<ParamId> sources {
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_MATCH_BSM", 1},
        {"EW_SCALE", 1},
        {ParameterType::BSM, "MASS", 37}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		complex_t C4charg_1 = 0.;

		std::shared_ptr<susy_parameters> sus_param;
		for (int ie = 0; ie < 2; ie++) {
			for (int ae = 0; ae < 6; ae++) {
				C4charg_1+= pow(src.at({ParameterType::SM, "MASS", 24})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)
							*(src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()*h40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
			}
		}
		C4charg_1*=src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		complex_t C4H_1 = EH(src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val(),src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val()); //7 -> lu; 1 -> yt
		complex_t coeff_temp = C4H_1+C4charg_1;

        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 6153, 1, 1)}, sources, func);

    // complex_t C4charg_1 = 0.;
    // for (int ie = 0; ie < 2; ie++) {

	// 	for (int ae = 0; ae < 6; ae++) {
	// 		C4charg_1+= pow(sm("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
    //     }
    // }
    // C4charg_1*=(*sus_param).kappa;
    // complex_t C4H_1 = EH(sus_param->yt,sus_param->lu);
    // complex_t coeff_temp = C4H_1+C4charg_1;
    // this->set_WilsonCoeffMatching("NLO", coeff_temp);
    // // return coeff_temp;
}

void C4_susy::NNLO_calculation() {

	std::unordered_set<ParamId> sources {
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_MATCH_BSM", 1},
        {"EW_SCALE", 1},
        {ParameterType::BSM, "MASS", 37}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		complex_t C4charg_2{};
		complex_t C4four_2{};
		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		for(int ie = 0; ie < 2; ie++) {
			for(int ae = 0; ae < 6; ae++) {
				double ratio_mass_W_Mch = pow(mW/ src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double ratio_MsqU_Mch = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double log_mu_W_MsqU = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2.0));
				C4charg_2 += ratio_mass_W_Mch * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * h41(ratio_MsqU_Mch, log_mu_W_MsqU);
				for(int be = 0; be < 6; be++) {
					for(int ce = 0; ce < 6; ce++) {
						
						C4four_2 += ratio_mass_W_Mch * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, be}})->get_val() * src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,be, ce}})->get_val() *
									(1.0 + log_mu_W_MsqU) * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 2}})->get_val() *
									q61(ratio_MsqU_Mch, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0));
					}
				}
			}
		}

		C4charg_2 *= src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		C4four_2 *= src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		double C4H_2=G4H(yt,lu)+Delta4H(yt,lu)*log(pow(Q_match/mH,2.));

        dep_param->set_expected(C4charg_2+C4four_2+C4H_2);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 6153, 2, 1)}, sources, func);

    // complex_t C4charg_2{};
	// complex_t C4four_2{};

    // for(int ie = 0; ie < 2; ie++) {
	// 	for(int ae = 0; ae < 6; ae++) {
	// 		double ratio_mass_W_Mch = std::pow(sm("MASS",24)/ (*sus_param).Mch[ie], 2.0);
	// 		double ratio_MsqU_Mch = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
	// 		double log_mu_W_MsqU = std::log(std::pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2.0));
	// 		C4charg_2 += ratio_mass_W_Mch * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ae][2] * h41(ratio_MsqU_Mch, log_mu_W_MsqU);
	// 		for(int be = 0; be < 6; be++) {
	// 			for(int ce = 0; ce < 6; ce++) {
	// 				C4four_2 += ratio_mass_W_Mch * (*sus_param).P_U[ae][be] * (*sus_param).MsqU[be] / (*sus_param).Mch[ie] * (*sus_param).P_U[be][ce] *
	// 							(1.0 + log_mu_W_MsqU) * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ce][2] *
	// 							q61(ratio_MsqU_Mch, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0));
	// 			}
	// 		}
	// 	}
	// }

    // C4charg_2 *= (*sus_param).kappa;
	// C4four_2 *= (*sus_param).kappa;
	// double C4H_2=G4H(sus_param->yt,sus_param->lu)+Delta4H(sus_param->yt,sus_param->lu)*log(pow(this->get_Q_match()/(*susy)("MASS",37),2.));
    // this->set_WilsonCoeffMatching("NNLO", C4charg_2+C4four_2+C4H_2);
    // // return C4charg_2+C4four_2+C4H_2;
}

void C5_susy::NNLO_calculation() {
	
	std::unordered_set<ParamId> sources {
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_MATCH_BSM", 1},
        {"EW_SCALE", 1},
        {ParameterType::BSM, "MASS", 37}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		complex_t C3charg_2{};
		complex_t C4charg_1{};

		for(int ie = 0; ie < 2; ie++) {
			for(int ae = 0; ae < 6; ae++) {
				double ratio_mass_W_Mch = pow(mW/ src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double ratio_MsqU_Mch = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double log_mu_W_MsqU = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2.0));

				C3charg_2 += ratio_mass_W_Mch * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * h71(ratio_MsqU_Mch, log_mu_W_MsqU);

				C4charg_1+= pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*(src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val()
							*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()*h40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
			}
		}
		C3charg_2 *= src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		C4charg_1 *= src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		complex_t C5charg_2 = -C3charg_2 / 10.0 + 2.0 / 15.0 * C4charg_1;
		double C4H_1=EH(yt,lu);
		double C3H_2=G3H(yt,lu)+Delta3H(yt,lu)*log(pow(Q_match/mH,2.));
		double C5H_2=-C3H_2/10.+2./15.*C4H_1;

        dep_param->set_expected(C5charg_2+C5H_2);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 4536, 2, 1)}, sources, func);

    // complex_t C3charg_2{};
	// complex_t C4charg_1{};

	// for(int ie = 0; ie < 2; ie++) {
	// 	for(int ae = 0; ae < 6; ae++) {
	// 		double ratio_mass_W_Mch = std::pow(sm("MASS",24)/ (*sus_param).Mch[ie], 2.0);
	// 		double ratio_MsqU_Mch = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
	// 		double log_mu_W_MsqU = std::log(std::pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2.0));

	// 		C3charg_2 += ratio_mass_W_Mch * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ae][2] * h71(ratio_MsqU_Mch, log_mu_W_MsqU);

	// 		C4charg_1+= pow(sm("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
	// 	}
	// }
	// C3charg_2 *= (*sus_param).kappa;
	// C4charg_1 *= (*sus_param).kappa;
	// complex_t C5charg_2 = -C3charg_2 / 10.0 + 2.0 / 15.0 * C4charg_1;
	// double C4H_1=EH(sus_param->yt,sus_param->lu);
    // double C3H_2=G3H(sus_param->yt,sus_param->lu)+Delta3H(sus_param->yt,sus_param->lu)*log(pow(this->get_Q_match()/(*susy)("MASS",37),2.));
	// double C5H_2=-C3H_2/10.+2./15.*C4H_1;
    // this->set_WilsonCoeffMatching("NNLO", C5charg_2+C5H_2);
    // // return C5charg_2+C5H_2;
}

void C6_susy::NNLO_calculation() {

	std::unordered_set<ParamId> sources {
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_MATCH_BSM", 1},
        {"EW_SCALE", 1},
        {ParameterType::BSM, "MASS", 37}
    };

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		complex_t C3charg_2{};
		complex_t C4charg_1{};

		for(int ie = 0; ie < 2; ie++) {
			for(int ae = 0; ae < 6; ae++) {
				double ratio_mass_W_Mch = pow(mW/ src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double ratio_MsqU_Mch = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double log_mu_W_MsqU = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2.0));

				C3charg_2 += ratio_mass_W_Mch * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * h71(ratio_MsqU_Mch, log_mu_W_MsqU);

				C4charg_1+= pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*(src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()*h40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
			}
		}
		C3charg_2 *= src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		C4charg_1 *= src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		complex_t C6charg_2 = -3.0 / 16.0 * C3charg_2 + 1.0 / 4.0 * C4charg_1;
		complex_t C4H_1=EH(yt,lu);
		complex_t C3H_2=G3H(yt,lu)+Delta3H(yt,lu)*log(pow(Q_match/mH,2.));
		complex_t C6H_2=-3./16.*C3H_2+1./4.*C4H_1;

        dep_param->set_expected(C6charg_2+C6H_2);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3050707, 6556, 2, 1)}, sources, func);

    // complex_t C3charg_2{};
	// complex_t C4charg_1{};

	// for(int ie = 0; ie < 2; ie++) {
	// 	for(int ae = 0; ae < 6; ae++) {
	// 		double ratio_mass_W_Mch = std::pow(sm("MASS",24)/ (*sus_param).Mch[ie], 2.0);
	// 		double ratio_MsqU_Mch = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
	// 		double log_mu_W_MsqU = std::log(std::pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2.0));

	// 		C3charg_2 += ratio_mass_W_Mch * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ae][2] * h71(ratio_MsqU_Mch, log_mu_W_MsqU);

	// 		C4charg_1+= pow(sm("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
	// 	}
	// }
	// C3charg_2 *= (*sus_param).kappa;
	// C4charg_1 *= (*sus_param).kappa;
    // complex_t C6charg_2 = -3.0 / 16.0 * C3charg_2 + 1.0 / 4.0 * C4charg_1;
	// complex_t C4H_1=EH(sus_param->yt,sus_param->lu);
    // complex_t C3H_2=G3H(sus_param->yt,sus_param->lu)+Delta3H(sus_param->yt,sus_param->lu)*log(pow(this->get_Q_match()/(*susy)("MASS",37),2.));
	// complex_t C6H_2=-3./16.*C3H_2+1./4.*C4H_1;
    // this->set_WilsonCoeffMatching("NNLO", C6charg_2+C6H_2);
    // // return C6charg_2+C6H_2;
}


void C7_susy::LO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::SM, "MASS", 24},                          // utilisé directement
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},              // pour kappa
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},              // pour lu
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13, 0}},        // utilisé dans boucle sur ie
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13, 1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 0}},        // utilisé dans boucle sur ae
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14, 5}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},         // matrice pour ie/ae
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}            // pour yt
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
		double alpha = src.at({ParameterType::BSM, "ALPHA", 0})->get_val();
		double mH0 = src.at({ParameterType::BSM, "MASS", 35})->get_val();
		double mh = src.at({ParameterType::BSM, "MASS", 25})->get_val();
		double mA0 = src.at({ParameterType::BSM, "MASS", 36})->get_val();

		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();

		double epsilon0 = src.at({ParameterType::WILSON, "EPSILON_SUSY", {0,1}})->get_val();
		double epsilon0p = src.at({ParameterType::WILSON, "EPSILON_SUSY", {0,2}})->get_val();
		double epsilon1p = src.at({ParameterType::WILSON, "EPSILON_SUSY", 1})->get_val();
		double epsilon2 = src.at({ParameterType::WILSON, "EPSILON_SUSY", 2})->get_val();
		double epsilonb = src.at({ParameterType::WILSON, "EPSILON_SUSY", 3})->get_val();
		double epsilonbp = src.at({ParameterType::WILSON, "EPSILON_SUSY", 4})->get_val();
		src.at({ParameterType::BSM, "AMIX", 00})->get_val();
		double H03  = 0.; //TODO, wtf
		double A02 = 0.; //TODO, WTF
		complex_t C7SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*tanb)*tanb*F7_2(xt);
		complex_t C7Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*tanb)*tanb*F7_2(yt);

		complex_t C7Heps2_0=0.;

		if((A02==0.)&&(H03==0.)) {
			C7Heps2_0=-epsilon2*epsilon1p*pow(tanb,2.)/(1.+epsilonb*tanb)/(1.+epsilon0*tanb)*F7_2(yt);
			C7Heps2_0+=epsilon2/pow(1.+epsilonb*tanb,2.)*(1.+pow(tanb,2.))/(1.+epsilon0*tanb)/72.		*((cos(alpha)+sin(alpha)*tanb)*(-sin(alpha)+epsilonb*cos(alpha))*pow(mass_b_muW/mh,2.)
			+(sin(alpha)-cos(alpha)*tanb)*(cos(alpha)+epsilonb*sin(alpha))*pow(mass_b_muW/mH0,2.)			+(-cos(atan(tanb))-sin(atan(tanb))*tanb)*(sin(atan(tanb))-epsilonb*cos(atan(tanb)))*pow(mass_b_muW/mA0,2.));

		}
		else {		
			C7Heps2_0=-epsilon2*epsilon1p*pow(tanb,2.)/(1.+epsilonb*tanb)/(1.+epsilon0*tanb)*F7_2(yt);
			C7Heps2_0+=epsilon2/pow(1.+epsilonb*tanb,2.)*(1.+pow(tanb,2.))/(1.+epsilon0*tanb)/72.	*((src.at({ParameterType::BSM, "HMIX", 00})->get_val()+src.at({ParameterType::BSM, "HMIX", 01})->get_val()*tanb)*(-src.at({ParameterType::BSM, "HMIX", 01})->get_val()+epsilonb*src.at({ParameterType::BSM, "HMIX", 00})->get_val())*pow(mass_b_muW/mh,2.)
			+(src.at({ParameterType::BSM, "HMIX", 10})->get_val()+src.at({ParameterType::BSM, "HMIX", 11})->get_val()*tanb)*(-src.at({ParameterType::BSM, "HMIX", 11})->get_val()+epsilonb*src.at({ParameterType::BSM, "HMIX", 10})->get_val())*pow(mass_b_muW/mH0,2.)
			+(src.at({ParameterType::BSM, "HMIX", 20})->get_val()+src.at({ParameterType::BSM, "HMIX", 21})->get_val()*tanb)*(-src.at({ParameterType::BSM, "HMIX", 21})->get_val()+epsilonb*src.at({ParameterType::BSM, "HMIX", 20})->get_val())*pow(mass_b_muW/H03,2.)

			+(src.at({ParameterType::BSM, "AMIX", 00})->get_val()+src.at({ParameterType::BSM, "AMIX", 01})->get_val()*tanb)*(-src.at({ParameterType::BSM, "AMIX", 01})->get_val()+epsilonb*src.at({ParameterType::BSM, "AMIX", 00})->get_val())*pow(mass_b_muW/mA0,2.) //mass_A0 = 36 ? = HO3 ?
			+(src.at({ParameterType::BSM, "AMIX", 10})->get_val()+src.at({ParameterType::BSM, "AMIX", 11})->get_val()*tanb)*(-src.at({ParameterType::BSM, "AMIX", 11})->get_val()+epsilonb*src.at({ParameterType::BSM, "AMIX", 10})->get_val())*pow(mass_b_muW/A02,2.));

		}

		complex_t C7charg_0 = 0.0;
		complex_t C7_chargeps_0 = 0.0;

		std::function<double(std::function<double(double)>, int, int, int, int, bool)> calculateContribution = [src, mW, tanb, epsilonb](auto hFunc, int X, int X2, int ie, int ae, bool isChargeps) -> double {
            complex_t ratio = std::pow(mW / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2);
            double msqOverMchSquared = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
            double factor = isChargeps ? (-epsilonb / (1.0 + epsilonb * tanb) * tanb) : 1.0;
            return ratio * (
			src.at({ParameterType::WILSON, "MATRIX_BSM", {X,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {X2,ie, ae, 2}})->get_val() * hFunc(msqOverMchSquared)) * src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 18})->get_val() * factor; //18->kappafactor
        };

		for (int ie = 0; ie < 2; ++ie) {
			for (int ae = 0; ae < 6; ++ae) {
				
				C7charg_0 += calculateContribution(h10, 3,3, ie, ae, false) + src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/mass_b_muW *calculateContribution(h20, 3, 4, ie, ae, false);
				C7_chargeps_0 += src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/mass_b_muW *calculateContribution(h20, 3, 4, ie, ae, true);

			}
		}

		complex_t C7H_0 = 1./3.*lu*lu*F7_1(yt) - lu*ld*F7_2(yt);

        dep_param->set_expected(C7SMeps_0 + C7Heps_0 + C7Heps2_0 + C7charg_0 + C7_chargeps_0+C7H_0);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 4422, 0, 1)}, sources, func);

    // complex_t C7SMeps_0= ((*sus_param).epsilonb-(*sus_param).epsilonbp)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))*(*susy)("HMIX",2)*F7_2((*sus_param).xt);
    // complex_t C7Heps_0=(-(*sus_param).epsilon0p-(*sus_param).epsilonb)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))*(*susy)("HMIX",2)*F7_2((*sus_param).yt);

    // complex_t C7Heps2_0=0.;

    // if(((*sus_param).mass_A02==0.)&&((*sus_param).mass_H03==0.)) {
    //     C7Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("HMIX",2),2.)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))*F7_2((*sus_param).yt);
    //     C7Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("HMIX",2),2.)*(1.+pow((*susy)("HMIX",2),2.))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))/72.		*((cos((*susy)("ALPHA",0))+sin((*susy)("ALPHA",0))*(*susy)("HMIX",2))*(-sin((*susy)("ALPHA",0))+(*sus_param).epsilonb*cos((*susy)("ALPHA",0)))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",25),2.)
    //     +(sin((*susy)("ALPHA",0))-cos((*susy)("ALPHA",0))*(*susy)("HMIX",2))*(cos((*susy)("ALPHA",0))+(*sus_param).epsilonb*sin((*susy)("ALPHA",0)))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",35),2.)			+(-cos(atan((*susy)("HMIX",2)))-sin(atan((*susy)("HMIX",2)))*(*susy)("HMIX",2))*(sin(atan((*susy)("HMIX",2)))-(*sus_param).epsilonb*cos(atan((*susy)("HMIX",2))))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",36),2.));

    // }
	// else {		
	// 	C7Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("HMIX",2),2.)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))*F7_2((*sus_param).yt);
	// 	C7Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("HMIX",2),2.)*(1.+pow((*susy)("HMIX",2),2.))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))/72.	*(((*susy)("HMIX", 00)+(*susy)("HMIX", 01)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 01)+(*sus_param).epsilonb*(*susy)("HMIX", 00))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",25),2.)
	// 	+((*susy)("HMIX", 10)+(*susy)("HMIX", 11)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 11)+(*sus_param).epsilonb*(*susy)("HMIX", 10))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",35),2.)
	// 	+((*susy)("HMIX", 20)+(*susy)("HMIX", 21)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 21)+(*sus_param).epsilonb*(*susy)("HMIX", 20))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).mass_H03,2.)

	// 	+((*susy)("AMIX", 00)+(*susy)("AMIX", 01)*(*susy)("HMIX",2))*(-(*susy)("AMIX", 01)+(*sus_param).epsilonb*(*susy)("AMIX", 00))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",36),2.) //mass_A0 = 36 ? = HO3 ?
	// 	+((*susy)("AMIX", 10)+(*susy)("AMIX", 11)*(*susy)("HMIX",2))*(-(*susy)("AMIX", 11)+(*sus_param).epsilonb*(*susy)("AMIX", 10))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).mass_A02,2.));

    // }

    // complex_t C7charg_0 = 0.0;
	// complex_t C7_chargeps_0 = 0.0;

	// for (int ie = 0; ie < 2; ++ie) {
	// 	for (int ae = 0; ae < 6; ++ae) {
			
	// 		C7charg_0 += calculateContribution(h10, (*sus_param).X_UL,(*sus_param).X_UL, ie, ae, false) + (*sus_param).Mch[ie]/wilson_p("WPARAM_MATCH_SM", {5,1}) *calculateContribution(h20, (*sus_param).X_UL, (*sus_param).X_UR, ie, ae, false);
	// 		C7_chargeps_0 += (*sus_param).Mch[ie]/wilson_p("WPARAM_MATCH_SM", {5,1}) *calculateContribution(h20, (*sus_param).X_UL, (*sus_param).X_UR, ie, ae, true);

	// 	}
	// }

    // complex_t C7H_0 = 1./3.*sus_param->lu*sus_param->lu*F7_1(sus_param->yt) - sus_param->lu*sus_param->ld*F7_2(sus_param->yt);

    // this->set_WilsonCoeffMatching("LO", C7SMeps_0 + C7Heps_0 + C7Heps2_0 + C7charg_0 + C7_chargeps_0+C7H_0);
    // // return C7SMeps_0 + C7Heps_0 + C7Heps2_0 + C7charg_0 + C7_chargeps_0 + C7H_0;

}

void C7_susy::NLO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},   // kappa
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},   // lu
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},   // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, // yt
		{ParameterType::WILSON, "EW_SCALE", 1},     // Q_match
		{ParameterType::SM, "MASS", 24},               // mW
		{ParameterType::BSM, "MASS", 25},              // mh
		{ParameterType::BSM, "MASS", 35},              // mH0
		{ParameterType::BSM, "MASS", 37},              // mH
	
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,5}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,5}}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double mH0 = src.at({ParameterType::BSM, "MASS", 35})->get_val();
		double mh = src.at({ParameterType::BSM, "MASS", 25})->get_val();

		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();

		double H03  = 0.; //TODO, wtf
		double A02 = 0.; //TODO, WTF
		complex_t C7charg_1{};
		complex_t C7four_1{};
		for (int ie = 0; ie < 2; ie++) {

			for (int ae = 0; ae < 6; ae++) {
				C7charg_1+=pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*(src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()
							*h11(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),log(pow(Q_match/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ae}})->get_val(),2.))) 
							+ src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/mass_b_muW*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()
							*h21(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),log(pow(Q_match/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ae}})->get_val(),2.))));

				for (int ce = 0; ce < 6; ce++) {
					double log_mu_W_MsqU_ce = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ce}})->get_val(), 2.0));
					double MsqU_be_Mch_ie_ratio = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
					
					for (int de = 0; de < 6; de++) {
						double ratio_MsqU_ae_Mch_ie = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
						double ratio_MsqU_de_Mch_ie = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, de}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);

						C7four_1 += std::pow(mW / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0) * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce}})->get_val() * MsqU_be_Mch_ie_ratio * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ce, de}})->get_val() * (1.0 + log_mu_W_MsqU_ce) *
									(src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 2}})->get_val() * (-q11(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) + 2.0 / 3.0 * q21(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie)) +
									src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/ mass_b_muW * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() 
									* src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, de, 2}})->get_val() * (-q31(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) + 2.0 / 3.0 * q41(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie)));

					}
				}
			}
		}
		C7charg_1*=-0.5*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val(); // 6 -> kappa
		complex_t C7H_1 = G7H(yt,lu,ld)+Delta7H(yt,lu,ld)*log(pow(Q_match/mH,2.))
		-4./9.*EH(yt,lu);

        dep_param->set_expected(C7charg_1+C7four_1+C7H_1);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 4422, 1, 1)}, sources, func);

    // complex_t C7charg_1{};
    // complex_t C7four_1{};
    // for (int ie = 0; ie < 2; ie++) {

	// 	for (int ae = 0; ae < 6; ae++) {
	// 		C7charg_1+=pow(sm("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h11(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(this->get_Q_match()/(*sus_param).MsqU[ae],2.))) + (*sus_param).Mch[ie]/sus_param->mass_b_muW*(*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h21(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(this->get_Q_match()/(*sus_param).MsqU[ae],2.))));

	// 		for (int ce = 0; ce < 6; ce++) {

	// 			double log_mu_W_MsqU_ce = std::log(std::pow(this->get_Q_match() / (*sus_param).MsqU[ce], 2.0));
	// 			double MsqU_be_Mch_ie_ratio = (*sus_param).MsqU[ce] / (*sus_param).Mch[ie];
				
	// 			for (int de = 0; de < 6; de++) {
	// 				double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
	// 				double ratio_MsqU_de_Mch_ie = std::pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0);

	// 				C7four_1 += std::pow(sm("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).P_U[ae][ce] * MsqU_be_Mch_ie_ratio * (*sus_param).P_U[ce][de] * (1.0 + log_mu_W_MsqU_ce) *
	// 							((*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][de][2] * (-q11(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) + 2.0 / 3.0 * q21(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie)) +
	// 							(*sus_param).Mch[ie] / wilson_p("WPARAM_MATCH_SM", {5,1}) * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UR[ie][de][2] * (-q31(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) + 2.0 / 3.0 * q41(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie)));

	// 			}
	// 		}
	// 	}
	// }
    // C7charg_1*=-0.5*(*sus_param).kappa;
    // complex_t C7H_1 = G7H(sus_param->yt,sus_param->lu,sus_param->ld)+Delta7H(sus_param->yt,sus_param->lu,sus_param->ld)*log(pow(this->get_Q_match()/sus_param->m_H,2.))
    // -4./9.*EH(sus_param->yt,sus_param->lu);

    // this->set_WilsonCoeffMatching("NLO",C7charg_1+C7four_1+C7H_1);
    // // return C7charg_1+C7four_1 + C7H_1;

}

void C7_susy::NNLO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},     // lu
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},     // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},  // yt
		{ParameterType::WILSON, "EW_SCALE", 1},          // Q_match
		{ParameterType::WILSON, "WPARAM_MATCH_SM", 6}    // mass_top_muW
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();

		double coeff_temp =C7H2(yt,lu,ld,log(pow(Q_match/mass_top_muW, 2.)));
        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 4422, 2, 1)}, sources, func);

	// double coeff_temp =C7H2(sus_param->yt,sus_param->lu,sus_param->ld,log(pow(this->get_Q_match()/sus_param->mass_top_muW, 2.)));
    // // return this->double_to_complex_save("NNLO", coeff_temp);
}

void C8_susy::LO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},      // lu
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},      // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},   // yt
		{ParameterType::WILSON, "EW_SCALE", 1},           // Q_match
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}}, // xt
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}}, // mass_b_muW
		{ParameterType::SM, "MASS", 24},                  // mW
		{ParameterType::BSM, "MASS", 25},                 // mh
		{ParameterType::BSM, "MASS", 35},                 // mH0
		{ParameterType::BSM, "MASS", 36},                 // mA0
		{ParameterType::BSM, "MASS", 37},                 // mH
		{ParameterType::BSM, "HMIX", 2},                  // tanb
		{ParameterType::BSM, "ALPHA", 0},                 // alpha
	
		// EPSILON_SUSY parameters
		{ParameterType::WILSON, "EPSILON_SUSY", {0,1}},
		{ParameterType::WILSON, "EPSILON_SUSY", {0,2}},
		{ParameterType::WILSON, "EPSILON_SUSY", 1},
		{ParameterType::WILSON, "EPSILON_SUSY", 2},
		{ParameterType::WILSON, "EPSILON_SUSY", 3},
		{ParameterType::WILSON, "EPSILON_SUSY", 4},
	
		// HMIX & AMIX matrix elements used in the else branch
		{ParameterType::BSM, "HMIX", 00},
		{ParameterType::BSM, "HMIX", 01},
		{ParameterType::BSM, "HMIX", 10},
		{ParameterType::BSM, "HMIX", 11},
		{ParameterType::BSM, "HMIX", 20},
		{ParameterType::BSM, "HMIX", 21},
		{ParameterType::BSM, "HMIX", 31},
		{ParameterType::BSM, "AMIX", 00},
		{ParameterType::BSM, "AMIX", 01},
		{ParameterType::BSM, "AMIX", 10},
		{ParameterType::BSM, "AMIX", 11},
	
		// Additional WPARAM_SI_BSM entries for calculateContribution
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 18},     // kappafactor
	
		// Matrices MATRIX_BSM (h50 and h60 contributions)
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		complex_t C4charg_2{};
		complex_t C4four_2{};
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
		double alpha = src.at({ParameterType::BSM, "ALPHA", 0})->get_val();
		double mH0 = src.at({ParameterType::BSM, "MASS", 35})->get_val();
		double mh = src.at({ParameterType::BSM, "MASS", 25})->get_val();
		double mA0 = src.at({ParameterType::BSM, "MASS", 36})->get_val();

		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();

		double epsilon0 = src.at({ParameterType::WILSON, "EPSILON_SUSY", {0,1}})->get_val();
		double epsilon0p = src.at({ParameterType::WILSON, "EPSILON_SUSY", {0,2}})->get_val();
		double epsilon1p = src.at({ParameterType::WILSON, "EPSILON_SUSY", 1})->get_val();
		double epsilon2 = src.at({ParameterType::WILSON, "EPSILON_SUSY", 2})->get_val();
		double epsilonb = src.at({ParameterType::WILSON, "EPSILON_SUSY", 3})->get_val();
		double epsilonbp = src.at({ParameterType::WILSON, "EPSILON_SUSY", 4})->get_val();

		double H03  = 0.; //TODO, wtf
		double A02 = 0.; //TODO, WTF
		complex_t C8SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*tanb)*tanb*F8_2(xt);
		complex_t C8Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*tanb)*tanb*F8_2(yt);
		complex_t C8Heps2_0=0.;

		if((A02==0.)&&(H03==0.)) {

			C8Heps2_0=-epsilon2*epsilon1p*pow(tanb,2.)/(1.+epsilonb*tanb)/(1.+epsilon0*tanb)*F8_2(yt);
			C8Heps2_0+=epsilon2/pow(1.+epsilonb*tanb,2.)*(1.+pow(tanb,2.))/(1.+epsilon0*tanb)/72.		*((cos(alpha)+sin(alpha)*tanb)*(-sin(alpha)
					+epsilonb*cos(alpha))*pow(mass_b_muW/mh,2.)
					+(sin(alpha)-cos(alpha)*tanb)*(cos(alpha)+epsilonb*sin(alpha))*pow(mass_b_muW/mH0,2.)
					+(-cos(atan(tanb))-sin(atan(tanb))*tanb)*(sin(atan(tanb))-epsilonb*cos(atan(tanb)))*pow(mass_b_muW/mA0,2.));
		}
		else {	
			C8Heps2_0=-epsilon2*epsilon1p*pow(tanb,2.)/(1.+epsilonb*tanb)/(1.+epsilon0*tanb)*F8_2(yt);
			C8Heps2_0+=-3.*epsilon2/pow(1.+epsilonb*tanb,2.)*(1.+pow(tanb,2.))/(1.+epsilon0*tanb)/72.
			*((src.at({ParameterType::BSM, "HMIX", 00})->get_val()+src.at({ParameterType::BSM, "HMIX", 01})->get_val()*tanb)*(-src.at({ParameterType::BSM, "HMIX", 01})->get_val()+epsilonb*src.at({ParameterType::BSM, "HMIX", 00})->get_val())*pow(mass_b_muW/mh,2.)
			+(src.at({ParameterType::BSM, "HMIX", 10})->get_val()+src.at({ParameterType::BSM, "HMIX", 11})->get_val()*tanb)*(-src.at({ParameterType::BSM, "HMIX", 11})->get_val()+epsilonb*src.at({ParameterType::BSM, "HMIX", 10})->get_val())*pow(mass_b_muW/mH0,2.)
			+(src.at({ParameterType::BSM, "HMIX", 31})->get_val()+src.at({ParameterType::BSM, "HMIX", 21})->get_val()*tanb)*(-src.at({ParameterType::BSM, "HMIX", 21})->get_val()+epsilonb*src.at({ParameterType::BSM, "HMIX", 20})->get_val())*pow(mass_b_muW/H03,2.)		

			+(src.at({ParameterType::BSM, "AMIX", 00})->get_val()+src.at({ParameterType::BSM, "AMIX", 01})->get_val()*tanb)*(-src.at({ParameterType::BSM, "AMIX", 01})->get_val()+epsilonb*src.at({ParameterType::BSM, "AMIX", 00})->get_val())*pow(mass_b_muW/mA0,2.)
			+(src.at({ParameterType::BSM, "AMIX", 10})->get_val()+src.at({ParameterType::BSM, "AMIX", 11})->get_val()*tanb)*(-src.at({ParameterType::BSM, "AMIX", 11})->get_val()+epsilonb*src.at({ParameterType::BSM, "AMIX", 10})->get_val())*pow(mass_b_muW/A02,2.));
			}

		std::function<double(std::function<double(double)>, int, int, int, int, bool)> calculateContribution = [src, mW, tanb, epsilonb](auto hFunc, int X, int X2, int ie, int ae, bool isChargeps) -> double {
			complex_t ratio = std::pow(mW / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2);
			double msqOverMchSquared = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
			double factor = isChargeps ? (-epsilonb / (1.0 + epsilonb * tanb) * tanb) : 1.0;
			return ratio * (
			src.at({ParameterType::WILSON, "MATRIX_BSM", {X,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {X2,ie, ae, 2}})->get_val() * hFunc(msqOverMchSquared)) * src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 18})->get_val() * factor; //18->kappafactor
		};

		complex_t C8charg_0 = 0.0;
		complex_t C8_chargeps_0 = 0.0;

		for (int ie = 0; ie < 2; ++ie) {
			for (int ae = 0; ae < 6; ++ae) {
				C8charg_0 += calculateContribution(h50, 3, 3, ie, ae, false) + src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/mass_b_muW*calculateContribution(h60,3,  4, ie, ae, false);
				C8_chargeps_0 += src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/mass_b_muW *calculateContribution(h60, 3,  4, ie, ae, true);
			}
		}
		complex_t C8H_0 = 1./3.*lu*lu*F8_1(yt) - lu*ld*F8_2(yt);

        dep_param->set_expected(C8SMeps_0 +  C8Heps_0 + C8Heps2_0 + C8charg_0 + C8_chargeps_0 + C8H_0);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 6421, 0, 1)}, sources, func);

// 	complex_t C8SMeps_0= ((*sus_param).epsilonb-(*sus_param).epsilonbp)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))*(*susy)("HMIX",2)*F8_2((*sus_param).xt);
// 	complex_t C8Heps_0=(-(*sus_param).epsilon0p-(*sus_param).epsilonb)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))*(*susy)("HMIX",2)*F8_2((*sus_param).yt);
// 	complex_t C8Heps2_0=0.;

// 	if(((*sus_param).mass_A02==0.)&&((*sus_param).mass_H03==0.)) {

// 		C8Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("HMIX",2),2.)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))*F8_2((*sus_param).yt);
// 		C8Heps2_0+=(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("HMIX",2),2.)*(1.+pow((*susy)("HMIX",2),2.))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))/72.		*((cos((*susy)("ALPHA",0))+sin((*susy)("ALPHA",0))*(*susy)("HMIX",2))*(-sin((*susy)("ALPHA",0))+(*sus_param).epsilonb*cos((*susy)("ALPHA",0)))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",25),2.)
// 		+(sin((*susy)("ALPHA",0))-cos((*susy)("ALPHA",0))*(*susy)("HMIX",2))*(cos((*susy)("ALPHA",0))+(*sus_param).epsilonb*sin((*susy)("ALPHA",0)))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",35),2.)			+(-cos(atan((*susy)("HMIX",2)))-sin(atan((*susy)("HMIX",2)))*(*susy)("HMIX",2))*(sin(atan((*susy)("HMIX",2)))-(*sus_param).epsilonb*cos(atan((*susy)("HMIX",2))))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",36),2.));
// 	}
// 	else {	
//         C8Heps2_0=-(*sus_param).epsilon2*(*sus_param).epsilon1p*pow((*susy)("HMIX",2),2.)/(1.+(*sus_param).epsilonb*(*susy)("HMIX",2))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))*F8_2((*sus_param).yt);
// 		C8Heps2_0+=-3.*(*sus_param).epsilon2/pow(1.+(*sus_param).epsilonb*(*susy)("HMIX",2),2.)*(1.+pow((*susy)("HMIX",2),2.))/(1.+(*sus_param).epsilon0*(*susy)("HMIX",2))/72.
// 		*(((*susy)("HMIX", 00)+(*susy)("HMIX", 01)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 01)+(*sus_param).epsilonb*(*susy)("HMIX", 00))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",25),2.)
// 		+((*susy)("HMIX", 10)+(*susy)("HMIX", 11)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 11)+(*sus_param).epsilonb*(*susy)("HMIX", 10))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",35),2.)
// 		+((*susy)("HMIX", 31)+(*susy)("HMIX", 21)*(*susy)("HMIX",2))*(-(*susy)("HMIX", 21)+(*sus_param).epsilonb*(*susy)("HMIX", 20))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).mass_H03,2.)		

// +((*susy)("AMIX", 00)+(*susy)("AMIX", 01)*(*susy)("HMIX",2))*(-(*susy)("AMIX", 01)+(*sus_param).epsilonb*(*susy)("AMIX", 00))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*susy)("MASS",36),2.)
// 		+((*susy)("AMIX", 10)+(*susy)("AMIX", 11)*(*susy)("HMIX",2))*(-(*susy)("AMIX", 11)+(*sus_param).epsilonb*(*susy)("AMIX", 10))*pow(wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).mass_A02,2.));
// 		}

// 	complex_t C8charg_0 = 0.0;
// 	complex_t C8_chargeps_0 = 0.0;

// 	for (int ie = 0; ie < 2; ++ie) {
// 		for (int ae = 0; ae < 6; ++ae) {
// 			C8charg_0 += calculateContribution(h50, (*sus_param).X_UL, (*sus_param).X_UL, ie, ae, false) + (*sus_param).Mch[ie]/wilson_p("WPARAM_MATCH_SM", {5,1})*calculateContribution(h60,(*sus_param).X_UL,  (*sus_param).X_UR, ie, ae, false);
//             C8_chargeps_0 += (*sus_param).Mch[ie]/wilson_p("WPARAM_MATCH_SM", {5,1}) *calculateContribution(h60, (*sus_param).X_UL,  (*sus_param).X_UR, ie, ae, true);
// 		}
// 	}
//     complex_t C8H_0 = 1./3.*sus_param->lu*sus_param->lu*F8_1(sus_param->yt) - sus_param->lu*sus_param->ld*F8_2(sus_param->yt);
//     this->set_WilsonCoeffMatching("LO",C8SMeps_0 +  C8Heps_0 + C8Heps2_0 + C8charg_0 + C8_chargeps_0 + C8H_0);
//     // return C8SMeps_0 +  C8Heps_0 + C8Heps2_0 + C8charg_0 + C8_chargeps_0 + C8H_0;
}

void C8_susy::NLO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},     // lu
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},     // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},  // yt
		{ParameterType::WILSON, "EW_SCALE", 1},          // Q_match
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}}, // xt
		{ParameterType::SM, "MASS", 24},                 // mW
		{ParameterType::BSM, "MASS", 37},                // mH
		{ParameterType::BSM, "HMIX", 2},                 // tanb
		{ParameterType::BSM, "ALPHA", 0},                // alpha
		{ParameterType::BSM, "MASS", 35},                // mH0
		{ParameterType::BSM, "MASS", 25},                // mh
		{ParameterType::BSM, "MASS", 36},                // mA0
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}}, // mass_b_muW
	
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,3,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,4,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,5,5}}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
		double alpha = src.at({ParameterType::BSM, "ALPHA", 0})->get_val();
		double mH0 = src.at({ParameterType::BSM, "MASS", 35})->get_val();
		double mh = src.at({ParameterType::BSM, "MASS", 25})->get_val();
		double mA0 = src.at({ParameterType::BSM, "MASS", 36})->get_val();

		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();

		complex_t C8charg_1{};
		complex_t C8four_1{};

		for (int ie = 0; ie < 2; ie++) {
			for (int ae = 0; ae < 6; ae++) {
				C8charg_1+=pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*(src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() *src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() 
							*h51(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),log(pow(Q_match/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.))) 
							+ src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/mass_b_muW*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() *src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val() 
							*h61(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),log(pow(Q_match/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.))));

				for (int ce = 0; ce < 6; ce++) {
					double log_mu_W_MsqU_ce = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2.0));
					double MsqU_be_Mch_ie_ratio = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
					
					for (int de = 0; de < 6; de++) {
						double ratio_MsqU_ae_Mch_ie = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
						double ratio_MsqU_de_Mch_ie = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
						
						C8four_1 += std::pow(mW / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0) * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce}})->get_val() 
									* MsqU_be_Mch_ie_ratio * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ce, de}})->get_val() * (1.0 + log_mu_W_MsqU_ce) *
									(src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val()  * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 2}})->get_val()  * q21(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) +
									src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val() / mass_b_muW * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() 
									 * src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, de, 2}})->get_val()  * q41(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie));
					}
				}
			}
		}
		complex_t C8H_1 = G8H(yt,lu,ld)+Delta8H(yt,lu,ld)*log(pow(Q_match/mH,2.))
		-1./6.*EH(yt,lu);

        dep_param->set_expected(C8charg_1 + C8four_1+C8H_1);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 6421, 1, 1)}, sources, func);


    // complex_t C8charg_1{};
    // complex_t C8four_1{};

    // for (int ie = 0; ie < 2; ie++) {
	// 	for (int ae = 0; ae < 6; ae++) {
	// 		C8charg_1+=pow(sm("MASS", 24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h51(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(this->get_Q_match()/(*sus_param).MsqU[ae],2.))) + (*sus_param).Mch[ie]/sus_param->mass_b_muW*(*sus_param).X_UL[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h61(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),log(pow(this->get_Q_match()/(*sus_param).MsqU[ae],2.))));

	// 		for (int ce = 0; ce < 6; ce++) {
	// 			double log_mu_W_MsqU_ce = std::log(std::pow(this->get_Q_match() / (*sus_param).MsqU[ce], 2.0));
	// 			double MsqU_be_Mch_ie_ratio = (*sus_param).MsqU[ce] / (*sus_param).Mch[ie];
				
	// 			for (int de = 0; de < 6; de++) {
	// 				double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
	// 				double ratio_MsqU_de_Mch_ie = std::pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0);
					
	// 				C8four_1 += std::pow(sm("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).P_U[ae][ce] * MsqU_be_Mch_ie_ratio * (*sus_param).P_U[ce][de] * (1.0 + log_mu_W_MsqU_ce) *
	// 							((*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][de][2] * q21(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie) +
	// 							(*sus_param).Mch[ie] / wilson_p("WPARAM_MATCH_SM", {5,1}) * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UR[ie][de][2] * q41(ratio_MsqU_ae_Mch_ie, ratio_MsqU_de_Mch_ie));
	// 			}
	// 		}
	// 	}
	// }
    // complex_t C8H_1 = G8H(sus_param->yt,sus_param->lu,sus_param->ld)+Delta8H(sus_param->yt,sus_param->lu,sus_param->ld)*log(pow(this->get_Q_match()/sus_param->m_H,2.))
    // -1./6.*EH(sus_param->yt,sus_param->lu);

    // this->set_WilsonCoeffMatching("NLO", C8charg_1 + C8four_1+C8H_1);
    // // return C8charg_1 + C8four_1 + C8H_1;
}

void C8_susy::NNLO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},    // lu
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},    // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, // yt
		{ParameterType::WILSON, "EW_SCALE", 1},         // Q_match
		{ParameterType::WILSON, "WPARAM_MATCH_SM", 6}   // mass_top_muW
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();

		double coeff_temp =C8H2(yt,lu,ld,log(pow(Q_match/mass_top_muW, 2.)));

        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 6421, 2, 1)}, sources, func);

	// double coeff_temp =C8H2(sus_param->yt,sus_param->lu,sus_param->ld,log(pow(this->get_Q_match()/sus_param->mass_top_muW, 2.)));
    // // return this->double_to_complex_save("NNLO", coeff_temp);
}

void C9_susy::LO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},        // lu
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},     // yt
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},  // xt
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},         // sw2
		{ParameterType::WILSON, "MATRIX_BSM", 13},          // B90c
		{ParameterType::WILSON, "MATRIX_BSM", 14},          // C90c
		{ParameterType::WILSON, "MATRIX_BSM", 15}           // D90c
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {

		
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();

		complex_t C90c = src.at({ParameterType::WILSON, "MATRIX_BSM", 14})->get_val();
		complex_t B90c = src.at({ParameterType::WILSON, "MATRIX_BSM", 13})->get_val();
		complex_t D90c = src.at({ParameterType::WILSON, "MATRIX_BSM", 15})->get_val();
		complex_t C9charg_0 = (1.0 - 4.0 * sw2) / sw2 * C90c - B90c / sw2 - D90c;
		double C9H_0 = (1.-4.*sw2)/sw2*C9llH0(xt,yt,lu)-D9H0(yt,lu);

        dep_param->set_expected(C9charg_0+C9H_0);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4133, 0, 1)}, sources, func);

    // complex_t C9charg_0 = (1.0 - 4.0 * (*sus_param).sw2) / (*sus_param).sw2 * (*sus_param).C90c - (*sus_param).B90c / (*sus_param).sw2 - (*sus_param).D90c;
    // double C9H_0 = (1.-4.*sus_param->sw2)/sus_param->sw2*C9llH0(sus_param->xt,sus_param->yt,sus_param->lu)-D9H0(sus_param->yt,sus_param->lu);
    // this->set_WilsonCoeffMatching("LO",C9charg_0+C9H_0);
    // // return C9charg_0+C9H_0;
}

void C9_susy::NLO_calculation() {

	std::unordered_set<ParamId> sources {
		// Paramètres simples
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},   // kappa
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},   // lu
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},   // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, // yt
		{ParameterType::WILSON, "EW_SCALE", 1},         // Q_match
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}}, // xt
		{ParameterType::SM, "MASS", 24},                // mW
		{ParameterType::BSM, "MASS", 37},               // mH
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},     // sw2
		{ParameterType::SM, "GAUGE", 2},                // g2
	
		// Pour les boucles : masses
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
	
		// Matrices utilisées
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,5}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,3}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,4}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,5}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
		
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,1}},
		
		// Mixages
		{ParameterType::BSM, "UMIX", 0},
		{ParameterType::BSM, "UMIX", 10},
		{ParameterType::BSM, "VMIX", 0},
		{ParameterType::BSM, "VMIX", 10}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double H03  = 0.; //TODO, wtf
		double A02 = 0.; //TODO, WTF
		complex_t C91f = 0.0;
		complex_t B1f1 = 0.0;
		complex_t B1f2 = 0.0;
		complex_t D91f = 0.0;
		complex_t B1c1=0.;
		complex_t B1c2=0.;
		complex_t D91c= 0.;
		complex_t C91c = 0.;

		for (int ie = 0; ie < 2; ie++) {

			for (int ae = 0; ae < 6; ae++) {
				double mass24_Mch_ie_squared = pow(mW / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double ratio_MsqU_ae_Mch_ie = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double log_mu_W_MsqU_ae = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2.0));

				D91c += std::pow(mW / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0) * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() *
					h31(ratio_MsqU_ae_Mch_ie, log_mu_W_MsqU_ae);

				

				for (int je = 0; je < 2; je++) {
					double ratio_Mch_je_ie = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
					double factor_abs = 2.0 * std::fabs(ratio_Mch_je_ie);
					double factor_f31_f30 = (f31(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
											4.0 * (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
											(f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
											f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);
					double factor_f41_f40 = (f41(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
											4.0 * (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
											(f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
											f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);

					C91c += src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() *
							((factor_abs * factor_f31_f30 * src.at({ParameterType::BSM, "UMIX",je*10+0})->get_val() * src.at({ParameterType::BSM, "UMIX",ie*10+0})->get_val()) -
							(factor_f41_f40 * src.at({ParameterType::BSM, "VMIX",je*10+0})->get_val() * src.at({ParameterType::BSM, "VMIX",ie*10+0})->get_val()));

					for (int de=0; de<6; de++) {
						for (int ke=0; ke<6; ke++) {
							C91f+= src.at({ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}})->get_val() * pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.) *  src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}})->get_val() * (1+log(pow(Q_match/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(),2.)))*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * (
								2.*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()) * f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.), pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.), pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.))*src.at({ParameterType::BSM, "UMIX",je*10+0})->get_val()*src.at({ParameterType::BSM, "UMIX",ie*10+0})->get_val()
								-f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.), pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.), pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.)) *  src.at({ParameterType::BSM, "VMIX",je*10+0})->get_val()*src.at({ParameterType::BSM, "VMIX",ie*10+0})->get_val());
							
						}
					}

					for (int be=0; be<3; be++) {
						double ratio_Mch = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
						double ratio_MsqU = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
						double ratio_Msn = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
						
						B1c1 += src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() / pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2) * 
								(0.5 * src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}})->get_val() * 
								(f81(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
								4 * (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
								(f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
								f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
								log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2))));

						B1c2 += src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() / pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2) * 
								(src.at({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}})->get_val() * fabs(ratio_Mch) * 
								(f91(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
								4 * (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
								(f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
								f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
								log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2))));
					}
				}
				for (int ce = 0; ce < 6; ce++) {
					for (int fe = 0; fe < 3; fe++) {
						for (int de = 0; de < 6; de++) {
							for (int ke = 0; ke < 6; ke++) {
								double MsqU_ke_Mch_ie_squared = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
								double log_scale_MsqU_ke = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(), 2.0));
								C91f +=  src.at({ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}})->get_val() * MsqU_ke_Mch_ie_squared *  src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}})->get_val() * 
										(1.0 + log_scale_MsqU_ke) * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * 
										f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) * 
										src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, fe}})->get_val();

								C91f +=  src.at({ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}})->get_val() * MsqU_ke_Mch_ie_squared *  src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ke, ce}})->get_val() * 
										(1.0 + log_scale_MsqU_ke) * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * 
										f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) * 
										src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val();
								
							}

							for (int je = 0; je < 2; je++) {
								double factor_common = mass24_Mch_ie_squared * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, de}})->get_val() *pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2)*src.at({ParameterType::WILSON, "MATRIX_BSM", {9,de, ce}})->get_val() *
								(1+log(pow(Q_match/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val(),2))) *  src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 2}})->get_val();

								B1f1 += factor_common * 0.5 * f90(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) *
										src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, fe, 1}})->get_val();

								B1f2 += factor_common * fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()) * 
										f100(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
												pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
												pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
												pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) *
										src.at({ParameterType::WILSON, "MATRIX_BSM", {6, ie, fe, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, fe, 1}})->get_val();
							}
						}
						C91c += src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() *
							(f51(ratio_MsqU_ae_Mch_ie, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) +
								4.0 * (f40(ratio_MsqU_ae_Mch_ie, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) + 
								(f40(ratio_MsqU_ae_Mch_ie, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)* 1.0001) - 
								f40(ratio_MsqU_ae_Mch_ie, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)* 0.9999)) / 0.0002
								+(f40(ratio_MsqU_ae_Mch_ie*1.0001, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))
								-f40(ratio_MsqU_ae_Mch_ie*0.9999, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)))/0.0002)* log_mu_W_MsqU_ae) *
							src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val();
					}

					double MsqU_ce_Mch_ie_squared = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
					double log_scale_MsqU_ce = 1.0 + log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val(), 2.0));

					
					for (int de = 0; de < 6; de++) {
						double ratio_MsqU_ae_Mch_ie = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
						double ratio_MsqU_de_Mch_ie = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
						
						D91f += mass24_Mch_ie_squared * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce}})->get_val() * MsqU_ce_Mch_ie_squared * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ce, de}})->get_val() * log_scale_MsqU_ce * 
							src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 2}})->get_val() * 
							q51(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), ratio_MsqU_de_Mch_ie);

					}
				}
			}
		}


		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		C91c *= -kappa / 8.0;
		D91c *= kappa;
		complex_t B91c = -(B1c1 - B1c2) * kappa * mW * mW / 2.0 / pow(g2, 2);

		C91f *= kappa / 6.0;
		D91f *= kappa;

		complex_t B91f = (B1f1 - B1f2) * 2.0 / 3.0 * kappa / pow(g2, 2);

		complex_t C9four_1 = (1. - 4. * sw2) / sw2 * C91f - B91f / sw2 - D91f;

		complex_t C9charg_1=(1.-4.*sw2)/sw2*C91c-B91c/sw2-D91c;
		complex_t C9H_1 = (1.-4.*sw2)/sw2*C9llH1(xt,yt,lu,log(pow(Q_match/mH,2.)))
		-D9H1(yt,lu,log(pow(Q_match/mH,2.)));

        dep_param->set_expected(C9four_1+C9charg_1 + C9H_1);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4133, 1, 1)}, sources, func);


	// complex_t C91f = 0.0;
	// complex_t B1f1 = 0.0;
	// complex_t B1f2 = 0.0;
	// complex_t D91f = 0.0;
	// complex_t B1c1=0.;
	// complex_t B1c2=0.;
	// complex_t D91c= 0.;
	// complex_t C91c = 0.;

	// for (int ie = 0; ie < 2; ie++) {

	// 	for (int ae = 0; ae < 6; ae++) {
	// 		double mass24_Mch_ie_squared = pow(sm("MASS", 24) / (*sus_param).Mch[ie], 2.0);
	// 		double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
	// 		double log_mu_W_MsqU_ae = std::log(std::pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2.0));

	// 		D91c += std::pow(sm("MASS", 24) / (*sus_param).Mch[ie], 2.0) * (*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][ae][2] *
    //             h31(ratio_MsqU_ae_Mch_ie, log_mu_W_MsqU_ae);

			

	// 		for (int je = 0; je < 2; je++) {
	// 			double ratio_Mch_je_ie = (*sus_param).Mch[je] / (*sus_param).Mch[ie];
	// 			double factor_abs = 2.0 * std::fabs(ratio_Mch_je_ie);
	// 			double factor_f31_f30 = (f31(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
	// 									4.0 * (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
	// 									(f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
	// 									f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);
	// 			double factor_f41_f40 = (f41(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
	// 									4.0 * (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
	// 									(f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
	// 									f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);

	// 			C91c += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] *
	// 					((factor_abs * factor_f31_f30 * (*susy)("UMIX", je*10+0) * (*susy)("UMIX", ie*10+0)) -
	// 					(factor_f41_f40 * (*susy)("VMIX", je*10+0) * (*susy)("VMIX", ie*10+0)));

	// 			for (int de=0; de<6; de++) {
	// 				for (int ke=0; ke<6; ke++) {
	// 					C91f+=(*sus_param).P_U[de][ke] * pow((*sus_param).MsqU[ke]/(*sus_param).Mch[ie],2.) * (*sus_param).P_U[ke][ae] * (1+log(pow(this->get_Q_match()/(*sus_param).MsqU[ke],2.)))*(*sus_param).X_UL[je][de][1]*(*sus_param).X_UL[ie][ae][2] * (
	// 						2.*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie]) * f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.), pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie], 2.))*(*susy)("UMIX", je*10+0)*(*susy)("UMIX", ie*10+0)
	// 						-f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie], 2.)) *  (*susy)("VMIX", je * 10+0)*(*susy)("VMIX", ie * 10+0));
						
	// 				}
	// 			}

	// 			for (int be=0; be<3; be++) {
	// 				double ratio_Mch = (*sus_param).Mch[je] / (*sus_param).Mch[ie];
	// 				double ratio_MsqU = (*sus_param).MsqU[ae] / (*sus_param).Mch[ie];
	// 				double ratio_Msn = (*sus_param).Msn[be] / (*sus_param).Mch[ie];
					
	// 				B1c1 += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] / pow((*sus_param).Mch[ie], 2) * 
	// 						(0.5 * (*sus_param).X_NL[ie][be][1] * (*sus_param).X_NL[je][be][1] * 
	// 						(f81(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
	// 						4 * (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
	// 						(f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
	// 						f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
	// 						log(pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2))));

	// 				B1c2 += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] / pow((*sus_param).Mch[ie], 2) * 
	// 						((*sus_param).X_NR[ie][be][1] * (*sus_param).X_NR[je][be][1] * fabs(ratio_Mch) * 
	// 						(f91(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
	// 						4 * (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
	// 						(f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
	// 						f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
	// 						log(pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2))));
	// 			}
	// 		}
	// 		for (int ce = 0; ce < 6; ce++) {
	// 			for (int fe = 0; fe < 3; fe++) {
	// 				for (int de = 0; de < 6; de++) {
	// 					for (int ke = 0; ke < 6; ke++) {
	// 						double MsqU_ke_Mch_ie_squared = pow((*sus_param).MsqU[ke] / (*sus_param).Mch[ie], 2.0);
	// 						double log_scale_MsqU_ke = log(pow(this->get_Q_match() / (*sus_param).MsqU[ke], 2.0));
	// 						C91f += (*sus_param).P_U[de][ke] * MsqU_ke_Mch_ie_squared * (*sus_param).P_U[ke][ae] * 
	// 								(1.0 + log_scale_MsqU_ke) * (*sus_param).X_UL[ie][ce][1] * (*sus_param).X_UL[ie][ae][2] * 
	// 								f50(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
	// 									pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0), 
	// 									pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) * 
	// 								(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[de][fe];

	// 						C91f += (*sus_param).P_U[de][ke] * MsqU_ke_Mch_ie_squared * (*sus_param).P_U[ke][ce] * 
	// 								(1.0 + log_scale_MsqU_ke) * (*sus_param).X_UL[ie][de][1] * (*sus_param).X_UL[ie][ae][2] * 
	// 								f50(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
	// 									pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
	// 									pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0)) * 
	// 								(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[ae][fe];
							
	// 					}

	// 					for (int je = 0; je < 2; je++) {
	// 						double factor_common = mass24_Mch_ie_squared * (*sus_param).P_U[ae][de] *pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2)*(*sus_param).P_U[de][ce] *
	// 						(1+log(pow(this->get_Q_match()/(*sus_param).MsqU[de],2))) *  (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ce][2];

	// 						B1f1 += factor_common * 0.5 * f90(pow((*sus_param).Mch[je] / (*sus_param).Mch[ie], 2.0), 
	// 									pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
	// 									pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
	// 									pow((*sus_param).Msn[fe] / (*sus_param).Mch[ie], 2.0)) *
	// 								(*sus_param).X_NL[ie][fe][1] * (*sus_param).X_NL[je][fe][1];

	// 						B1f2 += factor_common * fabs((*sus_param).Mch[je] / (*sus_param).Mch[ie]) * 
	// 								f100(pow((*sus_param).Mch[je] / (*sus_param).Mch[ie], 2.0), 
	// 										pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
	// 										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
	// 										pow((*sus_param).Msn[fe] / (*sus_param).Mch[ie], 2.0)) *
	// 								(*sus_param).X_NR[ie][fe][1] * (*sus_param).X_NR[je][fe][1];
	// 					}
	// 				}
	// 				C91c += (*sus_param).X_UL[ie][ce][1] * (*sus_param).X_UL[ie][ae][2] *
	// 					(f51(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) +
	// 						4.0 * (f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) + 
	// 						(f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)* 1.0001) - 
	// 						f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)* 0.9999)) / 0.0002
	// 						+(f40(ratio_MsqU_ae_Mch_ie*1.0001, std::pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.))
	// 						-f40(ratio_MsqU_ae_Mch_ie*0.9999, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)))/0.0002)* log_mu_W_MsqU_ae) *
	// 					(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[ae][fe];
	// 			}

	// 			double MsqU_ce_Mch_ie_squared = pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0);
	// 			double log_scale_MsqU_ce = 1.0 + log(pow(this->get_Q_match() / (*sus_param).MsqU[ce], 2.0));

				
	// 			for (int de = 0; de < 6; de++) {
	// 				double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
	// 				double ratio_MsqU_de_Mch_ie = std::pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0);
					
	// 				D91f += mass24_Mch_ie_squared * (*sus_param).P_U[ae][ce] * MsqU_ce_Mch_ie_squared * (*sus_param).P_U[ce][de] * log_scale_MsqU_ce * 
	// 					(*sus_param).X_UL[ie][ae][1] * (*sus_param).X_UL[ie][de][2] * 
	// 					q51(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), ratio_MsqU_de_Mch_ie);

	// 			}
	// 		}
	// 	}
	// }



	// C91c *= -(*sus_param).kappa / 8.0;
	// D91c *= (*sus_param).kappa;
	// complex_t B91c = -(B1c1 - B1c2) * (*sus_param).kappa * sm("MASS", 24) * sm("MASS", 24) / 2.0 / pow(sm("GAUGE", 2), 2);

	// C91f *= (*sus_param).kappa / 6.0;
	// D91f *= (*sus_param).kappa;

	// complex_t B91f = (B1f1 - B1f2) * 2.0 / 3.0 * (*sus_param).kappa / pow(sm("GAUGE", 2), 2);

    // complex_t C9four_1 = (1. - 4. * (*sus_param).sw2) / (*sus_param).sw2 * C91f - B91f / (*sus_param).sw2 - D91f;

	// complex_t C9charg_1=(1.-4.*(*sus_param).sw2)/(*sus_param).sw2*C91c-B91c/(*sus_param).sw2-D91c;
    // complex_t C9H_1 = (1.-4.*sus_param->sw2)/sus_param->sw2*C9llH1(sus_param->xt,sus_param->yt,sus_param->lu,log(pow(this->get_Q_match()/sus_param->m_H,2.)))
    // -D9H1(sus_param->yt,sus_param->lu,log(pow(this->get_Q_match()/sus_param->m_H,2.)));
    // this->set_WilsonCoeffMatching("NLO",C9four_1+C9charg_1 + C9H_1);
    // // return C9four_1+C9charg_1 + C9H_1;
}

void C10_susy::LO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},        // lu
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},     // yt
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},  // xt
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},         // sw2
		{ParameterType::WILSON, "MATRIX_BSM", 14},          // C90c
		{ParameterType::WILSON, "MATRIX_BSM", 16}           // B100c
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {

		
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();
		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();

		complex_t C90c = src.at({ParameterType::WILSON, "MATRIX_BSM", 14})->get_val();
		complex_t B100c = src.at({ParameterType::WILSON, "MATRIX_BSM", 16})->get_val();
		complex_t C10charg_0 = (B100c - C90c) / sw2;
		complex_t C10H_0 = -C9llH0(xt,yt,lu)/sw2;

        dep_param->set_expected(C10charg_0+C10H_0);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4137, 0, 1)}, sources, func);

    // complex_t C10charg_0 = ((*sus_param).B100c - (*sus_param).C90c) / (*sus_param).sw2;
    // complex_t C10H_0 = -C9llH0(sus_param->xt,sus_param->yt,sus_param->lu)/sus_param->sw2;
    // this->set_WilsonCoeffMatching("LO",C10charg_0+C10H_0);
    // // return C10charg_0+C10H_0;
}

void C10_susy::NLO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
		{ParameterType::WILSON, "EW_SCALE", 1},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::SM, "GAUGE", 2},
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},
		{ParameterType::BSM, "MASS", 37},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,3,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,4,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,5,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {9,2,2}},
		{ParameterType::BSM, "UMIX", 0},
		{ParameterType::BSM, "UMIX", 10},
		{ParameterType::BSM, "VMIX", 0},
		{ParameterType::BSM, "VMIX", 10}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();

		complex_t C91f = 0.0;
		complex_t B1f1 = 0.0;
		complex_t B1f2 = 0.0;
		complex_t B1c1=0.;
		complex_t B1c2=0.;
		complex_t C91c = 0.;

		for (int ie = 0; ie < 2; ie++) {
			for (int ae = 0; ae < 6; ae++) {
				double mass24_Mch_ie_squared = pow(mW / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double ratio_MsqU_ae_Mch_ie = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
				double log_mu_W_MsqU_ae = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2.0));

				for (int je = 0; je < 2; je++) {
					double ratio_Mch_je_ie = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
					double factor_abs = 2.0 * std::fabs(ratio_Mch_je_ie);
					double factor_f31_f30 = (f31(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
											4.0 * (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
											(f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
											f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);
					double factor_f41_f40 = (f41(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
											4.0 * (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
											(f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
											f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);

					C91c += src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() *
							((factor_abs * factor_f31_f30 * src.at({ParameterType::BSM, "UMIX",je*10+0})->get_val() * src.at({ParameterType::BSM, "UMIX",ie*10+0})->get_val()) -
							(factor_f41_f40 * src.at({ParameterType::BSM, "VMIX",je*10+0})->get_val() * src.at({ParameterType::BSM, "VMIX",ie*10+0})->get_val()));

					for (int de=0; de<6; de++) {
						for (int ke=0; ke<6; ke++) {
							C91f+=src.at({ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}})->get_val() * pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.) * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}})->get_val() * (1+log(pow(Q_match/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(),2.)))*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, de, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * (
								2.*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()) * f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.), pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.), pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.))*src.at({ParameterType::BSM, "UMIX",je*10+0})->get_val()*src.at({ParameterType::BSM, "UMIX",ie*10+0})->get_val()
								-f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.), pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.), pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.)) *  src.at({ParameterType::BSM, "VMIX",je*10+0})->get_val()*src.at({ParameterType::BSM, "VMIX",ie*10+0})->get_val());
							
						}
					}

					for (int be=0; be<3; be++) {
						double ratio_Mch = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
						double ratio_MsqU = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
						double ratio_Msn = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val();
						
						B1c1 += src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() / pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2) * 
								(0.5 * src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}})->get_val() * 
								(f81(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
								4 * (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
								(f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
								f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
								log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2))));

						B1c2 += src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() / pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2) * 
								(src.at({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}})->get_val() * fabs(ratio_Mch) * 
								(f91(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
								4 * (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
								(f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
								f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
								log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(), 2))));
					}
				}

				for (int ce = 0; ce < 6; ce++) {
					for (int fe = 0; fe < 3; fe++) {
						for (int de = 0; de < 6; de++) {
							for (int ke = 0; ke < 6; ke++) {
								double MsqU_ke_Mch_ie_squared = pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0);
								double log_scale_MsqU_ke = log(pow(Q_match / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(), 2.0));
								C91f += src.at({ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}})->get_val() * MsqU_ke_Mch_ie_squared * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ke, ae}})->get_val() * 
										(1.0 + log_scale_MsqU_ke) * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * 
										f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) * 
										src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, fe}})->get_val();

								C91f += src.at({ParameterType::WILSON, "MATRIX_BSM", {9,de, ke}})->get_val() * MsqU_ke_Mch_ie_squared * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ke, ce}})->get_val() * 
										(1.0 + log_scale_MsqU_ke) * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, de, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() * 
										f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) * 
										src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val();
								
							}

							for (int je = 0; je < 2; je++) {
								// LOG_INFO(std::to_string(src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1}})->get_val()));
								double factor_common = mass24_Mch_ie_squared * src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, de}})->get_val() *pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2)*src.at({ParameterType::WILSON, "MATRIX_BSM", {9,de, ce}})->get_val() *
								(1+log(pow(Q_match/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val(),2))) *  src.at({ParameterType::WILSON, "MATRIX_BSM", {3,je, ae, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 2}})->get_val();

								B1f1 += factor_common * 0.5 * f90(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
											pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) *
										src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, fe, 1}})->get_val();

								B1f2 += factor_common * fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()) * 
										f100(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
												pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
												pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0), 
												pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) *
										src.at({ParameterType::WILSON, "MATRIX_BSM", {6, ie, fe, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, fe, 1}})->get_val();

								
							}
						}

						C91c += src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ce, 1}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val() *
							(f51(ratio_MsqU_ae_Mch_ie, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) +
								4.0 * (f40(ratio_MsqU_ae_Mch_ie, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)) + 
								(f40(ratio_MsqU_ae_Mch_ie, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)* 1.0001) - 
								f40(ratio_MsqU_ae_Mch_ie, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)* 0.9999)) / 0.0002
								+(f40(ratio_MsqU_ae_Mch_ie*1.0001, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))
								-f40(ratio_MsqU_ae_Mch_ie*0.9999, pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val() / src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(), 2.0)))/0.0002)* log_mu_W_MsqU_ae) *
							src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val() * src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val();
					}

				}
			}
		}


		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		C91c *= -kappa / 8.0;

		complex_t B101c = (B1c1 + B1c2) * kappa * mW * mW / 2.0 / pow(g2, 2);

		C91f *= kappa / 6.0;

		complex_t B101f = -(B1f1 + B1f2) * 2.0 / 3.0 * kappa / pow(g2, 2);

		complex_t C10four_1 = (B101f - C91f) / sw2;	

		complex_t C10charg_1=(B101c-C91c)/sw2;
		double C10H_1 = -C9llH1(xt,yt,lu,log(pow(Q_match/mH,2.)))/sw2;
		// this->set_WilsonCoeffMatching("NLO",C10four_1+C10charg_1+C10H_1);

        dep_param->set_expected(C10four_1+C10charg_1+C10H_1);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4137, 1, 1)}, sources, func);


// 	complex_t C91f = 0.0;
// 	complex_t B1f1 = 0.0;
// 	complex_t B1f2 = 0.0;
// 	complex_t B1c1=0.;
// 	complex_t B1c2=0.;
// 	complex_t C91c = 0.;

// 	for (int ie = 0; ie < 2; ie++) {
// 		for (int ae = 0; ae < 6; ae++) {
// 			double mass24_Mch_ie_squared = pow(sm("MASS", 24) / (*sus_param).Mch[ie], 2.0);
// 			double ratio_MsqU_ae_Mch_ie = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
// 			double log_mu_W_MsqU_ae = std::log(std::pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2.0));

// 			for (int je = 0; je < 2; je++) {
// 				double ratio_Mch_je_ie = (*sus_param).Mch[je] / (*sus_param).Mch[ie];
// 				double factor_abs = 2.0 * std::fabs(ratio_Mch_je_ie);
// 				double factor_f31_f30 = (f31(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
// 										4.0 * (f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
// 										(f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
// 										f30(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);
// 				double factor_f41_f40 = (f41(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) +
// 										4.0 * (f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie) + 
// 										(f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 1.0001) - 
// 										f40(std::pow(ratio_Mch_je_ie, 2), ratio_MsqU_ae_Mch_ie * 0.9999)) / 0.0002) * log_mu_W_MsqU_ae);

// 				C91c += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] *
// 						((factor_abs * factor_f31_f30 * (*susy)("UMIX", je*10+0) * (*susy)("UMIX", ie*10+0)) -
// 						(factor_f41_f40 * (*susy)("VMIX", je*10+0) * (*susy)("VMIX", ie*10+0)));

// 				for (int de=0; de<6; de++) {
// 					for (int ke=0; ke<6; ke++) {
// 						C91f+=(*sus_param).P_U[de][ke] * pow((*sus_param).MsqU[ke]/(*sus_param).Mch[ie],2.) * (*sus_param).P_U[ke][ae] * (1+log(pow(this->get_Q_match()/(*sus_param).MsqU[ke],2.)))*(*sus_param).X_UL[je][de][1]*(*sus_param).X_UL[ie][ae][2] * (
// 							2.*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie]) * f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.), pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie], 2.))*(*susy)("UMIX", je*10+0)*(*susy)("UMIX", ie*10+0)
// 							-f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie], 2.), pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie], 2.)) *  (*susy)("VMIX", je * 10+0)*(*susy)("VMIX", ie * 10+0));
						
// 					}
// 				}

// 				for (int be=0; be<3; be++) {
// 					double ratio_Mch = (*sus_param).Mch[je] / (*sus_param).Mch[ie];
// 					double ratio_MsqU = (*sus_param).MsqU[ae] / (*sus_param).Mch[ie];
// 					double ratio_Msn = (*sus_param).Msn[be] / (*sus_param).Mch[ie];
					
// 					B1c1 += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] / pow((*sus_param).Mch[ie], 2) * 
// 							(0.5 * (*sus_param).X_NL[ie][be][1] * (*sus_param).X_NL[je][be][1] * 
// 							(f81(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
// 							4 * (f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
// 							(f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
// 							f50(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
// 							log(pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2))));

// 					B1c2 += (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ae][2] / pow((*sus_param).Mch[ie], 2) * 
// 							((*sus_param).X_NR[ie][be][1] * (*sus_param).X_NR[je][be][1] * fabs(ratio_Mch) * 
// 							(f91(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) + 
// 							4 * (f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2), pow(ratio_Msn, 2)) +
// 							(f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 1.0001, pow(ratio_Msn, 2)) - 
// 							f60(pow(ratio_Mch, 2), pow(ratio_MsqU, 2) * 0.9999, pow(ratio_Msn, 2))) / 0.0002) *
// 							log(pow(this->get_Q_match() / (*sus_param).MsqU[ae], 2))));
// 				}
// 			}

// 			for (int ce = 0; ce < 6; ce++) {
// 				for (int fe = 0; fe < 3; fe++) {
// 					for (int de = 0; de < 6; de++) {
// 						for (int ke = 0; ke < 6; ke++) {
// 							double MsqU_ke_Mch_ie_squared = pow((*sus_param).MsqU[ke] / (*sus_param).Mch[ie], 2.0);
// 							double log_scale_MsqU_ke = log(pow(this->get_Q_match() / (*sus_param).MsqU[ke], 2.0));
// 							C91f += (*sus_param).P_U[de][ke] * MsqU_ke_Mch_ie_squared * (*sus_param).P_U[ke][ae] * 
// 									(1.0 + log_scale_MsqU_ke) * (*sus_param).X_UL[ie][ce][1] * (*sus_param).X_UL[ie][ae][2] * 
// 									f50(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
// 										pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0), 
// 										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) * 
// 									(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[de][fe];

// 							C91f += (*sus_param).P_U[de][ke] * MsqU_ke_Mch_ie_squared * (*sus_param).P_U[ke][ce] * 
// 									(1.0 + log_scale_MsqU_ke) * (*sus_param).X_UL[ie][de][1] * (*sus_param).X_UL[ie][ae][2] * 
// 									f50(pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
// 										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
// 										pow((*sus_param).MsqU[de] / (*sus_param).Mch[ie], 2.0)) * 
// 									(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[ae][fe];
							
// 						}

// 						for (int je = 0; je < 2; je++) {
// 							// LOG_INFO(std::to_string((*sus_param).X_NL[ie][fe][1]));
// 							double factor_common = mass24_Mch_ie_squared * (*sus_param).P_U[ae][de] *pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2)*(*sus_param).P_U[de][ce] *
// 							(1+log(pow(this->get_Q_match()/(*sus_param).MsqU[de],2))) *  (*sus_param).X_UL[je][ae][1] * (*sus_param).X_UL[ie][ce][2];

// 							B1f1 += factor_common * 0.5 * f90(pow((*sus_param).Mch[je] / (*sus_param).Mch[ie], 2.0), 
// 										pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
// 										pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
// 										pow((*sus_param).Msn[fe] / (*sus_param).Mch[ie], 2.0)) *
// 									(*sus_param).X_NL[ie][fe][1] * (*sus_param).X_NL[je][fe][1];

// 							B1f2 += factor_common * fabs((*sus_param).Mch[je] / (*sus_param).Mch[ie]) * 
// 									f100(pow((*sus_param).Mch[je] / (*sus_param).Mch[ie], 2.0), 
// 											pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0), 
// 											pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0), 
// 											pow((*sus_param).Msn[fe] / (*sus_param).Mch[ie], 2.0)) *
// 									(*sus_param).X_NR[ie][fe][1] * (*sus_param).X_NR[je][fe][1];

							
// 						}
// 					}

// 					C91c += (*sus_param).X_UL[ie][ce][1] * (*sus_param).X_UL[ie][ae][2] *
// 						(f51(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) +
// 							4.0 * (f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)) + 
// 							(f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)* 1.0001) - 
// 							f40(ratio_MsqU_ae_Mch_ie, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)* 0.9999)) / 0.0002
// 							+(f40(ratio_MsqU_ae_Mch_ie*1.0001, std::pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.))
// 							-f40(ratio_MsqU_ae_Mch_ie*0.9999, std::pow((*sus_param).MsqU[ce] / (*sus_param).Mch[ie], 2.0)))/0.0002)* log_mu_W_MsqU_ae) *
// 						(*sus_param).Gamma_UL[ce][fe] * (*sus_param).Gamma_UL[ae][fe];
// 				}

// 			}
// 		}
	// }



// 	C91c *= -(*sus_param).kappa / 8.0;

// 	complex_t B101c = (B1c1 + B1c2) * (*sus_param).kappa * sm("MASS", 24) * sm("MASS", 24) / 2.0 / pow(sm("GAUGE", 2), 2);

// 	C91f *= (*sus_param).kappa / 6.0;

// 	complex_t B101f = -(B1f1 + B1f2) * 2.0 / 3.0 * (*sus_param).kappa / pow(sm("GAUGE", 2), 2);

//     complex_t C10four_1 = (B101f - C91f) / (*sus_param).sw2;	

// 	complex_t C10charg_1=(B101c-C91c)/(*sus_param).sw2;
//     double C10H_1 = -C9llH1(sus_param->xt,sus_param->yt,sus_param->lu,log(pow(this->get_Q_match()/sus_param->m_H,2.)))/sus_param->sw2;
//     this->set_WilsonCoeffMatching("NLO",C10four_1+C10charg_1+C10H_1);
//     // return C10four_1+C10charg_1+C10H_1;
}

void CP7_susy::LO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},        // kappa
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},        // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},     // yt
		{ParameterType::WILSON, "WPARAM_MATCH_SM", 6},      // mass_top_muW
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},  // mass_b_muW
		{ParameterType::SM, "MASS", 24},                    // mW
		{ParameterType::SM, "MASS", 3},                     // masse pour C7pH
	
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {

		
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double C7pH=src.at({ParameterType::SM, "MASS", 3})->get_val()*mass_b_muW/mass_top_muW/mass_top_muW*1./3.*ld*ld*F7_1(yt);
		

		double C7pcharg=0.;
		

		double B10pc=0.;
		double C9pc=0.;
		double B9pc=0.;
		double D9pc=0.; 

		double BQ1pc1=0.;
		double BQ1pc2=0.;

		double Dp{0},Dm{0};
		double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a1{0}, NQ1pc{0}, NQ2pc{0};

		for(int ie=0;ie<2;ie++) {
			for(int ae=0;ae<6;ae++) {
				C7pcharg+=pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*(src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1}})->get_val()
				*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()*h10(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)) 
				+ src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/mass_b_muW*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()
				*h20(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
			}
		} 
		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		C7pcharg*=-0.5*kappa;
		// return this->double_to_complex_save("LO", C7pH+C7pcharg);

        dep_param->set_expected(C7pH+C7pcharg);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 4322, 0, 1)}, sources, func);

	// // sus_param->reset_PrimeCQG(this->get_Q_match());
	// double C7pH=sm("MASS",3)*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW*1./3.*(*sus_param).ld*(*sus_param).ld*F7_1((*sus_param).yt);
	

	// double C7pcharg=0.;
	

	// double B10pc=0.;
	// double C9pc=0.;
	// double B9pc=0.;
	// double D9pc=0.; 

	// double BQ1pc1=0.;
	// double BQ1pc2=0.;

	// double Dp{0},Dm{0};
	// double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a1{0}, NQ1pc{0}, NQ2pc{0};

	// for(int ie=0;ie<2;ie++) {
	// 	for(int ae=0;ae<6;ae++) {
	// 		C7pcharg+=pow(sm("MASS",24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h10(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)) + (*sus_param).Mch[ie]/wilson_p("WPARAM_MATCH_SM", {5,1})*(*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h20(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
	// 	}
	// } 		
	// C7pcharg*=-0.5*(*sus_param).kappa;
    // // return this->double_to_complex_save("LO", C7pH+C7pcharg);
}

void CP8_susy::LO_calculation() {

	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},         // kappa
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},         // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},      // yt
		{ParameterType::WILSON, "WPARAM_MATCH_SM", 6},       // mass_top_muW
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},   // mass_b_muW
		{ParameterType::SM, "MASS", 24},                     // mW
		{ParameterType::SM, "MASS", 3},                      // masse utilisée deux fois
	
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,5,2}}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {

		
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double C7pH=src.at({ParameterType::SM, "MASS", 3})->get_val()*mass_b_muW/mass_top_muW/mass_top_muW*1./3.*ld*ld*F7_1(yt);
		

		double C8pH=src.at({ParameterType::SM, "MASS", 3})->get_val()*mass_b_muW/mass_top_muW/mass_top_muW*1./3.*ld*ld*F8_1(yt);
		double C8pcharg=0.;
		for(int ie=0;ie<2;ie++) {
			for(int ae=0;ae<6;ae++) {
				C8pcharg+=pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*(src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()*h50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)) 
						+ src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/mass_b_muW*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()*h60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
			}
		}

		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		C8pcharg*=-0.5*kappa; 

		// return this->double_to_complex_save("LO", C7pH+C7pcharg);

        dep_param->set_expected(C8pH+C8pcharg);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(305, 4321, 0, 1)}, sources, func);

    // double C8pH=sm("MASS",3)*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW*1./3.*(*sus_param).ld*(*sus_param).ld*F8_1((*sus_param).yt);
    // double C8pcharg=0.;
	// for(int ie=0;ie<2;ie++) {
	// 	for(int ae=0;ae<6;ae++) {
	// 		C8pcharg+=pow(sm("MASS",24)/(*sus_param).Mch[ie],2.)*((*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h50(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)) + (*sus_param).Mch[ie]/wilson_p("WPARAM_MATCH_SM", {5,1})*(*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UL[ie][ae][2]*h60(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)));
    //     }
    // }
    // C8pcharg*=-0.5*(*sus_param).kappa;
    // // return this->double_to_complex_save("LO", C8pH+C8pcharg);
}

void CP9_susy::LO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},   // kappa
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},   // lu
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},   // ld
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},   // yt
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},  // xt
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},  // mass_b_muW
		{ParameterType::WILSON, "WPARAM_MATCH_SM", 6},  // mass_top_muW
		{ParameterType::WILSON, "WPARAM_SI_SM", 3},     // pour C10pH
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},     // sw2
		{ParameterType::SM, "MASS", 3},                 // masse dans C9pH et C10pH
		{ParameterType::SM, "MASS", 24},                // mW
		{ParameterType::SM, "GAUGE", 2},                // g2
		{ParameterType::BSM, "MASS", 37},               // mH
		{ParameterType::BSM, "HMIX", 2},                // tanb
	
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,2,1}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,2,1}},
	
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,3,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,4,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,5,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,5,2}},
	
		{ParameterType::BSM, "UMIX", 0},
		{ParameterType::BSM, "UMIX", 10},
		{ParameterType::BSM, "VMIX", 0},
		{ParameterType::BSM, "VMIX", 10}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();

		double B10pc=0.;
		double C9pc=0.;
		double B9pc=0.;
		double D9pc=0.; 

		for(int ie=0;ie<2;ie++) {
			for(int ae=0;ae<6;ae++) {
				for(int je=0;je<2;je++) {
					C9pc+=src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()*(2.*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*src.at({ParameterType::BSM, "VMIX",je*10+0})->get_val()*src.at({ParameterType::BSM, "VMIX",ie*10+0})->get_val() -f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX",je*10+0})->get_val()*src.at({ParameterType::BSM, "UMIX",ie*10+0})->get_val());
					for(int be=0;be<6;be++) {
						if (be<3) {
						B10pc+=-src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(0.5*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1}})->get_val()*f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)) +src.at({ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
						B9pc+=src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(0.5*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1}})->get_val()*f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)) -src.at({ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
						}
					}

				}
				for(int be=0;be<6;be++) {
					for(int ce=0;ce<3;ce++) {
						C9pc+=-src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()*f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*src.at({ParameterType::WILSON, "MATRIX_BSM", {2,be, ce}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2,ae, ce}})->get_val();
						}
				}
				D9pc+=pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()*h30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.));
			}
		} 		
		
		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();

		B9pc*=kappa*(mW)*mW/2./g2/g2;	

		C9pc*=-kappa/8.;

		
		double C10pH = -mass_b_muW*(src.at({ParameterType::SM, "MASS", 3})->get_val())*(tanb*tanb/8./mW/mW
		+pow(src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*tanb*tanb/4./mW/mH,2.))*f20(yt)/sw2;
		
		double C9pH =(4.*sw2-1.)*C10pH - src.at({ParameterType::SM, "MASS", 3})->get_val()*mass_b_muW/mass_top_muW/mass_top_muW*D9H0(yt,ld);
		
		
		D9pc*=kappa;

		double C9pcharg=(1.-4.*sw2)/sw2*C9pc-B9pc/sw2-D9pc;

        dep_param->set_expected(C9pH+C9pcharg);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4233, 0, 1)}, sources, func);

	// // sus_param->reset_PrimeCQG(this->get_Q_match());
    // double B10pc=0.;
	// double C9pc=0.;
	// double B9pc=0.;
	// double D9pc=0.; 

	// for(int ie=0;ie<2;ie++) {
	// 	for(int ae=0;ae<6;ae++) {
	// 		for(int je=0;je<2;je++) {
	// 			C9pc+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UR[ie][ae][2]*(2.*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f30(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*(*susy)("VMIX", je*10+0)*(*susy)("VMIX", ie*10+0) -f40(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*(*susy)("UMIX", je*10+0)*(*susy)("UMIX", ie*10+0));
	// 			for(int be=0;be<6;be++) {
	// 				if (be<3) {
	// 				B10pc+=-(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*(0.5*(*sus_param).X_NR[ie][be][1]*(*sus_param).X_NR[je][be][1]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)) +(*sus_param).X_NL[ie][be][1]*(*sus_param).X_NL[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
	// 				B9pc+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*(0.5*(*sus_param).X_NR[ie][be][1]*(*sus_param).X_NR[je][be][1]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)) -(*sus_param).X_NL[ie][be][1]*(*sus_param).X_NL[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
	// 				}
	// 			}

	// 		}
	// 		for(int be=0;be<6;be++) {
	// 			for(int ce=0;ce<3;ce++) {
	// 				C9pc+=-(*sus_param).X_UR[ie][be][1]*(*sus_param).X_UR[ie][ae][2]*f40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.))*(*sus_param).Gamma_UR[be][ce]*(*sus_param).Gamma_UR[ae][ce];
	// 				}
	// 		}
	// 		D9pc+=pow(sm("MASS",24)/(*sus_param).Mch[ie],2.)*(*sus_param).X_UR[ie][ae][1]*(*sus_param).X_UR[ie][ae][2]*h30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.));
	// 	}
	// } 		

	// B9pc*=(*sus_param).kappa*(sm("MASS",24))*sm("MASS",24)/2./sm("GAUGE", 2)/sm("GAUGE", 2);	

	// C9pc*=-(*sus_param).kappa/8.;

	
	// double C10pH = -wilson_p("WPARAM_MATCH_SM", {5,1})*(sm("MASS",3))*((*susy)("HMIX",2)*(*susy)("HMIX",2)/8./sm("MASS",24)/sm("MASS",24)
	// +pow(wilson_p("WPARAM_SI_SM", 3)*(*susy)("HMIX",2)*(*susy)("HMIX",2)/4./sm("MASS",24)/(*susy)("MASS",37),2.))*f20((*sus_param).yt)/(*sus_param).sw2;
	
	// double C9pH =(4.*(*sus_param).sw2-1.)*C10pH - sm("MASS",3)*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW*D9H0((*sus_param).yt,(*sus_param).ld);
	
	
	// D9pc*=(*sus_param).kappa;

	// double C9pcharg=(1.-4.*(*sus_param).sw2)/(*sus_param).sw2*C9pc-B9pc/(*sus_param).sw2-D9pc;
	// std::cout << "C9pcharg" << C9pcharg << std::endl;
	// std::cout << "C9pH" << C9pH << std::endl;
	// std::cout << "C9pcharg" << C9pcharg << std::endl;
	// std::cout << "total" << C9pH+ C9pcharg<< std::endl;
    // // return this->double_to_complex_save("LO", C9pH+C9pcharg);	
}

void CP10_susy::LO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
		{ParameterType::WILSON, "WPARAM_SI_SM", 3},
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::SM, "MASS", 3},
		{ParameterType::SM, "GAUGE", 2},
		{ParameterType::BSM, "MASS", 37},
		{ParameterType::BSM, "HMIX", 2},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,3,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,4,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,5,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,3,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,4,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,5,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,2,1}},
		{ParameterType::BSM, "UMIX", 0},
		{ParameterType::BSM, "UMIX", 10},
		{ParameterType::BSM, "VMIX", 0},
		{ParameterType::BSM, "VMIX", 10}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {

		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();

		double B10pc=0.;
		double C9pc=0.;

		for(int ie=0;ie<2;ie++) {
			for(int ae=0;ae<6;ae++) {
				for(int je=0;je<2;je++) {
					C9pc+=src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 1}})->get_val()*(2.*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()
							/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())
							*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()
							/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))
							*src.at({ParameterType::BSM, "VMIX",je*10+0})->get_val()*src.at({ParameterType::BSM, "VMIX",ie*10+0})->get_val() -f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()
							/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))
							*src.at({ParameterType::BSM, "UMIX",je*10+0})->get_val()*src.at({ParameterType::BSM, "UMIX",ie*10+0})->get_val());
					for(int be=0;be<6;be++) {
						if (be<3) {
							B10pc+=-src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()
									/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()
									*(0.5*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1}})->get_val()
									*f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()
									/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)) 
									+src.at({ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1}})->get_val()
									*src.at({ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())
									*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()
									/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
						}
					}
				}
				for(int be=0;be<6;be++) {
					for(int ce=0;ce<3;ce++) {
						C9pc+=-src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()
								*f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val()
								/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*src.at({ParameterType::WILSON, "MATRIX_BSM", {2,be, ce}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2,ae, ce}})->get_val();
						}
				}
			}
		} 		
		
		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		B10pc*=kappa* (mW)*mW/2./g2/g2;

		C9pc*=-kappa/8.;
		
		double C10pH = -mass_b_muW*(src.at({ParameterType::SM, "MASS", 3})->get_val())*(tanb*tanb/8./mW/mW
		+pow(src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*tanb*tanb/4./mW/mH,2.))*f20(yt)/sw2;

		double C10pcharg=(B10pc-C9pc)/sw2;

        dep_param->set_expected(C10pH+C10pcharg);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 4234, 0, 1)}, sources, func);

	// // sus_param->reset_PrimeCQG(this->get_Q_match());
    // double B10pc=0.;
	// double C9pc=0.;

	// for(int ie=0;ie<2;ie++) {
	// 	for(int ae=0;ae<6;ae++) {
	// 		for(int je=0;je<2;je++) {
	// 			C9pc+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UR[ie][ae][2]*(2.*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f30(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*(*susy)("VMIX", je*10+0)*(*susy)("VMIX", ie*10+0) -f40(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*(*susy)("UMIX", je*10+0)*(*susy)("UMIX", ie*10+0));
	// 			for(int be=0;be<6;be++) {
	// 				if (be<3) {
	// 				    B10pc+=-(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*(0.5*(*sus_param).X_NR[ie][be][1]*(*sus_param).X_NR[je][be][1]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)) +(*sus_param).X_NL[ie][be][1]*(*sus_param).X_NL[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
	// 				}
	// 			}
	// 		}
	// 		for(int be=0;be<6;be++) {
	// 			for(int ce=0;ce<3;ce++) {
	// 				C9pc+=-(*sus_param).X_UR[ie][be][1]*(*sus_param).X_UR[ie][ae][2]*f40(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.))*(*sus_param).Gamma_UR[be][ce]*(*sus_param).Gamma_UR[ae][ce];
	// 				}
	// 		}
	// 	}
	// } 		
	
	// B10pc*=(*sus_param).kappa* (sm("MASS",24))*sm("MASS",24)/2./sm("GAUGE",2)/sm("GAUGE",2);

	// C9pc*=-(*sus_param).kappa/8.;
	
	// double C10pH = -wilson_p("WPARAM_MATCH_SM", {5,1})*(sm("MASS",3))*((*susy)("HMIX",2)*(*susy)("HMIX",2)/8./sm("MASS",24)/sm("MASS",24)
	// +pow(wilson_p("WPARAM_SI_SM", 3)*(*susy)("HMIX",2)*(*susy)("HMIX",2)/4./sm("MASS",24)/(*susy)("MASS",37),2.))*f20((*sus_param).yt)/(*sus_param).sw2;

	
	// double C10pcharg=(B10pc-C9pc)/(*sus_param).sw2;
	// // return this->double_to_complex_save("LO", C10pH+C10pcharg);

}

void CPQ1_susy::LO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 1},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 11},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,0}},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,2}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
		{ParameterType::WILSON, "WPARAM_SI_SM", 3},
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::SM, "MASS", 3},
		{ParameterType::SM, "GAUGE", 2},
		{ParameterType::BSM, "MASS", 37},
		{ParameterType::BSM, "HMIX", 1},
		{ParameterType::BSM, "HMIX", 2},
		{ParameterType::WILSON, "EPSILON_SUSY", 5},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,1,2}},
		{ParameterType::BSM, "UMIX", 1},
		{ParameterType::BSM, "UMIX", 11},
		{ParameterType::BSM, "VMIX", 0},
		{ParameterType::BSM, "VMIX", 10}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double z = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 1})->get_val();
		double aY = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 11})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();

		double muQ = src.at({ParameterType::BSM, "HMIX", 1})->get_val();
		double BQ1pc1=0.;
		double BQ1pc2=0.;

		double Dp{0},Dm{0};
		double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a1{0}, NQ1pc{0}, NQ2pc{0};

		for(int ie=0;ie<2;ie++) {
			for(int ae=0;ae<6;ae++) {

				for(int je=0;je<2;je++) {
					for(int be=0;be<6;be++) {
						if (be<3) {
							BQ1pc1+=src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()
									/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1}})->get_val()
									*f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()
									/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()
									/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))); 	
							BQ1pc2+=src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()
									/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1}})->get_val()
									*src.at({ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())
									*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()
									/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
						}
						for(int me=0;me<3;me++) {
							for(int ne=0;ne<3;ne++) {
								Dp=0.;
								Dm=0.;
								for(int fe=0;fe<3;fe++) { 
									Dp+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2,ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1,be, fe}})->get_val()+src.at({ParameterType::WILSON, "MATRIX_BSM", {1,ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2,be, fe}})->get_val());
									Dm+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2,ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1,be, fe}})->get_val()-src.at({ParameterType::WILSON, "MATRIX_BSM", {1,ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2,be, fe}})->get_val());
								}
								a0a=-(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val())*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX",ie*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX",je*10+0})->get_val())*kron(ae,be);
								
								a0b=-(f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX",je*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX",ie*10+0})->get_val())*kron(ae,be);
								a0c=1./mW*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je);
								a0Q1=a0a+a0b+Dp*a0c;
								a0Q2=-a0a+a0b+Dm*a0c;
								a1=src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/sqrt(2.)/mW*f80(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je)*kron(ae,be);
								
								NQ1pc+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1,be, me}})->get_val()*src.at({ParameterType::BSM, "UMIX",je*10+1})->get_val()*(a0Q1+a1*tanb);
								NQ2pc+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1,be, me}})->get_val()*src.at({ParameterType::BSM, "UMIX",je*10+1})->get_val()*(a0Q2+a1*tanb);
							}
						}
					}

				}
			}
		}		
		
		/* Wilson coefficients CQ1 and CQ2 prime */ 
		double NQ1pH=-src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/4./mW/mW*xt*f30(xt,z);
		
		double BQ1pH=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/4./mW/mW*f70(xt,z);
		
		complex_t CQ1pH=(NQ1pH+BQ1pH)*src.at({ParameterType::SM, "MASS", 3})->get_val()/sw2;
		
		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		double BQ1pc=(BQ1pc1+BQ1pc2)*kappa*(mW)*mW/2./g2/g2/sw2;

		NQ1pc*=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*(src.at({ParameterType::SM, "MASS", 3})->get_val())/sw2;


		complex_t CQ1pcharg=NQ1pc+BQ1pc;

		// this->set_WilsonCoeffMatching("LO", (CQ1pcharg+CQ1pH)/sus_param->epsfac);
		// return (CQ1pcharg+CQ1pH)/sus_param->epsfac;
		double epsfac = src.at({ParameterType::WILSON, "EPSILON_SUSY", 5})->get_val();

        dep_param->set_expected((CQ1pcharg+CQ1pH)/epsfac);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3130, 0, 1)}, sources, func);


	// // sus_param->reset_PrimeCQG(this->get_Q_match());
	// double BQ1pc1=0.;
	// double BQ1pc2=0.;

    // double Dp{0},Dm{0};
	// double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a1{0}, NQ1pc{0}, NQ2pc{0};

	// for(int ie=0;ie<2;ie++) {
	// 	for(int ae=0;ae<6;ae++) {

	// 		for(int je=0;je<2;je++) {
	// 			for(int be=0;be<6;be++) {
	// 				if (be<3) {
    //                     BQ1pc1+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UL[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][1]*(*sus_param).X_NR[je][be][1]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))); 	
    //                     BQ1pc2+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UL[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][1]*(*sus_param).X_NL[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
	// 				}
	// 				for(int me=0;me<3;me++) {
	// 					for(int ne=0;ne<3;ne++) {
	// 						Dp=0.;
	// 						Dm=0.;
	// 						for(int fe=0;fe<3;fe++) { 
	// 							Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[be][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[be][fe]);
	// 							Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[be][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[be][fe]);
	// 						}
	// 						a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,be);

	// 						a0b=-(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", je*10+1)*(*susy)("VMIX", ie*10+0))*kron(ae,be);
	// 						a0c=1./sm("MASS",24)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.))*kron(ie,je);
	// 						a0Q1=a0a+a0b+Dp*a0c;
	// 						a0Q2=-a0a+a0b+Dm*a0c;
	// 						a1=(*sus_param).Mch[ie]/sqrt(2.)/sm("MASS",24)*f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,be);
							
	// 						NQ1pc+=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[be][me]*(*susy)("UMIX", je*10+1)*(a0Q1+a1*(*susy)("HMIX",2));
	// 						NQ2pc+=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[be][me]*(*susy)("UMIX", je*10+1)*(a0Q2+a1*(*susy)("HMIX",2));
	// 					}
	// 				}
	// 			}

	// 		}
	// 	}
	// }		

	// /* Wilson coefficients CQ1 and CQ2 prime */ 
	// double NQ1pH=-wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./sm("MASS",24)/sm("MASS",24)*(*sus_param).xt*f30((*sus_param).xt,(*sus_param).z);
	
	// double BQ1pH=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./sm("MASS",24)/sm("MASS",24)*f70((*sus_param).xt,(*sus_param).z);
	
	// complex_t CQ1pH=(NQ1pH+BQ1pH)*sm("MASS",3)/(*sus_param).sw2;
	

	// double BQ1pc=(BQ1pc1+BQ1pc2)*(*sus_param).kappa*(sm("MASS",24))*sm("MASS",24)/2./sm("GAUGE", 2)/sm("GAUGE", 2)/(*sus_param).sw2;

	// NQ1pc*=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*(sm("MASS",3))/(*sus_param).sw2;


	// complex_t CQ1pcharg=NQ1pc+BQ1pc;

    // this->set_WilsonCoeffMatching("LO", (CQ1pcharg+CQ1pH)/sus_param->epsfac);
    // // return (CQ1pcharg+CQ1pH)/sus_param->epsfac;
}

void CPQ2_susy::LO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 1},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 11},
		{ParameterType::WILSON, "WPARAM_SI_SM", 3},
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,0}},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,2}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
		{ParameterType::SM, "MASS", 3},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::SM, "GAUGE", 2},
		{ParameterType::BSM, "MASS", 37},
		{ParameterType::BSM, "HMIX", 1},
		{ParameterType::BSM, "HMIX", 2},
		{ParameterType::WILSON, "EPSILON_SUSY", 5},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {2,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,0,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,1,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {3,1,2,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {4,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {5,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,0,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {6,1,2,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,0,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,0,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,0,2}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,1,0}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,1,1}},
		{ParameterType::WILSON, "MATRIX_BSM", {12,0,0,1,2}},
		{ParameterType::BSM, "UMIX", 1},
		{ParameterType::BSM, "UMIX", 11},
		{ParameterType::BSM, "VMIX", 0},
		{ParameterType::BSM, "VMIX", 10}
	};

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double z = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 1})->get_val();
		double aY = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 11})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();

		double muQ = src.at({ParameterType::BSM, "HMIX", 1})->get_val();

		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		double epsfac = src.at({ParameterType::WILSON, "EPSILON_SUSY", 5})->get_val();

		double BQ1pc1=0.;
		double BQ1pc2=0.;

		double Dp{0},Dm{0};
		double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a1{0}, NQ1pc{0}, NQ2pc{0};

		for(int ie=0;ie<2;ie++) {
			for(int ae=0;ae<6;ae++) {

				for(int je=0;je<2;je++) {

					for(int be=0;be<6;be++) {
						if (be<3) {

						BQ1pc1+=src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {5,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,je, be, 1}})->get_val()*f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))); 	
						BQ1pc2+=src.at({ParameterType::WILSON, "MATRIX_BSM", {4,je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {3,ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5,je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
						}
						for(int me=0;me<3;me++) {
							for(int ne=0;ne<3;ne++) {
								Dp=0.;
								Dm=0.;
								for(int fe=0;fe<3;fe++) { 
									Dp+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, be, fe}})->get_val()+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, be, fe}})->get_val());
									Dm+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, be, fe}})->get_val()-src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, be, fe}})->get_val());
								}
								a0a=-(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val())*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX",ie*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX",je*10+0})->get_val())*kron(ae,be);

								a0b=-(f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX",je*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX",ie*10+0})->get_val())*kron(ae,be);
								a0c=1./mW*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je);
								a0Q1=a0a+a0b+Dp*a0c;
								a0Q2=-a0a+a0b+Dm*a0c;
								a1=src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/sqrt(2.)/mW*f80(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je)*kron(ae,be);
								
								NQ1pc+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, be, me}})->get_val()*src.at({ParameterType::BSM, "UMIX",je*10+1})->get_val()*(a0Q1+a1*tanb);
								NQ2pc+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, be, me}})->get_val()*src.at({ParameterType::BSM, "UMIX",je*10+1})->get_val()*(a0Q2+a1*tanb);
							}
						}
					}

				}

			}
		} 		

		/* Wilson coefficients CQ1 and CQ2 prime */ 
		double NQ1pH=-src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/4./mW/mW*xt*f30(xt,z);
		
		double BQ1pH=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/4./mW/mW*f70(xt,z);
		
		complex_t CQ1pH=(NQ1pH+BQ1pH)*src.at({ParameterType::SM, "MASS", 3})->get_val()/sw2;
		
		complex_t CQ2pH=CQ1pH;

		
		double BQ2pc=(BQ1pc1-BQ1pc2)*kappa*(mW)*mW/2./g2/g2/sw2;

		NQ2pc*=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*(src.at({ParameterType::SM, "MASS", 3})->get_val())/sw2;


		complex_t CQ2pcharg=NQ2pc+BQ2pc;

        dep_param->set_expected((CQ2pH+CQ2pcharg)/epsfac);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3133, 0, 1)}, sources, func);

	// // sus_param->reset_PrimeCQG(this->get_Q_match());
	// double BQ1pc1=0.;
	// double BQ1pc2=0.;

    // double Dp{0},Dm{0};
	// double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a1{0}, NQ1pc{0}, NQ2pc{0};

	// for(int ie=0;ie<2;ie++) {
	// 	for(int ae=0;ae<6;ae++) {

	// 		for(int je=0;je<2;je++) {

	// 			for(int be=0;be<6;be++) {
	// 				if (be<3) {

	// 				BQ1pc1+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UL[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][1]*(*sus_param).X_NR[je][be][1]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))); 	
	// 				BQ1pc2+=(*sus_param).X_UR[je][ae][1]*(*sus_param).X_UL[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][1]*(*sus_param).X_NL[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
	// 				}
	// 				for(int me=0;me<3;me++) {
	// 					for(int ne=0;ne<3;ne++) {
	// 						Dp=0.;
	// 						Dm=0.;
	// 						for(int fe=0;fe<3;fe++) { 
	// 							Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[be][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[be][fe]);
	// 							Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[be][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[be][fe]);
	// 						}
	// 						a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,be);

	// 						a0b=-(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", je*10+1)*(*susy)("VMIX", ie*10+0))*kron(ae,be);
	// 						a0c=1./sm("MASS",24)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.))*kron(ie,je);
	// 						a0Q1=a0a+a0b+Dp*a0c;
	// 						a0Q2=-a0a+a0b+Dm*a0c;
	// 						a1=(*sus_param).Mch[ie]/sqrt(2.)/sm("MASS",24)*f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,be);
							
	// 						NQ1pc+=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[be][me]*(*susy)("UMIX", je*10+1)*(a0Q1+a1*(*susy)("HMIX",2));
	// 						NQ2pc+=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[be][me]*(*susy)("UMIX", je*10+1)*(a0Q2+a1*(*susy)("HMIX",2));
	// 					}
	// 				}
	// 			}

	// 		}

	// 	}
	// } 		

	// /* Wilson coefficients CQ1 and CQ2 prime */ 
	// double NQ1pH=-wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./sm("MASS",24)/sm("MASS",24)*(*sus_param).xt*f30((*sus_param).xt,(*sus_param).z);
	
	// double BQ1pH=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./sm("MASS",24)/sm("MASS",24)*f70((*sus_param).xt,(*sus_param).z);
	
	// complex_t CQ1pH=(NQ1pH+BQ1pH)*sm("MASS",3)/(*sus_param).sw2;
	
	// complex_t CQ2pH=CQ1pH;

	
	// double BQ2pc=(BQ1pc1-BQ1pc2)*(*sus_param).kappa*(sm("MASS",24))*sm("MASS",24)/2./sm("GAUGE", 2)/sm("GAUGE", 2)/(*sus_param).sw2;

	// NQ2pc*=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*(sm("MASS",3))/(*sus_param).sw2;


	// complex_t CQ2pcharg=NQ2pc+BQ2pc;
    // this->set_WilsonCoeffMatching("LO", (CQ2pH+CQ2pcharg)/sus_param->epsfac);
    // // return (CQ2pH+CQ2pcharg)/sus_param->epsfac;
}

void CQ1_susy::LO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 1},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 2},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 3},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 9},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 11},
		{ParameterType::WILSON, "WPARAM_SI_SM", 3},
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,0}},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,2}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", 6},
		{ParameterType::WILSON, "EW_SCALE", 1},
		{ParameterType::SM, "MASS", 3},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::SM, "GAUGE", 2},
		{ParameterType::BSM, "MASS", 37},
		{ParameterType::BSM, "HMIX", 1},
		{ParameterType::BSM, "HMIX", 2},
		{ParameterType::WILSON, "EPSILON_SUSY", 5},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {13,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,2}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,3}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,4}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {14,5}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,0}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,1}},
		{ParameterType::WILSON, "WPARAM_SI_BSM", {16,2}},
	};
	
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 3; ++j) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, i, j}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {2, i, j}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {9, i, j}});
		}
	}
	
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 6; ++j) {
			for (int k = 0; k < 2; ++k) {
				sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, i, j, k}});
				sources.insert({ParameterType::WILSON, "MATRIX_BSM", {4, i, j, k}});
			}
		}
	}
	
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j) {
			for (int k = 0; k < 3; ++k) {
				sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 0, j, k}});
				sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 1, j, k}});
			}
		}
	}
	
	for (int idx = 0; idx < 20; ++idx) {
		sources.insert({ParameterType::BSM, "UMIX", idx});
		sources.insert({ParameterType::BSM, "VMIX", idx});
	}

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double z = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 1})->get_val();
		double aY = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 11})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();

		double muQ = src.at({ParameterType::BSM, "HMIX", 1})->get_val();
		double xH = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 2})->get_val();
        double xH0 = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 3})->get_val();
		double alpha = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 9})->get_val();
        double beta = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();

		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		double epsfac = src.at({ParameterType::WILSON, "EPSILON_SUSY", 5})->get_val();

		complex_t BQ10c1=0.;
		complex_t BQ10c2=0.;

		double Dp, Dm;
		double a0a{0}, a0b{0}, a0c, a0Q1{0}, a0Q2{0};
		double a1{0};

		complex_t NQ10c=0.;

		for(int ie=0;ie<2;ie++) {
			for(int je=0;je<2;je++) {
				for(int ae=0;ae<6;ae++) {
					for(int be=0;be<3;be++) { 
						BQ10c1+=src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}})->get_val()*f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
						BQ10c2+=src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
						for(int me=0;me<6;me++) {
							for(int ne=0;ne<3;ne++) {
								Dp=0.;
								Dm=0.;
								for(int fe=0;fe<3;fe++) 
								{
									Dp+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, fe}})->get_val()+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}})->get_val());
									Dm+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, fe}})->get_val()-src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}})->get_val());
								}
								a0a=-(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val())*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX", ie*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX", je*10+0})->get_val())*kron(ae,me);
								a0b=-(f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX", ie*10+0})->get_val())*kron(ae,me);
								a0c=1./mW*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je);
								a0Q1=a0a+a0b+Dp*a0c;
								
								a1=src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/sqrt(2.)/mW*f80(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je)*kron(ae,me);
								NQ10c+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*(a0Q1+a1*tanb);
							}
						
						}
					}
				}
			}
		}
		complex_t BQ10c=(BQ10c1+BQ10c2)*kappa*mW*mW/2./g2/g2/sw2;

		NQ10c*=src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
		double le = -tanb;
		double G1=-3./4.+ld*lu*F4SP(xt,xH)+lu*lu*F5SP(xt,xH);
		double G2=ld*(ld*lu+1.)*F6SP(xt,xH)-ld*lu*lu*F7SP(xt,xH)
		+lu*lu*(ld*F8SP(xt,xH)+lu*F9SP(xt,xH)-lu*F10SP(xt,xH))+lu*F11SP(xt,xH)-lu*F12SP(xt,xH);

		double CSn_2HDM=xt*(F0SP(xt)+le*(ld*F1SP(xt,xH)+lu*F2SP(xt,xH))+le*lu*F3SP(xt,xH))
		+xt/2./xH*(sin(alpha-beta)+cos(alpha-beta)*le)*(sin(alpha-beta)*G1+cos(alpha-beta)*G2)
		+xt/2./xH0*(cos(alpha-beta)-sin(alpha-beta)*le)*(cos(alpha-beta)*G1-sin(alpha-beta)*G2);
		complex_t CQ1H_0=CSc_2HDM(xH,xt,lu,ld,le)+CSn_2HDM;
		CQ1H_0*=(src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*mass_b_muW/mW/mW)/sw2;

		complex_t CQ1charg_0=NQ10c+BQ10c;
		complex_t coeff_temp = (CQ1charg_0+CQ1H_0)/epsfac;

		/* NMSSM */

		double lambdaNMSSM = 1;
		double lambdaSNMSSM = 1;
		double AlambdaNSSM = 1;
		double kappaNMSSM = 1;
		double m_Bs = 1;
		double mass_nutl = 1;

		if(src.at({ParameterType::BSM, "MASS",  46})->get_val()!=0.||src.at({ParameterType::BSM, "MASS",  45})->get_val()!=0.)
		{
			LOG_INFO("NMSSM ? Doesn't exist, don't search for it.");

			double s=lambdaSNMSSM/lambdaNMSSM;
			double v=sqrt(1./sqrt(2.)/src.at({ParameterType::SM, "SMINPUTS",  2})->get_val());
			
			double v_deltam_s=v/s*(sqrt(2.)*AlambdaNSSM-2.*kappaNMSSM*s)/(sqrt(2.)*AlambdaNSSM+kappaNMSSM*s);
			
			double mH0[4],mA0[3],mstop[3];
			
			mstop[0]=src.at({ParameterType::BSM, "MASS",  2000013})->get_val(); //mass upr, is that right ?
			mstop[1]=src.at({ParameterType::BSM, "MASS",  1000006})->get_val();
			mstop[2]=src.at({ParameterType::BSM, "MASS",  2000006})->get_val();
			
			mH0[1]=src.at({ParameterType::BSM, "MASS",  25})->get_val();
			mH0[2]=src.at({ParameterType::BSM, "MASS",  35})->get_val();
			mH0[3]=src.at({ParameterType::BSM, "MASS",  36})->get_val();
			mA0[1]=src.at({ParameterType::BSM, "MASS",  36})->get_val();
			mA0[2]=src.at({ParameterType::BSM, "MASS",  36})->get_val(); //TODO WTF
			
			double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
			double T1[3][4][4],T2[4][4][4];
			std::array<std::array<double,4>,4> TU;
		
			TU[1][1]=1.;
			for(int ie=0;ie<2;ie++){
				for(int je=0;je<2;je++) {
					TU[ie+1][je+1]=src.at({ParameterType::BSM, "STOPMIX",  ie*10+je})->get_val();
				}
			}

			
			double vu=sqrt(pow(sin(atan(tanb)),2.)/sqrt(2.)/src.at({ParameterType::SM, "SMINPUTS",  2})->get_val());
			double vd=vu/tanb;

			for(int je=0;je<2;je++) {
				for(int le=0;le<2;le++) {
					for(int ae=0;ae<3;ae++) {
						if (ae <3 ){
							Ralj[ae][le][je]=-g2/sqrt(2.)*(src.at({ParameterType::BSM, "A0MIX", ae*10+1})->get_val()*src.at({ParameterType::BSM, "UMIX", 20+le})->get_val()*src.at({ParameterType::BSM, "VMIX", 20+je})->get_val()+src.at({ParameterType::BSM, "A0MIX", ae*10+2})->get_val()*src.at({ParameterType::BSM, "UMIX", 10+le})->get_val()*src.at({ParameterType::BSM, "VMIX", 20+je})->get_val())-lambdaNMSSM/sqrt(2.)*src.at({ParameterType::BSM, "A0MIX", ae*10+3})->get_val()*src.at({ParameterType::BSM, "UMIX", 20+le})->get_val()*src.at({ParameterType::BSM, "VMIX", 20+je})->get_val();
						}
						Qalj[ae][le][je]=g2/sqrt(2.)*(src.at({ParameterType::BSM, "H0MIX", ae*10+1})->get_val()*src.at({ParameterType::BSM, "UMIX", 20+le})->get_val()*src.at({ParameterType::BSM, "VMIX", 20+je})->get_val()+src.at({ParameterType::BSM, "H0MIX", ae*10+2})->get_val()*src.at({ParameterType::BSM, "UMIX", 10+le})->get_val()*src.at({ParameterType::BSM, "VMIX", 20+je})->get_val())-lambdaNMSSM/sqrt(2.)*src.at({ParameterType::BSM, "H0MIX", ae*10+3})->get_val()*src.at({ParameterType::BSM, "UMIX", 20+le})->get_val()*src.at({ParameterType::BSM, "VMIX", 20+je})->get_val();
						for(int ke=1;ke<=3;ke++) {
							G1[ae][ke][je][le]=(TU[ae][2]*TU[ke][2]-kron(ae,1)*kron(ke,1))*src.at({ParameterType::BSM, "VMIX", 10+le})->get_val()*src.at({ParameterType::BSM, "UMIX", 20+je})->get_val()-mass_top_muW/sqrt(2.)/sin(atan(tanb))/mW*TU[ae][3]*TU[ke][2]*src.at({ParameterType::BSM, "VMIX", 20+le})->get_val()*src.at({ParameterType::BSM, "UMIX", 20+je})->get_val();
						}
					}
				}
				for(int ie=0;ie<3;ie++) {
					for(int ke=0;ke<3;ke++) {
						T1[je][ie][ke]=(TU[ie][3]*TU[ke][2]-TU[ie][2]*TU[ke][3])*((lambdaNMSSM/sqrt(2.)*(vd*src.at({ParameterType::BSM, "A0MIX", je*10+3})->get_val()+s*src.at({ParameterType::BSM, "A0MIX", je*10+1})->get_val()))-src.at({ParameterType::BSM, "AU", 11})->get_val()*src.at({ParameterType::BSM, "A0MIX", je*10+2})->get_val());
					}
				}
			}

			complex_t CQ1H=0.;
			complex_t CQ2H=0.;
			complex_t CQ1c=0.;
			complex_t CQ2c=0.;
			complex_t CAc=0.;

			for(int ae=0;ae<3;ae++) {
				for(int ie=0;ie<3;ie++) {
					for(int ke=0;ke<3;ke++){
						T2[ae][ie][ke]=-mass_top_muW/2./mW*(2.*mass_top_muW*src.at({ParameterType::BSM, "H0MIX", ae*10+2})->get_val()*(TU[ie][2]*TU[ke][2]+TU[ie][3]*TU[ke][3])	+((lambdaNMSSM/sqrt(2.)*(vd*src.at({ParameterType::BSM, "H0MIX", ae*10+3})->get_val()+s*src.at({ParameterType::BSM, "H0MIX", ae*10+1})->get_val()))+src.at({ParameterType::BSM, "AU", 11})->get_val()*src.at({ParameterType::BSM, "H0MIX", ae*10+2})->get_val())*(TU[ie][3]*TU[ke][2]+TU[ie][2]*TU[ke][3]))	+src.at({ParameterType::SM, "MASS",  23})->get_val()/2./sqrt(1.-sw2)*(1.-4./3.*sw2)*src.at({ParameterType::BSM, "H0MIX", ae*10+2})->get_val()*(TU[ie][1]*TU[ke][1]+TU[ie][2]*TU[ke][2])+2./3.*mW*sw2/(1.-sw2)*src.at({ParameterType::BSM, "H0MIX", ae*10+2})->get_val()*TU[ie][3]*TU[ke][3];

						for(int je=0;je<2;je++) {
							for(int le=0;le<2;le++) {
								CQ1c+=G1[ie][ke][je][le]/mH0[ae]/mH0[ae]*( sqrt(2.)*src.at({ParameterType::BSM, "H0MIX", ae*10+1})->get_val()*src.at({ParameterType::BSM, "H0MIX", ae*10+1})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mW/cos(atan(tanb))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))
								-2.*sqrt(2.)*src.at({ParameterType::BSM, "H0MIX", ae*10+1})->get_val()/g2*kron(ie,ke)*(Qalj[ae][le][je]*f40(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))+src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()*Qalj[ae][je][le]*f30(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.)))		+2.*sqrt(2.)*src.at({ParameterType::BSM, "H0MIX", ae*10+1})->get_val()*T2[ae][ie][ke]*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mstop[ke-1],2.))
								+mH0[ae]*mH0[ae]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*kron(ie,ke)*(src.at({ParameterType::BSM, "UMIX", 20+je})->get_val()*src.at({ParameterType::BSM, "VMIX", 10+le})->get_val()*f50(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(mass_nutl/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))
								-src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*src.at({ParameterType::BSM, "UMIX", 20+le})->get_val()*src.at({ParameterType::BSM, "VMIX", 10+je})->get_val()* f60(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(mass_nutl/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))));

								CQ2c+=G1[ie][ke][je][le]/mA0[ae]/mA0[ae]*(sqrt(2.)*src.at({ParameterType::BSM, "A0MIX", ae*10+1})->get_val()*src.at({ParameterType::BSM, "A0MIX", ae*10+1})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mW/cos(atan(tanb))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))
								-2.*sqrt(2.)*src.at({ParameterType::BSM, "A0MIX", ae*10+1})->get_val()/g2*kron(ie,ke)*(-Ralj[ae][le][je]*f40(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))+src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()*Ralj[ae][je][le]*f30(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.)))			-sqrt(2.)*src.at({ParameterType::BSM, "A0MIX", ae*10+1})->get_val()*T1[ae][ie][ke]*mass_top_muW*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mstop[ke-1],2.))
								+mA0[ae]*mA0[ae]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*kron(ie,ke)*(src.at({ParameterType::BSM, "UMIX", 20+je})->get_val()*src.at({ParameterType::BSM, "VMIX", 10+le})->get_val()*f50(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(mass_nutl/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))
								-src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*src.at({ParameterType::BSM, "UMIX", 20+le})->get_val()*src.at({ParameterType::BSM, "VMIX", 10+je})->get_val()*f60(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(mass_nutl/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))));
							}		
						}

					}
				}
				CQ1H+=(mH*mH/mW/mW*src.at({ParameterType::BSM, "H0MIX", ae*10+1})->get_val()*src.at({ParameterType::BSM, "H0MIX", ae*10+1})->get_val()*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)	+mass_top_muW*mass_top_muW*mH0[ae]*mH0[ae]/mW/mW/mH/mH*f30(mass_top_muW*mass_top_muW/mH/mH,mass_top_muW*mass_top_muW/mW/mW))/mH0[ae]/mH0[ae];
				if (ae < 3) {
					CQ2H+=((mH*mH/mW/mW*src.at({ParameterType::BSM, "A0MIX", ae*10+1})->get_val()*src.at({ParameterType::BSM, "A0MIX", ae*10+1})->get_val()+kron(ae,2)*src.at({ParameterType::BSM, "A0MIX", ae*10+1})->get_val())*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)	+mass_top_muW*mass_top_muW*mA0[ae]*mA0[ae]/mW/mW/mH/mH*f30(mass_top_muW*mass_top_muW/mH/mH,mass_top_muW*mass_top_muW/mW/mW))/mA0[ae]/mA0[ae];
				}
				for(int je=0;je<2;je++) {
					for(int le=0;le<2;le++) {
						CAc = complex_t(CAc.real(), CAc.imag()+(tanb)/sqrt(2.)*G1[ae][ae][je][le]*(v_deltam_s*kron(le,je)*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mW)*f80(pow(mstop[ae-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))-(Ralj[1][je][le]*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val())*f30(pow(mstop[ae-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))-Ralj[1][le][je]*f40(pow(mstop[ae-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.)))));
						}
					}
			}
		
			CQ1H*=-src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()/4.*tanb*tanb;
			CQ2H*=src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()/4.*tanb*tanb;
			
			complex_t CAH={0,-lambdaNMSSM*AlambdaNSSM/g2/mW*tanb*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)};
		
				
			CQ1c*=src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()/4.*tanb*tanb;
			CQ2c*=-src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()/4.*tanb*tanb;		
		
			coeff_temp = (CQ1H+CQ1c)*mass_b_muW/sw2/epsfac;

		}

        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3230, 0, 1)}, sources, func);

	// // sus_param->reset_PrimeCQG(this->get_Q_match());
	// complex_t BQ10c1=0.;
	// complex_t BQ10c2=0.;

	// double Dp, Dm;
	// double a0a{0}, a0b{0}, a0c, a0Q1{0}, a0Q2{0};
	// double a1{0};

	// complex_t NQ10c=0.;

	// for(int ie=0;ie<2;ie++) {
	// 	for(int je=0;je<2;je++) {
	// 		for(int ae=0;ae<6;ae++) {
	// 			for(int be=0;be<3;be++) { 
	// 				BQ10c1+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][1]*(*sus_param).X_NL[je][be][1]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
	// 				BQ10c2+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][1]*(*sus_param).X_NR[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
	// 				for(int me=0;me<6;me++) {
	// 					for(int ne=0;ne<3;ne++) {
	// 						Dp=0.;
	// 						Dm=0.;
	// 						for(int fe=0;fe<3;fe++) 
	// 						{
	// 							Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
	// 							Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
	// 						}
	// 						a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,me);
	// 						a0b=-(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", je*10+1)*(*susy)("VMIX", ie*10+0))*kron(ae,me);
	// 						a0c=1./sm("MASS",24)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))*kron(ie,je);
	// 						a0Q1=a0a+a0b+Dp*a0c;
							
	// 						a1=(*sus_param).Mch[ie]/sqrt(2.)/sm("MASS",24)*f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me);
	// 						NQ10c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q1+a1*(*susy)("HMIX",2));
	// 					}
					
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// complex_t BQ10c=(BQ10c1+BQ10c2)*(*sus_param).kappa*sm("MASS",24)*sm("MASS",24)/2./sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;

	// NQ10c*=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;
    // double le = -(*susy)("HMIX", 2);
	// double G1=-3./4.+sus_param->ld*sus_param->lu*F4SP(sus_param->xt,sus_param->xH)+sus_param->lu*sus_param->lu*F5SP(sus_param->xt,sus_param->xH);
	// double G2=sus_param->ld*(sus_param->ld*sus_param->lu+1.)*F6SP(sus_param->xt,sus_param->xH)-sus_param->ld*sus_param->lu*sus_param->lu*F7SP(sus_param->xt,sus_param->xH)
	// +sus_param->lu*sus_param->lu*(sus_param->ld*F8SP(sus_param->xt,sus_param->xH)+sus_param->lu*F9SP(sus_param->xt,sus_param->xH)-sus_param->lu*F10SP(sus_param->xt,sus_param->xH))+sus_param->lu*F11SP(sus_param->xt,sus_param->xH)-sus_param->lu*F12SP(sus_param->xt,sus_param->xH);

	// double CSn_2HDM=sus_param->xt*(F0SP(sus_param->xt)+le*(sus_param->ld*F1SP(sus_param->xt,sus_param->xH)+sus_param->lu*F2SP(sus_param->xt,sus_param->xH))+le*sus_param->lu*F3SP(sus_param->xt,sus_param->xH))
	// +sus_param->xt/2./sus_param->xh*(sin(sus_param->alpha-sus_param->beta)+cos(sus_param->alpha-sus_param->beta)*le)*(sin(sus_param->alpha-sus_param->beta)*G1+cos(sus_param->alpha-sus_param->beta)*G2)
	// +sus_param->xt/2./sus_param->xH0*(cos(sus_param->alpha-sus_param->beta)-sin(sus_param->alpha-sus_param->beta)*le)*(cos(sus_param->alpha-sus_param->beta)*G1-sin(sus_param->alpha-sus_param->beta)*G2);
	// complex_t CQ1H_0=CSc_2HDM(sus_param->xH,sus_param->xt,sus_param->lu,sus_param->ld,le)+CSn_2HDM;
	// CQ1H_0*=(wilson_p("WPARAM_SI_SM", 3)*sus_param->mass_b_muW/sm("MASS",24)/sm("MASS",24))/sus_param->sw2;

	// complex_t CQ1charg_0=NQ10c+BQ10c;
    // complex_t coeff_temp = (CQ1charg_0+CQ1H_0)/sus_param->epsfac;
    // this->set_WilsonCoeffMatching("LO", coeff_temp);

    // /* NMSSM */

	// double lambdaNMSSM = 1;
	// double lambdaSNMSSM = 1;
	// double AlambdaNSSM = 1;
	// double kappaNMSSM = 1;
	// double m_Bs = 1;
	// double mass_nutl = 1;

    // if((*susy)("MASS",46)!=0.||(*susy)("MASS",45)!=0.)
	// {
	// 	LOG_INFO("NMSSM ? Doesn't exist, don't search for it.");

	// 	double s=lambdaSNMSSM/lambdaNMSSM;
	// 	double v=sqrt(1./sqrt(2.)/sm("SMINPUTS", 2));
		
	// 	double v_deltam_s=v/s*(sqrt(2.)*AlambdaNSSM-2.*kappaNMSSM*s)/(sqrt(2.)*AlambdaNSSM+kappaNMSSM*s);
		
	// 	double mH0[4],mA0[3],mstop[3];
		
	// 	mstop[0]=(*susy)("MASS", 2000013); //mass upr, is that right ?
	// 	mstop[1]=(*susy)("MASS", 1000006);
	// 	mstop[2]=(*susy)("MASS", 2000006);
		
	// 	mH0[1]=(*susy)("MASS", 25);
	// 	mH0[2]=(*susy)("MASS", 35);
	// 	mH0[3]=(*susy)("MASS", 36);
	// 	mA0[1]=(*susy)("MASS",36);
	// 	mA0[2]=(*susy)("MASS",36);
		
	// 	double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
	// 	double T1[3][4][4],T2[4][4][4];
	// 	std::array<std::array<double,4>,4> TU;
	
	// 	TU[1][1]=1.;
	// 	for(int ie=0;ie<2;ie++){
	// 		 for(int je=0;je<2;je++) {
	// 			TU[ie+1][je+1]=(*susy)("STOPMIX", ie*10+je);
	// 		}
	// 	}

		
	// 	double vu=sqrt(pow(sin(atan((*susy)("HMIX",2))),2.)/sqrt(2.)/sm("SMINPUTS", 2));
	// 	double vd=vu/(*susy)("HMIX",2);

	// 	for(int je=0;je<2;je++) {
	// 		for(int le=0;le<2;le++) {
	// 			 for(int ae=0;ae<3;ae++) {
	// 				if (ae <3 ){
	// 					Ralj[ae][le][je]=-sm("GAUGE",2)/sqrt(2.)*((*susy)("A0MIX",ae*10+1)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je)+(*susy)("A0MIX",ae*10+2)*(*susy)("UMIX",10+le)*(*susy)("VMIX",20+je))-lambdaNMSSM/sqrt(2.)*(*susy)("A0MIX",ae*10+3)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je);
	// 				}
	// 				Qalj[ae][le][je]=sm("GAUGE",2)/sqrt(2.)*((*susy)("H0MIX",ae*10+1)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je)+(*susy)("H0MIX",ae*10+2)*(*susy)("UMIX",10+le)*(*susy)("VMIX",20+je))-lambdaNMSSM/sqrt(2.)*(*susy)("H0MIX",ae*10+3)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je);
	// 				for(int ke=1;ke<=3;ke++) {
	// 					G1[ae][ke][je][le]=(TU[ae][2]*TU[ke][2]-kron(ae,1)*kron(ke,1))*(*susy)("VMIX",10+le)*(*susy)("UMIX",20+je)-(*sus_param).mass_top_muW/sqrt(2.)/sin(atan((*susy)("HMIX",2)))/sm("MASS",24)*TU[ae][3]*TU[ke][2]*(*susy)("VMIX",20+le)*(*susy)("UMIX",20+je);
	// 				}
	// 			}
	// 		}
	// 		for(int ie=0;ie<3;ie++) {
	// 			for(int ke=0;ke<3;ke++) {
	// 				T1[je][ie][ke]=(TU[ie][3]*TU[ke][2]-TU[ie][2]*TU[ke][3])*((lambdaNMSSM/sqrt(2.)*(vd*(*susy)("A0MIX",je*10+3)+s*(*susy)("A0MIX",je*10+1)))-(*susy)("AU", 11)*(*susy)("A0MIX",je*10+2));
	// 			}
	// 		}
	// 	}

	// 	complex_t CQ1H=0.;
	// 	complex_t CQ2H=0.;
	// 	complex_t CQ1c=0.;
	// 	complex_t CQ2c=0.;
	// 	complex_t CAc=0.;

	// 	for(int ae=0;ae<3;ae++) {
	// 		for(int ie=0;ie<3;ie++) {
	// 			for(int ke=0;ke<3;ke++){
	// 				T2[ae][ie][ke]=-(*sus_param).mass_top_muW/2./sm("MASS",24)*(2.*(*sus_param).mass_top_muW*(*susy)("H0MIX",ae*10+2)*(TU[ie][2]*TU[ke][2]+TU[ie][3]*TU[ke][3])	+((lambdaNMSSM/sqrt(2.)*(vd*(*susy)("H0MIX",ae*10+3)+s*(*susy)("H0MIX",ae*10+1)))+(*susy)("AU", 11)*(*susy)("H0MIX",ae*10+2))*(TU[ie][3]*TU[ke][2]+TU[ie][2]*TU[ke][3]))	+src.at({ParameterType::SM, "MASS",  23})->get_val()/2./sqrt(1.-(*sus_param).sw2)*(1.-4./3.*(*sus_param).sw2)*(*susy)("H0MIX",ae*10+2)*(TU[ie][1]*TU[ke][1]+TU[ie][2]*TU[ke][2])+2./3.*sm("MASS",24)*(*sus_param).sw2/(1.-(*sus_param).sw2)*(*susy)("H0MIX",ae*10+2)*TU[ie][3]*TU[ke][3];

	// 				for(int je=0;je<2;je++) {
	// 					for(int le=0;le<2;le++) {
	// 						CQ1c+=G1[ie][ke][je][le]/mH0[ae]/mH0[ae]*( sqrt(2.)*(*susy)("H0MIX",ae*10+1)*(*susy)("H0MIX",ae*10+1)*(*sus_param).Mch[je]/sm("MASS",24)/cos(atan((*susy)("HMIX",2)))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/(*sus_param).Mch[je],2.))
	// 						-2.*sqrt(2.)*(*susy)("H0MIX",ae*10+1)/sm("GAUGE",2)*kron(ie,ke)*(Qalj[ae][le][je]*f40(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))+(*sus_param).Mch[je]/(*sus_param).Mch[le]*Qalj[ae][je][le]*f30(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))		+2.*sqrt(2.)*(*susy)("H0MIX",ae*10+1)*T2[ae][ie][ke]*(*sus_param).Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow((*sus_param).Mch[je]/mstop[ke-1],2.))
	// 						+mH0[ae]*mH0[ae]/(*sus_param).Mch[je]/(*sus_param).Mch[je]*kron(ie,ke)*((*susy)("UMIX",20+je)*(*susy)("VMIX",10+le)*f50(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))
	// 						-(*sus_param).Mch[le]/(*sus_param).Mch[je]*(*susy)("UMIX",20+le)*(*susy)("VMIX",10+je)* f60(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))));

	// 						CQ2c+=G1[ie][ke][je][le]/mA0[ae]/mA0[ae]*(sqrt(2.)*(*susy)("A0MIX",ae*10+1)*(*susy)("A0MIX",ae*10+1)*(*sus_param).Mch[je]/sm("MASS",24)/cos(atan((*susy)("HMIX",2)))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/(*sus_param).Mch[je],2.))
	// 						-2.*sqrt(2.)*(*susy)("A0MIX",ae*10+1)/sm("GAUGE",2)*kron(ie,ke)*(-Ralj[ae][le][je]*f40(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))+(*sus_param).Mch[je]/(*sus_param).Mch[le]*Ralj[ae][je][le]*f30(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))			-sqrt(2.)*(*susy)("A0MIX",ae*10+1)*T1[ae][ie][ke]*(*sus_param).mass_top_muW*(*sus_param).Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow((*sus_param).Mch[je]/mstop[ke-1],2.))
	// 						+mA0[ae]*mA0[ae]/(*sus_param).Mch[je]/(*sus_param).Mch[je]*kron(ie,ke)*((*susy)("UMIX",20+je)*(*susy)("VMIX",10+le)*f50(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))
	// 						-(*sus_param).Mch[le]/(*sus_param).Mch[je]*(*susy)("UMIX",20+le)*(*susy)("VMIX",10+je)*f60(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))));
	// 					}		
	// 				}

	// 			}
	// 		}
	// 		CQ1H+=((*susy)("MASS",37)*(*susy)("MASS",37)/sm("MASS",24)/sm("MASS",24)*(*susy)("H0MIX",ae*10+1)*(*susy)("H0MIX",ae*10+1)*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,sm("MASS",24)*sm("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)	+(*sus_param).mass_top_muW*(*sus_param).mass_top_muW*mH0[ae]*mH0[ae]/sm("MASS",24)/sm("MASS",24)/(*susy)("MASS",37)/(*susy)("MASS",37)*f30((*sus_param).mass_top_muW*(*sus_param).mass_top_muW/(*susy)("MASS",37)/(*susy)("MASS",37),(*sus_param).mass_top_muW*(*sus_param).mass_top_muW/sm("MASS",24)/sm("MASS",24)))/mH0[ae]/mH0[ae];
	// 		if (ae < 3) {
	// 			CQ2H+=(((*susy)("MASS",37)*(*susy)("MASS",37)/sm("MASS",24)/sm("MASS",24)*(*susy)("A0MIX",ae*10+1)*(*susy)("A0MIX",ae*10+1)+kron(ae,2)*(*susy)("A0MIX",ae*10+1))*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,sm("MASS",24)*sm("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)	+(*sus_param).mass_top_muW*(*sus_param).mass_top_muW*mA0[ae]*mA0[ae]/sm("MASS",24)/sm("MASS",24)/(*susy)("MASS",37)/(*susy)("MASS",37)*f30((*sus_param).mass_top_muW*(*sus_param).mass_top_muW/(*susy)("MASS",37)/(*susy)("MASS",37),(*sus_param).mass_top_muW*(*sus_param).mass_top_muW/sm("MASS",24)/sm("MASS",24)))/mA0[ae]/mA0[ae];
	// 		}
	// 		for(int je=0;je<2;je++) {
	// 			for(int le=0;le<2;le++) {
	// 				CAc = complex_t(CAc.real(), CAc.imag()+((*susy)("HMIX",2))/sqrt(2.)*G1[ae][ae][je][le]*(v_deltam_s*kron(le,je)*fabs((*sus_param).Mch[je]/sm("MASS",24))*f80(pow(mstop[ae-1]/(*sus_param).Mch[je],2.))-(Ralj[1][je][le]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[le])*f30(pow(mstop[ae-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))-Ralj[1][le][je]*f40(pow(mstop[ae-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))));
	// 				}
	// 			}
	// 	}
	
	// 	CQ1H*=-wilson_p("WPARAM_SI_SM", 3)/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
	// 	CQ2H*=wilson_p("WPARAM_SI_SM", 3)/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
		
	// 	complex_t CAH={0,-lambdaNMSSM*AlambdaNSSM/sm("GAUGE",2)/sm("MASS",24)*(*susy)("HMIX",2)*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,sm("MASS",24)*sm("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)};
	
			
	// 	CQ1c*=wilson_p("WPARAM_SI_SM", 3)/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
	// 	CQ2c*=-wilson_p("WPARAM_SI_SM", 3)/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);		
	
	// 	coeff_temp = (CQ1H+CQ1c)*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2/(*sus_param).epsfac;

    //     this->set_WilsonCoeffMatching("LO", coeff_temp);
	// }

    // // return coeff_temp;

}

void CQ1_susy::NLO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources = {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 1},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 2},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 3},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 9},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 11},
		{ParameterType::WILSON, "WPARAM_SI_SM", 3},
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", 6},
		{ParameterType::WILSON, "EW_SCALE", 1},
		{ParameterType::SM, "MASS", 3},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::SM, "MASS", 23},
		{ParameterType::SM, "GAUGE", 2},
		{ParameterType::SM, "RECKM", 21},
		{ParameterType::SM, "RECKM", 22},
		{ParameterType::BSM, "MASS", 37},
		{ParameterType::BSM, "HMIX", 1},
		{ParameterType::BSM, "HMIX", 2},
		{ParameterType::WILSON, "EPSILON_SUSY", 5}
	};
	
	for (int i = 0; i < 6; ++i) {
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {13, i}});
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14, i}});
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16, i}});
		sources.insert({ParameterType::BSM, "UMIX", i});
		sources.insert({ParameterType::BSM, "UMIX", i+10});
		sources.insert({ParameterType::BSM, "VMIX", i});
		sources.insert({ParameterType::BSM, "VMIX", i+10});
	}
	
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 3; ++j) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, i, j}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {2, i, j}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {9, i, j}});
		}
	}
	
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 6; ++j) {
			for (int k = 0; k < 2; ++k) {
				sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, i, j, k}});
				sources.insert({ParameterType::WILSON, "MATRIX_BSM", {4, i, j, k}});
			}
		}
	}
	
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j) {
			for (int k = 0; k < 3; ++k) {
				sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 0, j, k}});
				sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, i, 1, j, k}});
			}
		}
	}
	
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			sources.insert({ParameterType::SM, "RECKM", i*10 + j});
		}
	}
	
	for (int fe = 0; fe < 3; ++fe) {
		sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}});
	}

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double z = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 1})->get_val();
		double aY = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 11})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();

		double muQ = src.at({ParameterType::BSM, "HMIX", 1})->get_val();
		double xH = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 2})->get_val();
        double xH0 = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 3})->get_val();
		double alpha = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 9})->get_val();
        double beta = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();

		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		double epsfac = src.at({ParameterType::WILSON, "EPSILON_SUSY", 5})->get_val();

		complex_t NQ11H=-src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*(tanb)*tanb/4./mW/mW*(f141(xt,z)+8.*xt*(f30(xt,z)+xt*(f30(xt*1.0001,z)-f30(xt*0.9999,z))/0.0002)*log(Q_match*Q_match/mass_top_muW/mass_top_muW));
		complex_t BQ11H=src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*(tanb)*tanb/4./mW/mW*(f111(xt,z)+8.*(f70(xt*1.0001,z)-f70(xt*0.9999,z))/0.0002*log(Q_match*Q_match/mass_top_muW/mass_top_muW));
		complex_t CQ1H_1=(NQ11H+BQ11H)*mass_b_muW/sw2;
		complex_t CQ2H_1=-CQ1H_1;

		
		complex_t BQ11c1=0.;
		complex_t BQ11c2=0.;
		complex_t NQ11c=0.;
		complex_t NQ21c=0.;
		complex_t BQ11f1=0.;
		complex_t BQ11f2=0.;
		complex_t NQ11f=0.;
		complex_t NQ21f=0.;

		double Dp{0}, Dm{0}, temp{0}, temp2{0};
		double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a0p{0}, a1{0},a2p{0};

		for(int ie=0;ie<2;ie++) {
			for(int ae=0;ae<6;ae++){
				for(int je=0;je<2;je++)  {
					for(int be=0;be<3;be++) {
						BQ11c1+=src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}})->get_val()*(f121(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+4.*(f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))-f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.))));
						BQ11c2+=src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*(f131(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+4.*(f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))-f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.))));

						for(int me=0;me<6;me++){ 
							for(int ne=0;ne<3;ne++) {
								Dp=0.;
								Dm=0.;
								for(int fe=1;fe<=3;fe++) { 	
									Dp+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, fe}})->get_val()+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}})->get_val());
									Dm+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, fe}})->get_val()-src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}})->get_val());
								}
								a0a=-(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val())*(f181(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))+4.*(f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.)*1.0001)-f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.)))*src.at({ParameterType::BSM, "UMIX", ie*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX", je*10+0})->get_val())*kron(ae,me);
								a0b=-((f191(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))+4.*(f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.)*1.0001)-f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.)))*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*src.at({ParameterType::BSM, "UMIX", ie*10+0})->get_val())*kron(ae,me);
								a0c=1./mW*(f171(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+4.*(f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+(f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))-f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)))/0.0002+(f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001)-f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999))/0.0002)*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.)))*kron(ie,je);
							
								a0Q1=a0a+a0b+Dp*a0c;
								a0Q2=-a0a+a0b+Dm*a0c;
								a0p=4.*src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()/mW/(src.at({ParameterType::SM, "RECKM", be*10+2})->get_val()*src.at({ParameterType::SM, "RECKM", ne*10+1})->get_val()/src.at({ParameterType::SM, "RECKM", 22})->get_val()/src.at({ParameterType::SM, "RECKM", 21})->get_val())/src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*f151(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
								a1=src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/sqrt(2.)/mW*(f161(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+4.*(f80(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001)-f80(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.)))*kron(ie,je)*kron(ae,me);
								a2p=src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*(src.at({ParameterType::SM, "RECKM", be*10+2})->get_val()*src.at({ParameterType::SM, "RECKM", ne*10+1})->get_val()/src.at({ParameterType::SM, "RECKM", 22})->get_val()/src.at({ParameterType::SM, "RECKM", 21})->get_val())*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()/2./mW*f151(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
								
								NQ11c+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*(a0Q1+a1*tanb)
								+src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*a0p
								+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*a2p*pow(src.at({ParameterType::SM, "MASS", 3})->get_val()*tanb,2.);	
								NQ21c+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*(a0Q2+a1*tanb)
								+src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*a0p
								+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*a2p*pow(src.at({ParameterType::SM, "MASS", 3})->get_val()*tanb,2.);
							}
						}
					
					}
					for(int be=0;be<6;be++) {
						for(int ce=0;ce<6;ce++) {
							for(int fe=0;fe<3;fe++) {
								BQ11f1+=-src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()*pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ae, ce}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ce, be}})->get_val()*(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val(),2.)))	*(f90(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*src.at({ParameterType::WILSON, "MATRIX_BSM", {6, ie, fe, 2}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, fe, 2}})->get_val());
								BQ11f2+=-src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()*pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ae, ce}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ce, be}})->get_val()*(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val(),2.)))	*(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*f100(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, fe, 1}})->get_val());

							}
						}
					}
				}

				for(int me=0;me<3;me++) {
					for(int ne=0;ne<3;ne++) {
						for(int de=0;de<6;de++) {
							for(int ke=0;ke<6;ke++) {
								
								temp2 = src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, me}})->get_val()*src.at({ParameterType::BSM, "UMIX", ie*10+2})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ae, ke}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ke, de}})->get_val()*(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(),2.)))*tanb*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/sqrt(2.)*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.));
								NQ11f+=temp2;
								NQ21f+=temp2;
								for(int ce=0;ce<6;ce++) {
									Dp=0.;
									Dm=0.;
									for(int fe=0;fe<3;fe++) 
									{		
										Dp+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, de, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val()+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ce, fe}})->get_val()); 
										Dm+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, de, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val()-src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ce, fe}})->get_val()); 
									}
									temp=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, me}})->get_val()*src.at({ParameterType::BSM, "UMIX", ie*10+2})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ae, ke}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ke, de}})->get_val()*
									(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(),2.)))*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.));	
									NQ11f+=Dp*temp;
									NQ21f+=Dm*temp;
								}
								for(int je=0;je<2;je++) {
									temp=-src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, me}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+2})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ae, ke}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ke, de}})->get_val()*
									(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(),2.)))*mW*(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val())*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX", ie*10+2})->get_val()*src.at({ParameterType::BSM, "VMIX", je*10+1})->get_val()+
									f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX", je*10+2})->get_val()*src.at({ParameterType::BSM, "VMIX", ie*10+1})->get_val()); 
									NQ11f+=temp;
									NQ21f+=-temp;
								}

							}
						}
					}
						
				}
			}
		}
		
		complex_t BQ11c=(BQ11c1+BQ11c2)*kappa*mW*mW/2./g2/g2/sw2;
		complex_t BQ21c=-(BQ11c1-BQ11c2)*kappa*mW*mW/2./g2/g2/sw2;


		NQ11c*=src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
		NQ21c*=-src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
		
		complex_t CQ1charg_1=NQ11c+BQ11c;
		


		
		complex_t CQ2charg_1=NQ21c+BQ21c;

		complex_t BQ11f=(BQ11f1+BQ11f2)*2./3.*kappa/g2/g2/sw2;
		complex_t BQ21f=-(BQ11f1-BQ11f2)*2./3.*kappa/g2/g2/sw2;
		
		
		NQ11f*=-4./3.*src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*(tanb)*tanb/mW/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;

		NQ21f*=4./3.*src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*(tanb)*tanb/mW/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;


		complex_t CQ1four_1=NQ11f+BQ11f;
		complex_t coeff_temp = (CQ1H_1+CQ1charg_1)/epsfac+CQ1four_1;

        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3230, 1, 1)}, sources, func);

	// // sus_param->reset_PrimeCQG(this->get_Q_match());
	// complex_t NQ11H=-wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./sm("MASS",24)/sm("MASS",24)*(f141((*sus_param).xt,(*sus_param).z)+8.*(*sus_param).xt*(f30((*sus_param).xt,(*sus_param).z)+(*sus_param).xt*(f30((*sus_param).xt*1.0001,(*sus_param).z)-f30((*sus_param).xt*0.9999,(*sus_param).z))/0.0002)*log(this->get_Q_match()*this->get_Q_match()/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW));
	// complex_t BQ11H=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./sm("MASS",24)/sm("MASS",24)*(f111((*sus_param).xt,(*sus_param).z)+8.*(f70((*sus_param).xt*1.0001,(*sus_param).z)-f70((*sus_param).xt*0.9999,(*sus_param).z))/0.0002*log(this->get_Q_match()*this->get_Q_match()/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW));
	// complex_t CQ1H_1=(NQ11H+BQ11H)*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;
	// complex_t CQ2H_1=-CQ1H_1;

	
	// complex_t BQ11c1=0.;
	// complex_t BQ11c2=0.;
	// complex_t NQ11c=0.;
	// complex_t NQ21c=0.;
	// complex_t BQ11f1=0.;
	// complex_t BQ11f2=0.;
	// complex_t NQ11f=0.;
	// complex_t NQ21f=0.;

	// double Dp{0}, Dm{0}, temp{0}, temp2{0};
	// double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a0p{0}, a1{0},a2p{0};

	// for(int ie=0;ie<2;ie++) {
	// 	for(int ae=0;ae<6;ae++){
	// 		for(int je=0;je<2;je++)  {
	// 			for(int be=0;be<3;be++) {
	// 				BQ11c1+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][1]*(*sus_param).X_NL[je][be][1]*(f121(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))+4.*(f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))-f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.))));
	// 				BQ11c2+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][1]*(*sus_param).X_NR[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*(f131(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))+4.*(f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))-f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.))));

	// 				for(int me=0;me<6;me++){ 
	// 					for(int ne=0;ne<3;ne++) {
	// 						Dp=0.;
	// 						Dm=0.;
	// 						for(int fe=1;fe<=3;fe++) { 	
	// 							Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
	// 							Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
	// 						}
	// 						a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*(f181(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))+4.*(f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*1.0001)-f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,me);
	// 						a0b=-((f191(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))+4.*(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*1.0001)-f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*(*susy)("UMIX", je*10+1)*(*susy)("UMIX", ie*10+0))*kron(ae,me);
	// 						a0c=1./sm("MASS",24)*(f171(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))+4.*(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))+(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))-f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)))/0.0002+(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)*1.0001)-f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)*0.9999))/0.0002)*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*kron(ie,je);
						
	// 						a0Q1=a0a+a0b+Dp*a0c;
	// 						a0Q2=-a0a+a0b+Dm*a0c;
	// 						a0p=4.*(*sus_param).G_aimn[ae][ie][be][ne]/sm("MASS",24)/(src.at({ParameterType::SM, "RECKM", be*10+2})->get_val()*src.at({ParameterType::SM, "RECKM", ne*10+1})->get_val()/src.at({ParameterType::SM, "RECKM", 22})->get_val()/src.at({ParameterType::SM, "RECKM", 21})->get_val())/(*susy)("UMIX", je*10+1)*f151(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
	// 						a1=(*sus_param).Mch[ie]/sqrt(2.)/sm("MASS",24)*(f161(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))+4.*(f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001)-f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*kron(ie,je)*kron(ae,me);
	// 						a2p=(*sus_param).Gamma_UL[me][be]*(src.at({ParameterType::SM, "RECKM", be*10+2})->get_val()*src.at({ParameterType::SM, "RECKM", ne*10+1})->get_val()/src.at({ParameterType::SM, "RECKM", 22})->get_val()/src.at({ParameterType::SM, "RECKM", 21})->get_val())*(*susy)("UMIX", je*10+1)/2./sm("MASS",24)*f151(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
							
	// 						NQ11c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q1+a1*(*susy)("HMIX",2))
	// 						+(*sus_param).G_aimn[ae][ie][be][ne]*(*susy)("UMIX", je*10+1)*a0p
	// 						+(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*a2p*pow(sm("MASS",3)*(*susy)("HMIX",2),2.);	
	// 						NQ21c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q2+a1*(*susy)("HMIX",2))
	// 						+(*sus_param).G_aimn[ae][ie][be][ne]*(*susy)("UMIX", je*10+1)*a0p
	// 						+(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*a2p*pow(sm("MASS",3)*(*susy)("HMIX",2),2.);
	// 					}
	// 				}
				
	// 			}
	// 			for(int be=0;be<6;be++) {
	// 				for(int ce=0;ce<6;ce++) {
	// 					for(int fe=0;fe<3;fe++) {
	// 						BQ11f1+=-(*sus_param).X_UL[je][be][1]*(*sus_param).X_UR[ie][ae][2]*pow(sm("MASS",24)/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][ce]*(*sus_param).MsqU[ce]/(*sus_param).Mch[ie]*(*sus_param).P_U[ce][be]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ce],2.)))	*(f90(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[fe]/(*sus_param).Mch[ie],2.))*(*sus_param).X_NR[ie][fe][2]*(*sus_param).X_NL[je][fe][2]);
	// 						BQ11f2+=-(*sus_param).X_UL[je][be][1]*(*sus_param).X_UR[ie][ae][2]*pow(sm("MASS",24)/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][ce]*(*sus_param).MsqU[ce]/(*sus_param).Mch[ie]*(*sus_param).P_U[ce][be]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ce],2.)))	*(fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f100(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[fe]/(*sus_param).Mch[ie],2.))*(*sus_param).X_NL[ie][fe][1]*(*sus_param).X_NR[je][fe][1]);

	// 					}
	// 				}
	// 			}
	// 		}

	// 		for(int me=0;me<3;me++) {
	// 			for(int ne=0;ne<3;ne++) {
	// 				for(int de=0;de<6;de++) {
	// 					for(int ke=0;ke<6;ke++) {
							
	// 						temp2 = (*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[de][me]*(*susy)("UMIX",ie*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[ie]*(*sus_param).P_U[ke][de]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*(*susy)("HMIX",2)*(*sus_param).Mch[ie]/sqrt(2.)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2.));
	// 						NQ11f+=temp2;
	// 						NQ21f+=temp2;
	// 						for(int ce=0;ce<6;ce++) {
	// 							Dp=0.;
	// 							Dm=0.;
	// 							for(int fe=0;fe<3;fe++) 
	// 							{		
	// 								Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[de][fe]*(*sus_param).Gamma_UL[ce][fe]+(*sus_param).Gamma_UL[de][fe]*(*sus_param).Gamma_UR[ce][fe]); 
	// 								Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[de][fe]*(*sus_param).Gamma_UL[ce][fe]-(*sus_param).Gamma_UL[de][fe]*(*sus_param).Gamma_UR[ce][fe]); 
	// 							}
	// 							temp=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[ce][me]*(*susy)("UMIX",ie*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[ie]*(*sus_param).P_U[ke][de]*
	// 							(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*f60(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.));	
	// 							NQ11f+=Dp*temp;
	// 							NQ21f+=Dm*temp;
	// 						}
	// 						for(int je=0;je<2;je++) {
	// 							temp=-(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[de][me]*(*susy)("UMIX",je*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[je]*(*sus_param).P_U[ke][de]*
	// 							(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*sm("MASS",24)*(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f60(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[je],2.))*(*susy)("UMIX",ie*10+2)*(*susy)("VMIX",je*10+1)+
	// 							f50(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[je],2.))*(*susy)("UMIX",je*10+2)*(*susy)("VMIX",ie*10+1)); 
	// 							NQ11f+=temp;
	// 							NQ21f+=-temp;
	// 						}

	// 					}
	// 				}
	// 			}
					
	// 		}
	// 	}
	// }
	
	// complex_t BQ11c=(BQ11c1+BQ11c2)*(*sus_param).kappa*sm("MASS",24)*sm("MASS",24)/2./sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;
	// complex_t BQ21c=-(BQ11c1-BQ11c2)*(*sus_param).kappa*sm("MASS",24)*sm("MASS",24)/2./sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;


	// NQ11c*=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;
	// NQ21c*=-wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;
	
	// complex_t CQ1charg_1=NQ11c+BQ11c;
	


	
	// complex_t CQ2charg_1=NQ21c+BQ21c;

	// complex_t BQ11f=(BQ11f1+BQ11f2)*2./3.*(*sus_param).kappa/sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;
	// complex_t BQ21f=-(BQ11f1-BQ11f2)*2./3.*(*sus_param).kappa/sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;
	
	
	// NQ11f*=-4./3.*wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;

	// NQ21f*=4./3.*wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;


	// complex_t CQ1four_1=NQ11f+BQ11f;
	// complex_t coeff_temp = (CQ1H_1+CQ1charg_1)/sus_param->epsfac+CQ1four_1;

    // this->set_WilsonCoeffMatching("NLO", coeff_temp);
	// // return coeff_temp;
}

void CQ2_susy::LO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources = {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 1},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 2},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 3},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 4},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 9},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 11},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
		{ParameterType::WILSON, "EW_SCALE", 1},
		{ParameterType::BSM, "MASS", 37},
		{ParameterType::BSM, "HMIX", 1},
		{ParameterType::BSM, "HMIX", 2},
		{ParameterType::WILSON, "EPSILON_SUSY", 5},
		{ParameterType::WILSON, "WPARAM_SI_SM", 3},
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},
		{ParameterType::WILSON, "WPARAM_SI_SM", 5},
		{ParameterType::WILSON, "WPARAM_RUN_SM", 1},
		{ParameterType::SM, "MASS", 23},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::SM, "GAUGE", 2},
		{ParameterType::SM, "QCD", LhaID(5,1)},
		{ParameterType::SM, "SMINPUTS", 2}
	};
	
	// Ajouter les éléments qui dépendent de boucles :
	
	// MATRIX_BSM {3, je, ae, 1} , {4, ie, ae, 2} , {5, ie, be, 1}, {6, je, be, 1}
	for (int je = 0; je < 2; ++je) {
		for (int ae = 0; ae < 6; ++ae) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {4, je, ae, 2}});
		}
	}
	for (int ie = 0; ie < 2; ++ie) {
		for (int be = 0; be < 3; ++be) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}});
		}
	}
	
	// WPARAM_SI_BSM {13, ie}, {14, ae}, {16, be}
	for (int ie = 0; ie < 2; ++ie) {
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}});
	}
	for (int ae = 0; ae < 6; ++ae) {
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}});
	}
	for (int be = 0; be < 3; ++be) {
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}});
	}
	
	for (int i = 0; i < 6; ++i) {
		sources.insert({ParameterType::BSM, "UMIX", i});
		sources.insert({ParameterType::BSM, "UMIX", i+10});
		sources.insert({ParameterType::BSM, "VMIX", i});
		sources.insert({ParameterType::BSM, "VMIX", i+10});
	}
	
	sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}});
	sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}});
	sources.insert({ParameterType::WILSON, "WPARAM_MATCH_SM", 6});
	
	for (int fe = 0; fe < 3; ++fe) {
		sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}});
	}
	
	for (int ae = 0; ae < 6; ++ae) {
		for (int fe = 0; fe < 3; ++fe) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}});
		}
	}
	for (int me = 0; me < 6; ++me) {
		for (int be = 0; be < 3; ++be) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}});
		}
	}
	
	for (int ae = 0; ae < 6; ++ae) {
		for (int ie = 0; ie < 2; ++ie) {
			for (int be = 0; be < 3; ++be) {
				for (int ne = 0; ne < 3; ++ne) {
					sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, ae, ie, be, ne}});
				}
			}
		}
	}

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double z = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 1})->get_val();
		double aY = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 11})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();

		double muQ = src.at({ParameterType::BSM, "HMIX", 1})->get_val();
		double xH = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 2})->get_val();
        double xH0 = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 3})->get_val();
		double xA = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 4})->get_val();
		double alpha = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 9})->get_val();
        double beta = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();

		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		double epsfac = src.at({ParameterType::WILSON, "EPSILON_SUSY", 5})->get_val();

		
		complex_t BQ10c1=0.;
		complex_t BQ10c2=0.;

		double Dp, Dm;
		double a0a{0}, a0b{0}, a0c, a0Q1{0}, a0Q2{0};
		double a1{0};
		complex_t NQ20c=0.;

		for(int ie=0;ie<2;ie++) {
			for(int je=0;je<2;je++) {
				for(int ae=0;ae<6;ae++) {
					for(int be=0;be<3;be++) { 
						BQ10c1+=src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {6,ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}})->get_val()*f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));
						BQ10c2+=src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4,ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)));

						for(int me=0;me<6;me++) {
							for(int ne=0;ne<3;ne++) {
								Dp=0.;
								Dm=0.;
								for(int fe=0;fe<3;fe++) 
								{
									Dp+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, fe}})->get_val()+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}})->get_val());
									Dm+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, fe}})->get_val()-src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}})->get_val());
								}
								a0a=-(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val())*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX", ie*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX", je*10+0})->get_val())*kron(ae,me);
								a0b=-(f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX", ie*10+0})->get_val())*kron(ae,me);
								a0c=1./mW*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je);
								a0Q2=-a0a+a0b+Dm*a0c;
								
								a1=src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/sqrt(2.)/mW*f80(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je)*kron(ae,me);
								NQ20c+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*(a0Q2+a1*tanb);
								// LOG_INFO("(*sus_param).G_aimn",ae, ie, be, ne, src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val());
							}
						
						}
					}
				}
			}
		}
		complex_t BQ20c=-(BQ10c1-BQ10c2)*kappa*mW*mW/2./g2/g2/sw2;
		
		double le = -tanb;
		double G3=ld*(ld*lu+1.)*F6SP(xt,xH)+ld*lu*lu*F7SP(xt,xH)
		+lu*lu*(ld*F8SP(xt,xH)+lu*F9SP(xt,xH)+lu*F10SP(xt,xH))+lu*F11SP(xt,xH)+lu*F12SP(xt,xH);
		double CPn_2HDM=xt*(-le*(ld*F1SP(xt,xH)+lu*F2SP(xt,xH))+le*lu*F3SP(xt,xH))+xt/2./xA*(le)*G3;

		double CQ2H_0=CPc_2HDM(xH,xt,lu,ld,le,sw2)+CPn_2HDM;
		CQ2H_0*=(src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*mass_b_muW/mW/mW)/sw2;

		NQ20c*=-src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;

		complex_t CQ2charg_0=NQ20c+BQ20c;


		complex_t coeff_temp = (CQ2charg_0+CQ2H_0)/epsfac;

		/* NMSSM */

		double lambdaNMSSM = 1;
		double lambdaSNMSSM = 1;
		double AlambdaNSSM = 1;
		double kappaNMSSM = 1;
		double m_Bs = 1;
		double mass_nutl = 1;

		if(src.at({ParameterType::BSM, "MASS",  46})->get_val()!=0.||src.at({ParameterType::BSM, "MASS",  45})->get_val()!=0.)
		{
			LOG_INFO("NMSSM ? Doesn't exist, don't search for it.");

			double s=lambdaSNMSSM/lambdaNMSSM;
			double v=sqrt(1./sqrt(2.)/src.at({ParameterType::SM, "SMINPUTS",  2})->get_val());
			
			double v_deltam_s=v/s*(sqrt(2.)*AlambdaNSSM-2.*kappaNMSSM*s)/(sqrt(2.)*AlambdaNSSM+kappaNMSSM*s);
			
			double mH0[4],mA0[3],mstop[3];
			
			mstop[0]=src.at({ParameterType::BSM, "MASS",  2000013})->get_val(); //mass upr, is that right ?
			mstop[1]=src.at({ParameterType::BSM, "MASS",  1000006})->get_val();
			mstop[2]=src.at({ParameterType::BSM, "MASS",  2000006})->get_val();
			
			mH0[1]=src.at({ParameterType::BSM, "MASS",  25})->get_val();
			mH0[2]=src.at({ParameterType::BSM, "MASS",  35})->get_val();
			mH0[3]=src.at({ParameterType::BSM, "MASS",  36})->get_val();
			mA0[1]=src.at({ParameterType::BSM, "MASS",  36})->get_val();
			mA0[2]=src.at({ParameterType::BSM, "MASS",  36})->get_val();
			
			double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
			double T1[3][4][4],T2[4][4][4];
			std::array<std::array<double,4>,4> TU;
		
			TU[1][1]=1.;
			for(int ie=0;ie<2;ie++){
				for(int je=0;je<2;je++) {
					TU[ie+1][je+1]=src.at({ParameterType::BSM, "STOPMIX",  ie*10+je})->get_val();
				}
			}

			
			double vu=sqrt(pow(sin(atan(tanb)),2.)/sqrt(2.)/src.at({ParameterType::SM, "SMINPUTS",  2})->get_val());
			double vd=vu/tanb;

			for(int je=0;je<2;je++) {
				for(int le=0;le<2;le++) {
					for(int ae=0;ae<3;ae++) {
						if (ae <3 ){
							Ralj[ae][le][je]=-g2/sqrt(2.)*(src.at({ParameterType::BSM, "A0MIX",  ae*10+1})->get_val()*src.at({ParameterType::BSM, "UMIX",20+le})->get_val()*src.at({ParameterType::BSM, "VMIX",20+je})->get_val()+src.at({ParameterType::BSM, "A0MIX", ae*10+2})->get_val()*src.at({ParameterType::BSM, "UMIX", 10+le})->get_val()*src.at({ParameterType::BSM, "VMIX",20+je})->get_val())-lambdaNMSSM/sqrt(2.)*src.at({ParameterType::BSM, "A0MIX", ae*10+3})->get_val()*src.at({ParameterType::BSM, "UMIX",20+le})->get_val()*src.at({ParameterType::BSM, "VMIX",20+je})->get_val();
						}
						Qalj[ae][le][je]=g2/sqrt(2.)*(src.at({ParameterType::BSM, "H0MIX",  ae*10+1})->get_val()*src.at({ParameterType::BSM, "UMIX",20+le})->get_val()*src.at({ParameterType::BSM, "VMIX",20+je})->get_val()+src.at({ParameterType::BSM, "H0MIX",  ae*10+2})->get_val()*src.at({ParameterType::BSM, "UMIX", 10+le})->get_val()*src.at({ParameterType::BSM, "VMIX",20+je})->get_val())-lambdaNMSSM/sqrt(2.)*src.at({ParameterType::BSM, "H0MIX", ae*10+3})->get_val()*src.at({ParameterType::BSM, "UMIX",20+le})->get_val()*src.at({ParameterType::BSM, "VMIX",20+je})->get_val();
						for(int ke=1;ke<=3;ke++) {
							G1[ae][ke][je][le]=(TU[ae][2]*TU[ke][2]-kron(ae,1)*kron(ke,1))*src.at({ParameterType::BSM, "VMIX",10+le})->get_val()*src.at({ParameterType::BSM, "UMIX",20+je})->get_val()-mass_top_muW/sqrt(2.)/sin(atan(tanb))/mW*TU[ae][3]*TU[ke][2]*src.at({ParameterType::BSM, "VMIX",20+le})->get_val()*src.at({ParameterType::BSM, "UMIX",20+je})->get_val();
						}
					}
				}
				for(int ie=0;ie<3;ie++) {
					for(int ke=0;ke<3;ke++) {
						T1[je][ie][ke]=(TU[ie][3]*TU[ke][2]-TU[ie][2]*TU[ke][3])*((lambdaNMSSM/sqrt(2.)*(vd*src.at({ParameterType::BSM, "A0MIX", je*10+3})->get_val()+s*src.at({ParameterType::BSM, "A0MIX", je*10+1})->get_val()))-src.at({ParameterType::BSM, "AU",  11})->get_val()*src.at({ParameterType::BSM, "A0MIX", je*10+2})->get_val());
					}
				}
			}

			complex_t CQ1H=0.;
			complex_t CQ2H=0.;
			complex_t CQ1c=0.;
			complex_t CQ2c=0.;
			complex_t CAc=0.;

			for(int ae=0;ae<3;ae++) {
				for(int ie=0;ie<3;ie++) {
					for(int ke=0;ke<3;ke++){
						T2[ae][ie][ke]=-mass_top_muW/2./mW*(2.*mass_top_muW*src.at({ParameterType::BSM, "H0MIX",  ae*10+2})->get_val()*(TU[ie][2]*TU[ke][2]+TU[ie][3]*TU[ke][3])	+((lambdaNMSSM/sqrt(2.)*(vd*src.at({ParameterType::BSM, "H0MIX", ae*10+3})->get_val()+s*src.at({ParameterType::BSM, "H0MIX",  ae*10+1})->get_val()))+src.at({ParameterType::BSM, "AU",  11})->get_val()*src.at({ParameterType::BSM, "H0MIX",  ae*10+2})->get_val())*(TU[ie][3]*TU[ke][2]+TU[ie][2]*TU[ke][3]))	+src.at({ParameterType::SM, "MASS",  23})->get_val()/2./sqrt(1.-sw2)*(1.-4./3.*sw2)*src.at({ParameterType::BSM, "H0MIX",  ae*10+2})->get_val()*(TU[ie][1]*TU[ke][1]+TU[ie][2]*TU[ke][2])+2./3.*mW*sw2/(1.-sw2)*src.at({ParameterType::BSM, "H0MIX",  ae*10+2})->get_val()*TU[ie][3]*TU[ke][3];

						for(int je=0;je<2;je++) {
							for(int le=0;le<2;le++) {
								CQ1c+=G1[ie][ke][je][le]/mH0[ae]/mH0[ae]*( sqrt(2.)*src.at({ParameterType::BSM, "H0MIX",  ae*10+1})->get_val()*src.at({ParameterType::BSM, "H0MIX",  ae*10+1})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mW/cos(atan(tanb))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))
								-2.*sqrt(2.)*src.at({ParameterType::BSM, "H0MIX",  ae*10+1})->get_val()/g2*kron(ie,ke)*(Qalj[ae][le][je]*f40(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))+src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()*Qalj[ae][je][le]*f30(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.)))		+2.*sqrt(2.)*src.at({ParameterType::BSM, "H0MIX",  ae*10+1})->get_val()*T2[ae][ie][ke]*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mstop[ke-1],2.))
								+mH0[ae]*mH0[ae]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*kron(ie,ke)*(src.at({ParameterType::BSM, "UMIX",20+je})->get_val()*src.at({ParameterType::BSM, "VMIX",10+le})->get_val()*f50(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(mass_nutl/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))
								-src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*src.at({ParameterType::BSM, "UMIX",20+le})->get_val()*src.at({ParameterType::BSM, "VMIX",10+je})->get_val()* f60(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(mass_nutl/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))));

								CQ2c+=G1[ie][ke][je][le]/mA0[ae]/mA0[ae]*(sqrt(2.)*src.at({ParameterType::BSM, "A0MIX",  ae*10+1})->get_val()*src.at({ParameterType::BSM, "A0MIX",  ae*10+1})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mW/cos(atan(tanb))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))
								-2.*sqrt(2.)*src.at({ParameterType::BSM, "A0MIX",  ae*10+1})->get_val()/g2*kron(ie,ke)*(-Ralj[ae][le][je]*f40(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))+src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()*Ralj[ae][je][le]*f30(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.)))			-sqrt(2.)*src.at({ParameterType::BSM, "A0MIX",  ae*10+1})->get_val()*T1[ae][ie][ke]*mass_top_muW*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mstop[ke-1],2.))
								+mA0[ae]*mA0[ae]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*kron(ie,ke)*(src.at({ParameterType::BSM, "UMIX",20+je})->get_val()*src.at({ParameterType::BSM, "VMIX",10+le})->get_val()*f50(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(mass_nutl/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))
								-src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*src.at({ParameterType::BSM, "UMIX",20+le})->get_val()*src.at({ParameterType::BSM, "VMIX",10+je})->get_val()*f60(pow(mstop[ie-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(mass_nutl/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))));
							}		
						}

					}
				}
				CQ1H+=(mH*mH/mW/mW*src.at({ParameterType::BSM, "H0MIX",  ae*10+1})->get_val()*src.at({ParameterType::BSM, "H0MIX",  ae*10+1})->get_val()*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)	+mass_top_muW*mass_top_muW*mH0[ae]*mH0[ae]/mW/mW/mH/mH*f30(mass_top_muW*mass_top_muW/mH/mH,mass_top_muW*mass_top_muW/mW/mW))/mH0[ae]/mH0[ae];
				if (ae < 3) {
					CQ2H+=((mH*mH/mW/mW*src.at({ParameterType::BSM, "A0MIX",  ae*10+1})->get_val()*src.at({ParameterType::BSM, "A0MIX",  ae*10+1})->get_val()+kron(ae,2)*src.at({ParameterType::BSM, "A0MIX",  ae*10+1})->get_val())*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)	+mass_top_muW*mass_top_muW*mA0[ae]*mA0[ae]/mW/mW/mH/mH*f30(mass_top_muW*mass_top_muW/mH/mH,mass_top_muW*mass_top_muW/mW/mW))/mA0[ae]/mA0[ae];
				}
				for(int je=0;je<2;je++) {
					for(int le=0;le<2;le++) {
						CAc = complex_t(CAc.real(), CAc.imag()+(tanb)/sqrt(2.)*G1[ae][ae][je][le]*(v_deltam_s*kron(le,je)*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/mW)*f80(pow(mstop[ae-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))-(Ralj[1][je][le]*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val())*f30(pow(mstop[ae-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.))-Ralj[1][le][je]*f40(pow(mstop[ae-1]/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, le}})->get_val(),2.)))));
						}
					}
			}
		
			CQ1H*=-src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()/4.*tanb*tanb;
			CQ2H*=src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()/4.*tanb*tanb;
			
			complex_t CAH={0,-lambdaNMSSM*AlambdaNSSM/g2/mW*tanb*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)};
		
				
			CQ1c*=src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()/4.*tanb*tanb;
			CQ2c*=-src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()/4.*tanb*tanb;		
		
			coeff_temp = (CQ2H+CQ2c)*mass_b_muW/sw2/epsfac;

			complex_t CA=CAH+CAc;

			if(src.at({ParameterType::BSM, "MASS",  36})->get_val()>Q_match) coeff_temp+=-v_deltam_s/2.*mass_b_muW/sw2*src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*CA/src.at({ParameterType::BSM, "MASS",  36})->get_val()/src.at({ParameterType::BSM, "MASS",  36})->get_val();

			
			if((src.at({ParameterType::BSM, "MASS",  36})->get_val()>(*Parameters::GetInstance())("QCD", LhaID(5, 1)))&&(src.at({ParameterType::BSM, "MASS",  36})->get_val()<Q_match ))
			{	
				double alphas_Ma1  = QCDHelper::alpha_s(src.at({ParameterType::BSM, "MASS",  36})->get_val());	
				double mass_b_ma1=QCDHelper::msbar_mass(5, src.at({ParameterType::BSM, "MASS",  36})->get_val());
				coeff_temp+=-v_deltam_s/2.*mass_b_ma1/sw2*src.at({ParameterType::WILSON, "WPARAM_SI_SM",  3})->get_val()*CA/src.at({ParameterType::BSM, "MASS",  36})->get_val()/src.at({ParameterType::BSM, "MASS",  36})->get_val() *pow(alphas_Ma1/src.at({ParameterType::WILSON, "WPARAM_RUN_SM", 1})->get_val(),-4./src.at({ParameterType::WILSON, "WPARAM_SI_SM", 5})->get_val());
			}
			

		}

        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3233, 0, 1)}, sources, func);


	// // sus_param->reset_PrimeCQG(this->get_Q_match());
    // complex_t BQ10c1=0.;
	// complex_t BQ10c2=0.;

	// double Dp, Dm;
	// double a0a{0}, a0b{0}, a0c, a0Q1{0}, a0Q2{0};
	// double a1{0};
	// complex_t NQ20c=0.;

	// for(int ie=0;ie<2;ie++) {
	// 	for(int je=0;je<2;je++) {
	// 		for(int ae=0;ae<6;ae++) {
	// 			for(int be=0;be<3;be++) { 
	// 				BQ10c1+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][1]*(*sus_param).X_NL[je][be][1]*f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));
	// 				BQ10c2+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][1]*(*sus_param).X_NR[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)));

	// 				for(int me=0;me<6;me++) {
	// 					for(int ne=0;ne<3;ne++) {
	// 						Dp=0.;
	// 						Dm=0.;
	// 						for(int fe=0;fe<3;fe++) 
	// 						{
	// 							Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
	// 							Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
	// 						}
	// 						a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,me);
	// 						a0b=-(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))*(*susy)("UMIX", je*10+1)*(*susy)("VMIX", ie*10+0))*kron(ae,me);
	// 						a0c=1./sm("MASS",24)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))*kron(ie,je);
	// 						a0Q2=-a0a+a0b+Dm*a0c;
							
	// 						a1=(*sus_param).Mch[ie]/sqrt(2.)/sm("MASS",24)*f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me);
	// 						NQ20c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q2+a1*(*susy)("HMIX",2));
	// 						// LOG_INFO("(*sus_param).G_aimn",ae, ie, be, ne, (*sus_param).G_aimn[ae][ie][be][ne]);
	// 					}
					
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// complex_t BQ20c=-(BQ10c1-BQ10c2)*(*sus_param).kappa*sm("MASS",24)*sm("MASS",24)/2./sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;
	
	// double le = -(*susy)("HMIX", 2);
    // double G3=sus_param->ld*(sus_param->ld*sus_param->lu+1.)*F6SP(sus_param->xt,sus_param->xH)+sus_param->ld*sus_param->lu*sus_param->lu*F7SP(sus_param->xt,sus_param->xH)
	// +sus_param->lu*sus_param->lu*(sus_param->ld*F8SP(sus_param->xt,sus_param->xH)+sus_param->lu*F9SP(sus_param->xt,sus_param->xH)+sus_param->lu*F10SP(sus_param->xt,sus_param->xH))+sus_param->lu*F11SP(sus_param->xt,sus_param->xH)+sus_param->lu*F12SP(sus_param->xt,sus_param->xH);
    // double CPn_2HDM=sus_param->xt*(-le*(sus_param->ld*F1SP(sus_param->xt,sus_param->xH)+sus_param->lu*F2SP(sus_param->xt,sus_param->xH))+le*sus_param->lu*F3SP(sus_param->xt,sus_param->xH))+sus_param->xt/2./sus_param->xA*(le)*G3;

    // double CQ2H_0=CPc_2HDM(sus_param->xH,sus_param->xt,sus_param->lu,sus_param->ld,le,sus_param->sw2)+CPn_2HDM;
    // CQ2H_0*=(wilson_p("WPARAM_SI_SM", 3)*sus_param->mass_b_muW/sm("MASS",24)/sm("MASS",24))/sus_param->sw2;

	// NQ20c*=-wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;

	// complex_t CQ2charg_0=NQ20c+BQ20c;


    // complex_t coeff_temp = (CQ2charg_0+CQ2H_0)/sus_param->epsfac;
    // this->set_WilsonCoeffMatching("LO", coeff_temp);

    // /* NMSSM */

	// double lambdaNMSSM = 1;
	// double lambdaSNMSSM = 1;
	// double AlambdaNSSM = 1;
	// double kappaNMSSM = 1;
	// double m_Bs = 1;
	// double mass_nutl = 1;

    // if((*susy)("MASS",46)!=0.||(*susy)("MASS",45)!=0.)
	// {
	// 	LOG_INFO("NMSSM ? Doesn't exist, don't search for it.");

	// 	double s=lambdaSNMSSM/lambdaNMSSM;
	// 	double v=sqrt(1./sqrt(2.)/sm("SMINPUTS", 2));
		
	// 	double v_deltam_s=v/s*(sqrt(2.)*AlambdaNSSM-2.*kappaNMSSM*s)/(sqrt(2.)*AlambdaNSSM+kappaNMSSM*s);
		
	// 	double mH0[4],mA0[3],mstop[3];
		
	// 	mstop[0]=(*susy)("MASS", 2000013); //mass upr, is that right ?
	// 	mstop[1]=(*susy)("MASS", 1000006);
	// 	mstop[2]=(*susy)("MASS", 2000006);
		
	// 	mH0[1]=(*susy)("MASS", 25);
	// 	mH0[2]=(*susy)("MASS", 35);
	// 	mH0[3]=(*susy)("MASS", 36);
	// 	mA0[1]=(*susy)("MASS",36);
	// 	mA0[2]=(*susy)("MASS",36);
		
	// 	double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
	// 	double T1[3][4][4],T2[4][4][4];
	// 	std::array<std::array<double,4>,4> TU;
	
	// 	TU[1][1]=1.;
	// 	for(int ie=0;ie<2;ie++){
	// 		 for(int je=0;je<2;je++) {
	// 			TU[ie+1][je+1]=(*susy)("STOPMIX", ie*10+je);
	// 		}
	// 	}

		
	// 	double vu=sqrt(pow(sin(atan((*susy)("HMIX",2))),2.)/sqrt(2.)/sm("SMINPUTS", 2));
	// 	double vd=vu/(*susy)("HMIX",2);

	// 	for(int je=0;je<2;je++) {
	// 		for(int le=0;le<2;le++) {
	// 			 for(int ae=0;ae<3;ae++) {
	// 				if (ae <3 ){
	// 					Ralj[ae][le][je]=-sm("GAUGE",2)/sqrt(2.)*((*susy)("A0MIX",ae*10+1)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je)+(*susy)("A0MIX",ae*10+2)*(*susy)("UMIX",10+le)*(*susy)("VMIX",20+je))-lambdaNMSSM/sqrt(2.)*(*susy)("A0MIX",ae*10+3)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je);
	// 				}
	// 				Qalj[ae][le][je]=sm("GAUGE",2)/sqrt(2.)*((*susy)("H0MIX",ae*10+1)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je)+(*susy)("H0MIX",ae*10+2)*(*susy)("UMIX",10+le)*(*susy)("VMIX",20+je))-lambdaNMSSM/sqrt(2.)*(*susy)("H0MIX",ae*10+3)*(*susy)("UMIX",20+le)*(*susy)("VMIX",20+je);
	// 				for(int ke=1;ke<=3;ke++) {
	// 					G1[ae][ke][je][le]=(TU[ae][2]*TU[ke][2]-kron(ae,1)*kron(ke,1))*(*susy)("VMIX",10+le)*(*susy)("UMIX",20+je)-(*sus_param).mass_top_muW/sqrt(2.)/sin(atan((*susy)("HMIX",2)))/sm("MASS",24)*TU[ae][3]*TU[ke][2]*(*susy)("VMIX",20+le)*(*susy)("UMIX",20+je);
	// 				}
	// 			}
	// 		}
	// 		for(int ie=0;ie<3;ie++) {
	// 			for(int ke=0;ke<3;ke++) {
	// 				T1[je][ie][ke]=(TU[ie][3]*TU[ke][2]-TU[ie][2]*TU[ke][3])*((lambdaNMSSM/sqrt(2.)*(vd*(*susy)("A0MIX",je*10+3)+s*(*susy)("A0MIX",je*10+1)))-(*susy)("AU", 11)*(*susy)("A0MIX",je*10+2));
	// 			}
	// 		}
	// 	}

	// 	complex_t CQ1H=0.;
	// 	complex_t CQ2H=0.;
	// 	complex_t CQ1c=0.;
	// 	complex_t CQ2c=0.;
	// 	complex_t CAc=0.;

	// 	for(int ae=0;ae<3;ae++) {
	// 		for(int ie=0;ie<3;ie++) {
	// 			for(int ke=0;ke<3;ke++){
	// 				T2[ae][ie][ke]=-(*sus_param).mass_top_muW/2./sm("MASS",24)*(2.*(*sus_param).mass_top_muW*(*susy)("H0MIX",ae*10+2)*(TU[ie][2]*TU[ke][2]+TU[ie][3]*TU[ke][3])	+((lambdaNMSSM/sqrt(2.)*(vd*(*susy)("H0MIX",ae*10+3)+s*(*susy)("H0MIX",ae*10+1)))+(*susy)("AU", 11)*(*susy)("H0MIX",ae*10+2))*(TU[ie][3]*TU[ke][2]+TU[ie][2]*TU[ke][3]))	+src.at({ParameterType::SM, "MASS",  23})->get_val()/2./sqrt(1.-(*sus_param).sw2)*(1.-4./3.*(*sus_param).sw2)*(*susy)("H0MIX",ae*10+2)*(TU[ie][1]*TU[ke][1]+TU[ie][2]*TU[ke][2])+2./3.*sm("MASS",24)*(*sus_param).sw2/(1.-(*sus_param).sw2)*(*susy)("H0MIX",ae*10+2)*TU[ie][3]*TU[ke][3];

	// 				for(int je=0;je<2;je++) {
	// 					for(int le=0;le<2;le++) {
	// 						CQ1c+=G1[ie][ke][je][le]/mH0[ae]/mH0[ae]*( sqrt(2.)*(*susy)("H0MIX",ae*10+1)*(*susy)("H0MIX",ae*10+1)*(*sus_param).Mch[je]/sm("MASS",24)/cos(atan((*susy)("HMIX",2)))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/(*sus_param).Mch[je],2.))
	// 						-2.*sqrt(2.)*(*susy)("H0MIX",ae*10+1)/sm("GAUGE",2)*kron(ie,ke)*(Qalj[ae][le][je]*f40(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))+(*sus_param).Mch[je]/(*sus_param).Mch[le]*Qalj[ae][je][le]*f30(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))		+2.*sqrt(2.)*(*susy)("H0MIX",ae*10+1)*T2[ae][ie][ke]*(*sus_param).Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow((*sus_param).Mch[je]/mstop[ke-1],2.))
	// 						+mH0[ae]*mH0[ae]/(*sus_param).Mch[je]/(*sus_param).Mch[je]*kron(ie,ke)*((*susy)("UMIX",20+je)*(*susy)("VMIX",10+le)*f50(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))
	// 						-(*sus_param).Mch[le]/(*sus_param).Mch[je]*(*susy)("UMIX",20+le)*(*susy)("VMIX",10+je)* f60(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))));

	// 						CQ2c+=G1[ie][ke][je][le]/mA0[ae]/mA0[ae]*(sqrt(2.)*(*susy)("A0MIX",ae*10+1)*(*susy)("A0MIX",ae*10+1)*(*sus_param).Mch[je]/sm("MASS",24)/cos(atan((*susy)("HMIX",2)))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/(*sus_param).Mch[je],2.))
	// 						-2.*sqrt(2.)*(*susy)("A0MIX",ae*10+1)/sm("GAUGE",2)*kron(ie,ke)*(-Ralj[ae][le][je]*f40(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))+(*sus_param).Mch[je]/(*sus_param).Mch[le]*Ralj[ae][je][le]*f30(pow(mstop[ie-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))			-sqrt(2.)*(*susy)("A0MIX",ae*10+1)*T1[ae][ie][ke]*(*sus_param).mass_top_muW*(*sus_param).Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow((*sus_param).Mch[je]/mstop[ke-1],2.))
	// 						+mA0[ae]*mA0[ae]/(*sus_param).Mch[je]/(*sus_param).Mch[je]*kron(ie,ke)*((*susy)("UMIX",20+je)*(*susy)("VMIX",10+le)*f50(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))
	// 						-(*sus_param).Mch[le]/(*sus_param).Mch[je]*(*susy)("UMIX",20+le)*(*susy)("VMIX",10+je)*f60(pow(mstop[ie-1]/(*sus_param).Mch[je],2.),pow((*sus_param).Mch[le]/(*sus_param).Mch[je],2.),pow(mass_nutl/(*sus_param).Mch[le],2.))));
	// 					}		
	// 				}

	// 			}
	// 		}
	// 		CQ1H+=((*susy)("MASS",37)*(*susy)("MASS",37)/sm("MASS",24)/sm("MASS",24)*(*susy)("H0MIX",ae*10+1)*(*susy)("H0MIX",ae*10+1)*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,sm("MASS",24)*sm("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)	+(*sus_param).mass_top_muW*(*sus_param).mass_top_muW*mH0[ae]*mH0[ae]/sm("MASS",24)/sm("MASS",24)/(*susy)("MASS",37)/(*susy)("MASS",37)*f30((*sus_param).mass_top_muW*(*sus_param).mass_top_muW/(*susy)("MASS",37)/(*susy)("MASS",37),(*sus_param).mass_top_muW*(*sus_param).mass_top_muW/sm("MASS",24)/sm("MASS",24)))/mH0[ae]/mH0[ae];
	// 		if (ae < 3) {
	// 			CQ2H+=(((*susy)("MASS",37)*(*susy)("MASS",37)/sm("MASS",24)/sm("MASS",24)*(*susy)("A0MIX",ae*10+1)*(*susy)("A0MIX",ae*10+1)+kron(ae,2)*(*susy)("A0MIX",ae*10+1))*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,sm("MASS",24)*sm("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)	+(*sus_param).mass_top_muW*(*sus_param).mass_top_muW*mA0[ae]*mA0[ae]/sm("MASS",24)/sm("MASS",24)/(*susy)("MASS",37)/(*susy)("MASS",37)*f30((*sus_param).mass_top_muW*(*sus_param).mass_top_muW/(*susy)("MASS",37)/(*susy)("MASS",37),(*sus_param).mass_top_muW*(*sus_param).mass_top_muW/sm("MASS",24)/sm("MASS",24)))/mA0[ae]/mA0[ae];
	// 		}
	// 		for(int je=0;je<2;je++) {
	// 			for(int le=0;le<2;le++) {
	// 				CAc = complex_t(CAc.real(), CAc.imag()+((*susy)("HMIX",2))/sqrt(2.)*G1[ae][ae][je][le]*(v_deltam_s*kron(le,je)*fabs((*sus_param).Mch[je]/sm("MASS",24))*f80(pow(mstop[ae-1]/(*sus_param).Mch[je],2.))-(Ralj[1][je][le]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[le])*f30(pow(mstop[ae-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.))-Ralj[1][le][je]*f40(pow(mstop[ae-1]/(*sus_param).Mch[le],2.),pow((*sus_param).Mch[je]/(*sus_param).Mch[le],2.)))));
	// 				}
	// 			}
	// 	}
	
	// 	CQ1H*=-wilson_p("WPARAM_SI_SM", 3)/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
	// 	CQ2H*=wilson_p("WPARAM_SI_SM", 3)/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
		
	// 	complex_t CAH={0,-lambdaNMSSM*AlambdaNSSM/sm("GAUGE",2)/sm("MASS",24)*(*susy)("HMIX",2)*f30((*susy)("MASS",37)*(*susy)("MASS",37)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW,sm("MASS",24)*sm("MASS",24)/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW)};
	
			
	// 	CQ1c*=wilson_p("WPARAM_SI_SM", 3)/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);
	// 	CQ2c*=-wilson_p("WPARAM_SI_SM", 3)/4.*(*susy)("HMIX",2)*(*susy)("HMIX",2);		
	
	// 	coeff_temp = (CQ2H+CQ2c)*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2/(*sus_param).epsfac;

	// 	complex_t CA=CAH+CAc;

	// 	if((*susy)("MASS",36)>this->get_Q_match()) coeff_temp+=-v_deltam_s/2.*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2*wilson_p("WPARAM_SI_SM", 3)*CA/(*susy)("MASS",36)/(*susy)("MASS",36);

		
	// 	if(((*susy)("MASS",36)>QCDHelper::mass_b_pole())&&((*susy)("MASS",36)<this->get_Q_match() ))
	// 	{	
	// 		double alphas_Ma1  = QCDHelper::alpha_s((*susy)("MASS",36));	
	// 		double mass_b_ma1=QCDHelper::msbar_mass(5, (*susy)("MASS",36));
	// 		coeff_temp+=-v_deltam_s/2.*mass_b_ma1/(*sus_param).sw2*wilson_p("WPARAM_SI_SM", 3)*CA/(*susy)("MASS",36)/(*susy)("MASS",36) *pow(alphas_Ma1/wilson_p("WPARAM_RUN_SM", 1),-4./wilson_p("WPARAM_SI_SM", 5));
	// 	}
		

	// }

    // // return CQ2charg_0/sus_param->epsfac;

}

void CQ2_susy::NLO_calculation() {
	// sus_param->reset_PrimeCQG(this->get_Q_match());
	std::unordered_set<ParamId> sources = {
		{ParameterType::WILSON, "WPARAM_SI_BSM", 7},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 8},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 1},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 2},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 3},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 6},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 9},
		{ParameterType::WILSON, "WPARAM_SI_BSM", 11},
		{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}},
		{ParameterType::WILSON, "WPARAM_MATCH_SM", 6},
		{ParameterType::WILSON, "EW_SCALE", 1},
		{ParameterType::WILSON, "EPSILON_SUSY", 5},
		{ParameterType::WILSON, "WPARAM_SI_SM", 3},
		{ParameterType::WILSON, "WPARAM_SI_SM", 4},
		{ParameterType::SM, "GAUGE", 2},
		{ParameterType::SM, "MASS", 24},
		{ParameterType::SM, "MASS", 3},  // important pour "MASS(3)" (b-quark)
		{ParameterType::BSM, "MASS", 37},
		{ParameterType::BSM, "HMIX", 1},
		{ParameterType::BSM, "HMIX", 2},
	};
	
	// + Paramètres utilisés dans des boucles
	
	// WPARAM_SI_BSM {13, ie}, {14, ae}, {16, be}
	for (int ie = 0; ie < 2; ++ie) {
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}});
	}
	for (int ae = 0; ae < 6; ++ae) {
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}});
	}
	for (int be = 0; be < 3; ++be) {
		sources.insert({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}});
	}
	
	// UMIX, VMIX
	for (int i = 0; i < 6; ++i) {
		sources.insert({ParameterType::BSM, "UMIX", i});
		sources.insert({ParameterType::BSM, "UMIX", i + 10});
		sources.insert({ParameterType::BSM, "VMIX", i});
		sources.insert({ParameterType::BSM, "VMIX", i + 10});
	}
	
	// MATRIX_BSM
	for (int je = 0; je < 2; ++je) {
		for (int ae = 0; ae < 6; ++ae) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {4, je, ae, 2}});
		}
	}
	for (int ie = 0; ie < 2; ++ie) {
		for (int be = 0; be < 3; ++be) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 2}});
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 2}});
		}
	}
	
	for (int ae = 0; ae < 6; ++ae) {
		for (int fe = 0; fe < 3; ++fe) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}});
		}
	}
	for (int me = 0; me < 6; ++me) {
		for (int be = 0; be < 3; ++be) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}});
		}
	}
	
	// MATRIX_BSM 12
	for (int ae = 0; ae < 6; ++ae) {
		for (int ie = 0; ie < 2; ++ie) {
			for (int me = 0; me < 3; ++me) {
				for (int ne = 0; ne < 3; ++ne) {
					sources.insert({ParameterType::WILSON, "MATRIX_BSM", {12, ae, ie, me, ne}});
				}
			}
		}
	}
	
	// MATRIX_BSM 9
	for (int ae = 0; ae < 6; ++ae) {
		for (int ce = 0; ce < 6; ++ce) {
			sources.insert({ParameterType::WILSON, "MATRIX_BSM", {9, ae, ce}});
		}
	}
	
	// WPARAM_MATCH_BSM {2,fe}
	for (int fe = 1; fe <= 3; ++fe) {
		sources.insert({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}});
	}
	
	// RECKM
	for (int i = 0; i < 30; ++i) { //TODO real RECKM
		sources.insert({ParameterType::SM, "RECKM", i});
	}

    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
		double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2,1}})->get_val();

		double lu = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
		double ld = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
		double z = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 1})->get_val();
		double aY = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 11})->get_val();
		double yt = src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", 1})->get_val();
		double Q_match = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
		double mW = src.at({ParameterType::SM, "MASS", 24})->get_val();
		double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double sw2 = src.at({ParameterType::WILSON, "WPARAM_SI_SM", 4})->get_val();
		double mass_b_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5,1}})->get_val();
		double mass_top_muW = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();
		double g2 = src.at({ParameterType::SM, "GAUGE", 2})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();

		double muQ = src.at({ParameterType::BSM, "HMIX", 1})->get_val();
		double xH = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 2})->get_val();
        double xH0 = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 3})->get_val();
		double alpha = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 9})->get_val();
        double beta = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();

		double kappa = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 6})->get_val();
		double epsfac = src.at({ParameterType::WILSON, "EPSILON_SUSY", 5})->get_val();

		complex_t NQ11H=-src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/4./mW/mW*(f141(xt,z)+8.*xt*(f30(xt,z)+xt*(f30(xt*1.0001,z)-f30(xt*0.9999,z))/0.0002)*log(Q_match*Q_match/mass_top_muW/mass_top_muW));
		complex_t BQ11H=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/4./mW/mW*(f111(xt,z)+8.*(f70(xt*1.0001,z)-f70(xt*0.9999,z))/0.0002*log(Q_match*Q_match/mass_top_muW/mass_top_muW));
		complex_t CQ1H_1=(NQ11H+BQ11H)*mass_b_muW/sw2;
		complex_t CQ2H_1=-CQ1H_1;

		
		complex_t BQ11c1=0.;
		complex_t BQ11c2=0.;
		complex_t NQ11c=0.;
		complex_t NQ21c=0.;
		complex_t BQ11f1=0.;
		complex_t BQ11f2=0.;
		complex_t NQ11f=0.;
		complex_t NQ21f=0.;

		double Dp{0}, Dm{0}, temp{0}, temp2{0};
		double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a0p{0}, a1{0},a2p{0};

		for(int ie=0;ie<2;ie++) {
			for(int ae=0;ae<6;ae++){
				for(int je=0;je<2;je++)  {
					for(int be=0;be<3;be++) {
						BQ11c1+=src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {6, ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, be, 1}})->get_val()*(f121(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+4.*(f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))-f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.))));
						BQ11c2+=src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, ae, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*(src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, be, 1}})->get_val()*fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*(f131(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+4.*(f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))-f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.))));

						for(int me=0;me<6;me++){ 
							for(int ne=0;ne<3;ne++) {
								Dp=0.;
								Dm=0.;
								for(int fe=1;fe<=3;fe++) { 	
									Dp+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, fe}})->get_val()+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}})->get_val());
									Dm+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, fe}})->get_val()-src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ae, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, me, fe}})->get_val());
								}
								a0a=-(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val())*(f181(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))+4.*(f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.)*1.0001)-f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.)))*src.at({ParameterType::BSM, "UMIX", ie*10+1})->get_val()*src.at({ParameterType::BSM, "VMIX", je*10+0})->get_val())*kron(ae,me);
								a0b=-((f191(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))+4.*(f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.)*1.0001)-f40(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.)))*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*src.at({ParameterType::BSM, "UMIX", ie*10+0})->get_val())*kron(ae,me);
								a0c=1./mW*(f171(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+4.*(f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+(f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))-f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999,pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)))/0.0002+(f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001)-f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, me}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999))/0.0002)*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.)))*kron(ie,je);
							
								a0Q1=a0a+a0b+Dp*a0c;
								a0Q2=-a0a+a0b+Dm*a0c;
								a0p=4.*src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()/mW/(src.at({ParameterType::SM, "RECKM", be*10+2})->get_val()*src.at({ParameterType::SM, "RECKM", ne*10+1})->get_val()/src.at({ParameterType::SM, "RECKM", 22})->get_val()/src.at({ParameterType::SM, "RECKM", 21})->get_val())/src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*f151(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
								a1=src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/sqrt(2.)/mW*(f161(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))+4.*(f80(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*1.0001)-f80(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*0.9999))/0.0002*log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val(),2.)))*kron(ie,je)*kron(ae,me);
								a2p=src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*(src.at({ParameterType::SM, "RECKM", be*10+2})->get_val()*src.at({ParameterType::SM, "RECKM", ne*10+1})->get_val()/src.at({ParameterType::SM, "RECKM", 22})->get_val()/src.at({ParameterType::SM, "RECKM", 21})->get_val())*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()/2./mW*f151(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
								
								NQ11c+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*(a0Q1+a1*tanb)
								+src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*a0p
								+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*a2p*pow(src.at({ParameterType::SM, "MASS", 3})->get_val()*tanb,2.);	
								NQ21c+=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*(a0Q2+a1*tanb)
								+src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, be, ne}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*a0p
								+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, me, be}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+1})->get_val()*a2p*pow(src.at({ParameterType::SM, "MASS", 3})->get_val()*tanb,2.);
							}
						}
					
					}
					for(int be=0;be<6;be++) {
						for(int ce=0;ce<6;ce++) {
							for(int fe=0;fe<3;fe++) {
								BQ11f1+=-src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()*pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ce, be}})->get_val()*(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val(),2.)))	*(f90(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*src.at({ParameterType::WILSON, "MATRIX_BSM", {6,ie, fe, 2}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {5, je, fe, 2}})->get_val());
								BQ11f2+=-src.at({ParameterType::WILSON, "MATRIX_BSM", {3, je, be, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {4, ie, ae, 2}})->get_val()*pow(mW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.)*src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, ce}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ce, be}})->get_val()*(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val(),2.)))	*(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val())*f100(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, be}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {16, fe}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.))*src.at({ParameterType::WILSON, "MATRIX_BSM", {5, ie, fe, 1}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {6, je, fe, 1}})->get_val());

							}
						}
					}
				}

				for(int me=0;me<3;me++) {
					for(int ne=0;ne<3;ne++) {
						for(int de=0;de<6;de++) {
							for(int ke=0;ke<6;ke++) {
								
								temp2 = src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, me}})->get_val()*src.at({ParameterType::BSM, "UMIX", ie*10+2})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, ke}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ke, de}})->get_val()*(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(),2.)))*tanb*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/sqrt(2.)*f30(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.));
								NQ11f+=temp2;
								NQ21f+=temp2;
								for(int ce=0;ce<6;ce++) {
									Dp=0.;
									Dm=0.;
									for(int fe=0;fe<3;fe++) 
									{		
										Dp+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, de, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val()+src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ce, fe}})->get_val()); 
										Dm+=src.at({ParameterType::WILSON, "WPARAM_MATCH_BSM", {2, fe}})->get_val()/sqrt(2.)/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*muQ*(src.at({ParameterType::WILSON, "MATRIX_BSM", {2, de, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, fe}})->get_val()-src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, fe}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {2, ce, fe}})->get_val()); 
									}
									temp=src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, ce, me}})->get_val()*src.at({ParameterType::BSM, "UMIX", ie*10+2})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, ke}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ke, de}})->get_val()*
									(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(),2.)))*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ce}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val(),2.));	
									NQ11f+=Dp*temp;
									NQ21f+=Dm*temp;
								}
								for(int je=0;je<2;je++) {
									temp=-src.at({ParameterType::WILSON, "MATRIX_BSM", {12,ae, ie, me, ne}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {1, de, me}})->get_val()*src.at({ParameterType::BSM, "UMIX", je*10+2})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9,ae, ke}})->get_val()*src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val()*src.at({ParameterType::WILSON, "MATRIX_BSM", {9, ke, de}})->get_val()*
									(1.+log(pow(mass_top_muW/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ke}})->get_val(),2.)))*mW*(fabs(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val())*f60(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX", ie*10+2})->get_val()*src.at({ParameterType::BSM, "VMIX", je*10+1})->get_val()+
									f50(pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, ie}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, ae}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.),pow(src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {14, de}})->get_val()/src.at({ParameterType::WILSON, "WPARAM_SI_BSM", {13, je}})->get_val(),2.))*src.at({ParameterType::BSM, "UMIX", je*10+2})->get_val()*src.at({ParameterType::BSM, "VMIX", ie*10+1})->get_val()); 
									NQ11f+=temp;
									NQ21f+=-temp;
								}

							}
						}
					}
						
				}
			}
		}
		
		
		complex_t BQ11c=(BQ11c1+BQ11c2)*kappa*mW*mW/2./g2/g2/sw2;
		complex_t BQ21c=-(BQ11c1-BQ11c2)*kappa*mW*mW/2./g2/g2/sw2;


		NQ11c*=src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
		NQ21c*=-src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;
		
		complex_t CQ2charg_1=NQ21c+BQ21c;
			
		complex_t BQ11f=(BQ11f1+BQ11f2)*2./3.*kappa/g2/g2/sw2;
		complex_t BQ21f=-(BQ11f1-BQ11f2)*2./3.*kappa/g2/g2/sw2;
		
		
		NQ11f*=-4./3.*src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/mW/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;

		NQ21f*=4./3.*src.at({ParameterType::WILSON, "WPARAM_SI_SM", 3})->get_val()*(tanb)*tanb/mW/mW/(mH*mH-mW*mW)*aY*mass_b_muW/sw2;

		complex_t CQ2four_1=NQ21f+BQ21f;

		complex_t coeff_temp = (CQ2H_1+CQ2charg_1)/epsfac+CQ2four_1;
		// this->set_WilsonCoeffMatching("NLO", coeff_temp);
		std::cout << coeff_temp << std::endl;

        dep_param->set_expected(coeff_temp);
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3233, 1, 1)}, sources, func);


	// // sus_param->reset_PrimeCQG(this->get_Q_match());

	// complex_t NQ11H=-wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./sm("MASS",24)/sm("MASS",24)*(f141((*sus_param).xt,(*sus_param).z)+8.*(*sus_param).xt*(f30((*sus_param).xt,(*sus_param).z)+(*sus_param).xt*(f30((*sus_param).xt*1.0001,(*sus_param).z)-f30((*sus_param).xt*0.9999,(*sus_param).z))/0.0002)*log(this->get_Q_match()*this->get_Q_match()/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW));
	// complex_t BQ11H=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/4./sm("MASS",24)/sm("MASS",24)*(f111((*sus_param).xt,(*sus_param).z)+8.*(f70((*sus_param).xt*1.0001,(*sus_param).z)-f70((*sus_param).xt*0.9999,(*sus_param).z))/0.0002*log(this->get_Q_match()*this->get_Q_match()/(*sus_param).mass_top_muW/(*sus_param).mass_top_muW));
	// complex_t CQ1H_1=(NQ11H+BQ11H)*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;
	// complex_t CQ2H_1=-CQ1H_1;

	
	// complex_t BQ11c1=0.;
	// complex_t BQ11c2=0.;
	// complex_t NQ11c=0.;
	// complex_t NQ21c=0.;
	// complex_t BQ11f1=0.;
	// complex_t BQ11f2=0.;
	// complex_t NQ11f=0.;
	// complex_t NQ21f=0.;

	// double Dp{0}, Dm{0}, temp{0}, temp2{0};
	// double a0a{0}, a0b{0}, a0c{0}, a0Q1{0}, a0Q2{0}, a0p{0}, a1{0},a2p{0};

	// for(int ie=0;ie<2;ie++) {
	// 	for(int ae=0;ae<6;ae++){
	// 		for(int je=0;je<2;je++)  {
	// 			for(int be=0;be<3;be++) {
	// 				BQ11c1+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NR[ie][be][1]*(*sus_param).X_NL[je][be][1]*(f121(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))+4.*(f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))-f50(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.))));
	// 				BQ11c2+=(*sus_param).X_UL[je][ae][1]*(*sus_param).X_UR[ie][ae][2]/(*sus_param).Mch[ie]/(*sus_param).Mch[ie]*((*sus_param).X_NL[ie][be][1]*(*sus_param).X_NR[je][be][1]*fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*(f131(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))+4.*(f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.))-f60(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).Msn[be]/(*sus_param).Mch[ie],2.)))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.))));

	// 				for(int me=0;me<6;me++){ 
	// 					for(int ne=0;ne<3;ne++) {
	// 						Dp=0.;
	// 						Dm=0.;
	// 						for(int fe=1;fe<=3;fe++) { 	
	// 							Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]+(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
	// 							Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[ae][fe]*(*sus_param).Gamma_UL[me][fe]-(*sus_param).Gamma_UL[ae][fe]*(*sus_param).Gamma_UR[me][fe]);
	// 						}
	// 						a0a=-(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*(f181(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))+4.*(f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*1.0001)-f30(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*(*susy)("UMIX", ie*10+1)*(*susy)("VMIX", je*10+0))*kron(ae,me);
	// 						a0b=-((f191(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.))+4.*(f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*1.0001)-f40(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*(*susy)("UMIX", je*10+1)*(*susy)("UMIX", ie*10+0))*kron(ae,me);
	// 						a0c=1./sm("MASS",24)*(f171(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))+4.*(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))+(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001,pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.))-f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999,pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)))/0.0002+(f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)*1.0001)-f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[me]/(*sus_param).Mch[ie],2.)*0.9999))/0.0002)*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*kron(ie,je);
						
	// 						a0Q1=a0a+a0b+Dp*a0c;
	// 						a0Q2=-a0a+a0b+Dm*a0c;
	// 						a0p=4.*(*sus_param).G_aimn[ae][ie][be][ne]/sm("MASS",24)/(src.at({ParameterType::SM, "RECKM", be*10+2})->get_val()*src.at({ParameterType::SM, "RECKM", ne*10+1})->get_val()/src.at({ParameterType::SM, "RECKM", 22})->get_val()/src.at({ParameterType::SM, "RECKM", 21})->get_val())/(*susy)("UMIX", je*10+1)*f151(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
	// 						a1=(*sus_param).Mch[ie]/sqrt(2.)/sm("MASS",24)*(f161(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))+4.*(f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*1.0001)-f80(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.)*0.9999))/0.0002*log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ae],2.)))*kron(ie,je)*kron(ae,me);
	// 						a2p=(*sus_param).Gamma_UL[me][be]*(src.at({ParameterType::SM, "RECKM", be*10+2})->get_val()*src.at({ParameterType::SM, "RECKM", ne*10+1})->get_val()/src.at({ParameterType::SM, "RECKM", 22})->get_val()/src.at({ParameterType::SM, "RECKM", 21})->get_val())*(*susy)("UMIX", je*10+1)/2./sm("MASS",24)*f151(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.))*kron(ie,je)*kron(ae,me)*kron(be,ne);
							
	// 						NQ11c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q1+a1*(*susy)("HMIX",2))
	// 						+(*sus_param).G_aimn[ae][ie][be][ne]*(*susy)("UMIX", je*10+1)*a0p
	// 						+(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*a2p*pow(sm("MASS",3)*(*susy)("HMIX",2),2.);	
	// 						NQ21c+=(*sus_param).G_aimn[ae][ie][be][ne]*(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*(a0Q2+a1*(*susy)("HMIX",2))
	// 						+(*sus_param).G_aimn[ae][ie][be][ne]*(*susy)("UMIX", je*10+1)*a0p
	// 						+(*sus_param).Gamma_UL[me][be]*(*susy)("UMIX", je*10+1)*a2p*pow(sm("MASS",3)*(*susy)("HMIX",2),2.);
	// 					}
	// 				}
				
	// 			}
	// 			for(int be=0;be<6;be++) {
	// 				for(int ce=0;ce<6;ce++) {
	// 					for(int fe=0;fe<3;fe++) {
	// 						BQ11f1+=-(*sus_param).X_UL[je][be][1]*(*sus_param).X_UR[ie][ae][2]*pow(sm("MASS",24)/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][ce]*(*sus_param).MsqU[ce]/(*sus_param).Mch[ie]*(*sus_param).P_U[ce][be]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ce],2.)))	*(f90(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[fe]/(*sus_param).Mch[ie],2.))*(*sus_param).X_NR[ie][fe][2]*(*sus_param).X_NL[je][fe][2]);
	// 						BQ11f2+=-(*sus_param).X_UL[je][be][1]*(*sus_param).X_UR[ie][ae][2]*pow(sm("MASS",24)/(*sus_param).Mch[ie],2.)*(*sus_param).P_U[ae][ce]*(*sus_param).MsqU[ce]/(*sus_param).Mch[ie]*(*sus_param).P_U[ce][be]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ce],2.)))	*(fabs((*sus_param).Mch[je]/(*sus_param).Mch[ie])*f100(pow((*sus_param).Mch[je]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[be]/(*sus_param).Mch[ie],2.),pow((*sus_param).Msn[fe]/(*sus_param).Mch[ie],2.))*(*sus_param).X_NL[ie][fe][1]*(*sus_param).X_NR[je][fe][1]);

	// 					}
	// 				}
	// 			}
	// 		}

	// 		for(int me=0;me<3;me++) {
	// 			for(int ne=0;ne<3;ne++) {
	// 				for(int de=0;de<6;de++) {
	// 					for(int ke=0;ke<6;ke++) {
							
	// 						temp2 = (*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[de][me]*(*susy)("UMIX",ie*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[ie]*(*sus_param).P_U[ke][de]*(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*(*susy)("HMIX",2)*(*sus_param).Mch[ie]/sqrt(2.)*f30(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2.));
	// 						NQ11f+=temp2;
	// 						NQ21f+=temp2;
	// 						for(int ce=0;ce<6;ce++) {
	// 							Dp=0.;
	// 							Dm=0.;
	// 							for(int fe=0;fe<3;fe++) 
	// 							{		
	// 								Dp+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[de][fe]*(*sus_param).Gamma_UL[ce][fe]+(*sus_param).Gamma_UL[de][fe]*(*sus_param).Gamma_UR[ce][fe]); 
	// 								Dm+=(*sus_param).MU[fe]/sqrt(2.)/(*sus_param).Mch[ie]*(*susy)("HMIX",1)*((*sus_param).Gamma_UR[de][fe]*(*sus_param).Gamma_UL[ce][fe]-(*sus_param).Gamma_UL[de][fe]*(*sus_param).Gamma_UR[ce][fe]); 
	// 							}
	// 							temp=(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[ce][me]*(*susy)("UMIX",ie*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[ie]*(*sus_param).P_U[ke][de]*
	// 							(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*f60(pow((*sus_param).MsqU[ae]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[ie],2.),pow((*sus_param).MsqU[ce]/(*sus_param).Mch[ie],2.));	
	// 							NQ11f+=Dp*temp;
	// 							NQ21f+=Dm*temp;
	// 						}
	// 						for(int je=0;je<2;je++) {
	// 							temp=-(*sus_param).G_aimn[ae][ie][me][ne]*(*sus_param).Gamma_UL[de][me]*(*susy)("UMIX",je*10+2)*(*sus_param).P_U[ae][ke]*(*sus_param).MsqU[ke]/(*sus_param).Mch[je]*(*sus_param).P_U[ke][de]*
	// 							(1.+log(pow((*sus_param).mass_top_muW/(*sus_param).MsqU[ke],2.)))*sm("MASS",24)*(fabs((*sus_param).Mch[ie]/(*sus_param).Mch[je])*f60(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[je],2.))*(*susy)("UMIX",ie*10+2)*(*susy)("VMIX",je*10+1)+
	// 							f50(pow((*sus_param).Mch[ie]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[ae]/(*sus_param).Mch[je],2.),pow((*sus_param).MsqU[de]/(*sus_param).Mch[je],2.))*(*susy)("UMIX",je*10+2)*(*susy)("VMIX",ie*10+1)); 
	// 							NQ11f+=temp;
	// 							NQ21f+=-temp;
	// 						}

	// 					}
	// 				}
	// 			}
					
	// 		}
	// 	}
	// }
	
	
	// complex_t BQ11c=(BQ11c1+BQ11c2)*(*sus_param).kappa*sm("MASS",24)*sm("MASS",24)/2./sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;
	// complex_t BQ21c=-(BQ11c1-BQ11c2)*(*sus_param).kappa*sm("MASS",24)*sm("MASS",24)/2./sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;


	// NQ11c*=wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;
	// NQ21c*=-wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;
	
	// complex_t CQ2charg_1=NQ21c+BQ21c;
		
	// complex_t BQ11f=(BQ11f1+BQ11f2)*2./3.*(*sus_param).kappa/sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;
	// complex_t BQ21f=-(BQ11f1-BQ11f2)*2./3.*(*sus_param).kappa/sm("GAUGE",2)/sm("GAUGE",2)/(*sus_param).sw2;
	
	
	// NQ11f*=-4./3.*wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;

	// NQ21f*=4./3.*wilson_p("WPARAM_SI_SM", 3)*((*susy)("HMIX",2))*(*susy)("HMIX",2)/sm("MASS",24)/sm("MASS",24)/((*susy)("MASS",37)*(*susy)("MASS",37)-sm("MASS",24)*sm("MASS",24))*(*sus_param).aY*wilson_p("WPARAM_MATCH_SM", {5,1})/(*sus_param).sw2;

	// complex_t CQ2four_1=NQ21f+BQ21f;

    // complex_t coeff_temp = (CQ2H_1+CQ2charg_1)/sus_param->epsfac+CQ2four_1;
    // this->set_WilsonCoeffMatching("NLO", coeff_temp);
	// std::cout << coeff_temp << std::endl;
	// return coeff_temp;
}

void BScalarCoefficientGroup_susy::set_base_1_LO() {

	
	LOG_INFO("In BScalarCoefficientGroup::set_base_1_LO");

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {"B_SCALAR_MATCH", "WPARAM_RUN_SM", "WPARAM_SI_SM"}},
    };

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double g2 = src.at("GAUGE")->retrieve(2)->get_val();
		double tanb = src.at("HMIX")->retrieve(2)->get_val();
		double mW = src.at("MASS")->retrieve(24)->get_val();
		
		double mH = src.at("MASS")->retrieve(37)->get_val();
		double mA = src.at("MASS")->retrieve(36)->get_val();

		double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        double beta_0 = src.at("WPARAM_SI_SM")->retrieve(5)->get_val();

		double mass_top_muW = src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val();
		double sw2 = src.at("WPARAM_SI_SM")->retrieve(4)->get_val();

		auto ensure_coef = [src] (const LhaID& id) -> complex_t {
            return src.at("B_SCALAR_MATCH")->contains(id) ? src.at("B_SCALAR_MATCH")->retrieve(id)->get_val() : complex_t(0);
		};


		complex_t coeff_temp = ensure_coef(WCoefMapper::flha_full(WCoef::CQ1, QCDOrder::LO, ContributionType::BSM)) * pow(eta,-4./beta_0);

		ParamId pid {ParameterType::WILSON, "B_HADRONIC", WCoefMapper::flha_full(WCoef::CQ1, QCDOrder::LO, ContributionType::BSM)}; //PID TO CHANGE
        dep_block->store_or_assign(1, std::make_shared<Parameter>(pid, coeff_temp, 0., 0.));


		complex_t coeff_temp2= ensure_coef(WCoefMapper::flha_full(WCoef::CQ2, QCDOrder::LO, ContributionType::BSM)) * pow(eta,-4./beta_0);
		
		
		if(src.at("MASS")->retrieve(46)->get_val()!=0.||src.at("MASS")->retrieve(45)->get_val()!=0.) {
			if(mA < (*Parameters::GetInstance())("QCD", LhaID(5, 2))) {	
				double lambdaNMSSM = 1;
				double lambdaSNMSSM = 1;
				double AlambdaNSSM = 1;
				double kappaNMSSM = 1;
				double m_Bs = 1;
				double mass_nutl = 1;

				double mH0[4],mA0[3],mstop[3];
			
				mstop[0]=src.at("MASS")->retrieve(2000013)->get_val(); //mass upr, is that right ?
				mstop[1]=src.at("MASS")->retrieve(1000006)->get_val();
				mstop[2]=src.at("MASS")->retrieve(2000006)->get_val();

				complex_t CAH={0,-lambdaNMSSM*AlambdaNSSM/g2/mW*tanb*f30(mH*mH/mass_top_muW/mass_top_muW,mW*mW/mass_top_muW/mass_top_muW)};
				complex_t CAc{};
				double s=lambdaSNMSSM/lambdaNMSSM;
				double v=sqrt(1./sqrt(2.)/src.at("SMINPUTS")->retrieve(2)->get_val());
				double v_deltam_s=v/s*(sqrt(2.)*AlambdaNSSM-2.*kappaNMSSM*s)/(sqrt(2.)*AlambdaNSSM+kappaNMSSM*s);

				double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
				double T2[4][4][4];
				std::array<std::array<double,4>,4> TU;
				double vu=sqrt(pow(sin(atan(tanb)),2.)/sqrt(2.)/src.at("SMINPUTS")->retrieve(2)->get_val());
				double vd=vu/tanb;

				TU[1][1]=1.;
				for(int ie=0;ie<2;ie++){
					for(int je=0;je<2;je++) {
						TU[ie+1][je+1]=src.at("STOPMIX")->retrieve(ie*10+je)->get_val();
					}
				}

				for(int je=0;je<2;je++) {
					for(int le=0;le<2;le++) {
						for(int ae=0;ae<3;ae++) {
							if (ae <3 ){
								Ralj[ae][le][je]=-g2/sqrt(2.)*(src.at("A0MIX")->retrieve(ae*10+1)->get_val()*src.at("UMIX")->retrieve(20+le)->get_val()*src.at("VMIX")->retrieve(20+je)->get_val()+src.at("A0MIX")->retrieve(ae*10+2)->get_val()*src.at("UMIX")->retrieve(10+le)->get_val()*src.at("VMIX")->retrieve(20+je)->get_val())-lambdaNMSSM/sqrt(2.)*src.at("A0MIX")->retrieve(ae*10+3)->get_val()*src.at("UMIX")->retrieve(20+le)->get_val()*src.at("VMIX")->retrieve(20+je)->get_val();
							}
							Qalj[ae][le][je]=g2/sqrt(2.)*(src.at("H0MIX")->retrieve(ae*10+1)->get_val()*src.at("UMIX")->retrieve(20+le)->get_val()*src.at("VMIX")->retrieve(20+je)->get_val()+src.at("H0MIX")->retrieve(ae*10+2)->get_val()*src.at("UMIX")->retrieve(10+le)->get_val()*src.at("VMIX")->retrieve(20+je)->get_val())-lambdaNMSSM/sqrt(2.)*src.at("H0MIX")->retrieve(ae*10+3)->get_val()*src.at("UMIX")->retrieve(20+le)->get_val()*src.at("VMIX")->retrieve(20+je)->get_val();
							for(int ke=1;ke<=3;ke++) {
								G1[ae][ke][je][le]=(TU[ae][2]*TU[ke][2]-kron(ae,1)*kron(ke,1))*src.at("VMIX")->retrieve(10+le)->get_val()*src.at("UMIX")->retrieve(20+je)->get_val()-mass_top_muW/sqrt(2.)/sin(atan(tanb))/mW*TU[ae][3]*TU[ke][2]*src.at("VMIX")->retrieve(20+le)->get_val()*src.at("UMIX")->retrieve(20+je)->get_val();
							}
						}
					}
				}
				for(int ae=0;ae<3;ae++) {
					for(int je=0;je<2;je++) {
						for(int le=0;le<2;le++) {
							CAc = complex_t(CAc.real(), CAc.imag()+(tanb)/sqrt(2.)*G1[ae][ae][je][le]*(v_deltam_s*kron(le,je)*fabs(src.at("WPARAM_SI_BSM")->retrieve(je)->get_val()/mW)*f80(pow(mstop[ae-1]/src.at("WPARAM_SI_BSM")->retrieve(je)->get_val(),2.))-(Ralj[1][je][le]*fabs(src.at("WPARAM_SI_BSM")->retrieve(je)->get_val()/src.at("WPARAM_SI_BSM")->retrieve(le)->get_val())*f30(pow(mstop[ae-1]/src.at("WPARAM_SI_BSM")->retrieve(le)->get_val(),2.),pow(src.at("WPARAM_SI_BSM")->retrieve(je)->get_val()/src.at("WPARAM_SI_BSM")->retrieve(le)->get_val(),2.))-Ralj[1][le][je]*f40(pow(mstop[ae-1]/src.at("WPARAM_SI_BSM")->retrieve(le)->get_val(),2.),pow(src.at("WPARAM_SI_BSM")->retrieve(je)->get_val()/src.at("WPARAM_SI_BSM")->retrieve(le)->get_val(),2.)))));
						}
					}
				}
				complex_t CA=CAH+CAc;
				double width_A0=1.e-6;
				coeff_temp2+=complex_t{v_deltam_s/2.*(*Parameters::GetInstance())("QCD", LhaID(5, 1))/sw2*src.at("WPARAM_SI_SM")->retrieve(3)->get_val()*CA/(m_Bs*m_Bs-mA*mA,mA*width_A0)};
			}
		}
		ParamId pid2 {ParameterType::WILSON, "B_HADRONIC", WCoefMapper::flha_full(WCoef::CQ2, QCDOrder::LO, ContributionType::BSM)}; //TODO CHECK
        dep_block->store_or_assign(2, std::make_shared<Parameter>(pid2, coeff_temp2, 0., 0.));

    };

    WilsonParamComposer().compose_block("B_SCALAR_HADRONIC", src, func);

}

void C_Blnu_P_SUSY::LO_calculation() {

	std::unordered_set<ParamId> sources {
        {"WPARAM_SI_SM", 4},
        {"WPARAM_SI_BSM", 4},
        {"WPARAM_SI_BSM", 6},
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_SI_BSM", 8},
        {"WPARAM_SI_BSM", 10},
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", {5,1}},
        {ParameterType::SM, "MASS", 24}
    };


    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
        double m_b = (*Parameters::GetInstance())("QCD", LhaID(5, 1));
        double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
		double epsilon0 = src.at({ParameterType::WILSON, "EPSILON_SUSY", {0,1}})->get_val();
        dep_param->set_expected(m_b * m_tau * std::pow(tanb / mH, 2) / (1 + epsilon0 * tanb));
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3230, 0, 1)}, sources, func);

	// double m_b = QCDHelper::mass_b_msbar();
    // double m_tau = sm("MASS", 15);
	// double tanb = (*susy)("HMIX", 2);
    // // return this->double_to_complex_save("LO", m_b * m_tau * std::pow(tanb / sus_param->m_H, 2) / (1 + sus_param->epsilon0 * tanb));
}

void C_S1_SUSY::LO_calculation() {

	std::unordered_set<ParamId> sources {
        {"WPARAM_SI_SM", 4},
        {"WPARAM_SI_BSM", 4},
        {"WPARAM_SI_BSM", 6},
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_SI_BSM", 8},
        {"WPARAM_SI_BSM", 10},
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", {5,1}},
        {ParameterType::SM, "MASS", 24}
    };


    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
        double m_b = (*Parameters::GetInstance())("QCD", LhaID(5, 1));
        double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
		double epsilon0 = src.at({ParameterType::WILSON, "EPSILON_SUSY", {0,1}})->get_val();
        dep_param->set_expected(m_b * m_tau * std::pow(tanb / mH, 2) / (1 + epsilon0 * tanb));
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3230, 0, 1)}, sources, func);

    // double m_b = QCDHelper::mass_b_msbar();
    // double m_tau = sm("MASS", 15);
	// double tanb = (*susy)("HMIX", 2);
    // // return this->double_to_complex_save("LO", m_b * m_tau * std::pow(tanb / sus_param->m_H, 2) / (1 + sus_param->epsilon0 * tanb));
}

void C_S2_SUSY::LO_calculation() {

	std::unordered_set<ParamId> sources {
        {"WPARAM_SI_SM", 4},
        {"WPARAM_SI_BSM", 4},
        {"WPARAM_SI_BSM", 6},
        {"WPARAM_SI_BSM", 7},
        {"WPARAM_SI_BSM", 8},
        {"WPARAM_SI_BSM", 10},
        {"WPARAM_MATCH_SM", {2,1}},
        {"WPARAM_MATCH_SM", {5,1}},
        {ParameterType::SM, "MASS", 24}
    };


    auto func = [] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double mH = src.at({ParameterType::BSM, "MASS", 37})->get_val();
		double tanb = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
        double m_c =  src.at({ParameterType::SM, "MASS", 4})->get_val();
        double m_tau = src.at({ParameterType::SM, "MASS", 15})->get_val();
		double epsilon0 = src.at({ParameterType::WILSON, "EPSILON_SUSY", {0,1}})->get_val();
        dep_param->set_expected(m_c * m_tau * std::pow(tanb / mH, 2) / (1 + epsilon0 * tanb));
    };

    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, LhaID(3051313, 3230, 0, 1)}, sources, func);

    // double m_c = sm("MASS", 4);
    // double m_tau = sm("MASS", 15);
	// double tanb = (*susy)("HMIX", 2);
    // // return this->double_to_complex_save("LO", m_c * m_tau * std::pow(tanb / sus_param->m_H, 2) / (1 + sus_param->epsilon0 * tanb));
}

void WilsonCoefficient_susy::init(QCDOrder order) {
	if (!is_owned) {
        susy_parameters::init();
    }

    WilsonCoefficient::init(order);
}
