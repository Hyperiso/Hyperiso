#include "WilsonGroup.h"

void BCoefficientGroup::set_base_1_LO() {
    LOG_INFO("In BCoefficientGroup::set_base_1_LO");

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {"B_MATCH", "WPARAM_RUN_SM", "WPARAM_MATCH_SM", "SCALES", "U_MATRIX"}},
        {ParameterType::SM, {"SMINPUTS", "MASS"}}
    };

    auto func = [this] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        auto ensure_coef = [src] (const LhaID& id) -> complex_t {
            return src.at("B_MATCH")->contains(id) ? src.at("B_MATCH")->retrieve(id)->get_val() : complex_t(0);
        };

        auto U0 = [src] (int k, int l) -> double {
            return src.at("U_MATRIX")->retrieve(LhaID(0, k, l))->get_val();
        };

        std::array<LhaID, 10> ids = {
            LhaID(3040405, 6161, 0, 0),
            LhaID(3040405, 4141, 0, 0),
            LhaID(3050707, 4133, 0, 0),
            LhaID(3050707, 6153, 0, 0),
            LhaID(3050707, 4536, 0, 0),
            LhaID(3050707, 6556, 0, 0),
            LhaID(305, 4422, 0, 0),
            LhaID(305, 6421, 0, 0),
            LhaID(3051313, 4133, 0, 0),
            LhaID(3051313, 4137, 0, 0)
        };

        std::array<complex_t, 10> Ci_match = {};
        for (size_t k = 0; k < 10; k++) {
            Ci_match[k] = ensure_coef(ids[k]);
        }

        Ci_match[6] += -1. / 3. * Ci_match[2] - 4. / 9. * Ci_match[3] - 20. / 3. * Ci_match[4] - 80./9. * Ci_match[5]; 
        Ci_match[7] += Ci_match[2] - 1 / 6. * Ci_match[3] + 20. * Ci_match[4] - 10. / 3. * Ci_match[5];

        std::array<complex_t, 10> Ci_run {};

        // C1 - C9
        for (size_t k = 0; k < 9; k++) {
            for (size_t l = 0; l < 9; l++) {
                Ci_run[k] += U0(k, l) * Ci_match[l];
            }
            LOG_INFO("C_run_", k + 1, "=", Ci_run[k]);
        }

        double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
        Ci_run[9] *= fact;

        // C10
        // Ask Nazila : Pourquoi les N et NNLO dans le running de C10 à LO et pas séparé en plusieurs ordres ?
        double alpha_ew = 1 / src.at("SMINPUTS")->retrieve(1)->get_val();
        double m_h = src.at("MASS")->retrieve(25)->get_val();
        double m_t_muW = src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val();
        double sw2OS = src.at("SMINPUTS")->retrieve(LhaID(7, 2))->get_val();
        double mu_h = src.at("SCALES")->retrieve(2)->get_val();
        double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();

        complex_t C1_NLO = ensure_coef(LhaID(3040405, 6161, 1, 0));
        complex_t C4_NLO = ensure_coef(LhaID(3050707, 6153, 1, 0));

        complex_t C10_02=0.;
        for(size_t ie = 0; ie < 8; ie++) {
            C10_02 += BWilsonRunningParameters::b[ie] * pow(eta, BWilsonRunningParameters::a[ie]) * Ci_match[1];
        } 
        
        complex_t C10_12 = -0.11060 * std::log(eta) / eta * Ci_match[1]  + (1 / eta - 1) * (0.26087 * Ci_match[8] + 1.15942 * Ci_match[9]);
	    for(int ie=0;ie<=7;ie++) {
            C10_12 += pow(eta, BWilsonRunningParameters::a[ie] + 1.) * (
                (BWilsonRunningParameters::d_2a[ie] / eta + BWilsonRunningParameters::d_2b[ie]) * Ci_match[1] 
                + BWilsonRunningParameters::d_1[ie] * C1_NLO
                + BWilsonRunningParameters::d_4[ie] * C4_NLO);
        }

        double Delta_alpha = 0.06; 
        double Delta_rhosw2 = -0.03; 
        double Delta_rem = 0.01;
        double Deltar = Delta_alpha + Delta_rhosw2 + Delta_rem; 
	
	    double Gmu1_Gmu0 = 4 * PI / alpha_ew * Deltar;
        double L = std::log(mu_h * mu_h);
        complex_t C1022 = (46.9287715663914 - 3.102350691200236 * L + 0.0992974073578769 * L * L + 0.175877 * (m_t_muW - 163.5) + 0.0173725 * (m_h - 125.9)) / sw2OS;
        C1022 += -Ci_match[9] * Gmu1_Gmu0;

        complex_t C9_NLO = ensure_coef(LhaID(3051313, 4133, 1, 0));
        complex_t C10_NLO = ensure_coef(LhaID(3051313, 4137, 1, 0));
        
        complex_t C10_22 = (0.27924 * C1_NLO + 0.33157 * C4_NLO + 2.35917 * Ci_match[8] + 3.29679 * Ci_match[9]) * log(eta) + (1 - eta) * (0.26087 * C9_NLO + 1.15942 * C10_NLO) + C1022;
        
        complex_t C1_NNLO = ensure_coef(LhaID(3040405, 6161, 2, 0));
        complex_t C2_NNLO = ensure_coef(LhaID(3040405, 4141, 2, 0));
        complex_t C3_NNLO = ensure_coef(LhaID(3050707, 4133, 2, 0));
        complex_t C4_NNLO = ensure_coef(LhaID(3050707, 6153, 2, 0));
        complex_t C5_NNLO = ensure_coef(LhaID(3050707, 4536, 2, 0));
        complex_t C6_NNLO = ensure_coef(LhaID(3050707, 6556, 2, 0));

        for(int ie = 0; ie <= 7; ie++) {
            C10_22 += pow(eta, BWilsonRunningParameters::a[ie] + 2) * (
                      (BWilsonRunningParameters::e_1a[ie] / eta + BWilsonRunningParameters::e_1b[ie]) * C1_NLO
                    + (BWilsonRunningParameters::e_4a[ie] / eta + BWilsonRunningParameters::e_4b[ie]) * C4_NLO
                    + BWilsonRunningParameters::e_1[ie] * C1_NNLO
                    + BWilsonRunningParameters::e_2[ie] * C2_NNLO
                    + BWilsonRunningParameters::e_3[ie] * C3_NNLO
                    + BWilsonRunningParameters::e_4[ie] * C4_NNLO
                    + BWilsonRunningParameters::e_5[ie] * C5_NNLO
                    + BWilsonRunningParameters::e_6[ie] * C6_NNLO);
        }

        double alpha_s_mu_h = src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
        Ci_run[9] = Ci_match[9] + alpha_ew / alpha_s_mu_h * (fact * C10_02 + C10_12) + alpha_ew / (4 * PI) * C10_22;
        
        // Store
        for (size_t k = 0; k < 10; k++) {
            dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "B_HADRONIC", ids[k]}, Ci_run[k], 0., 0.));
        }
    };

    WilsonParamComposer().compose_block("B_HADRONIC", src, func);
    this->base["LO"] = 1;
}

void BCoefficientGroup::set_base_2_LO() {
    // std::vector<std::string> coeffs {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    // std::vector<complex_t> coeffs_t {};
    // std::vector<complex_t> coeffs_b(10);
    // coeffs_t.push_back(this->at("C1")->get_CoefficientMatchingValue("LO")/2.);

    // coeffs_t.push_back(-this->at("C1")->get_CoefficientMatchingValue("LO")/6.+this->at("C2")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back((54.*this->at("C3")->get_CoefficientMatchingValue("LO")-9.*this->at("C4")->get_CoefficientMatchingValue("LO")+864.*this->at("C5")->get_CoefficientMatchingValue("LO")-144.*this->at("C6")->get_CoefficientMatchingValue("LO"))/54.); 
	// coeffs_t.push_back(this->at("C4")->get_CoefficientMatchingValue("LO")/2.+8.*this->at("C6")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back((54.*this->at("C3")->get_CoefficientMatchingValue("LO")-9.*this->at("C4")->get_CoefficientMatchingValue("LO")+216.*this->at("C5")->get_CoefficientMatchingValue("LO")-36.*this->at("C6")->get_CoefficientMatchingValue("LO"))/54.);
 	// coeffs_t.push_back(this->at("C4")->get_CoefficientMatchingValue("LO")/2.+2.*this->at("C6")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(this->at("C7")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(this->at("C8")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(this->at("C9")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(this->at("C10")->get_CoefficientMatchingValue("LO"));

    // complex_t C0t7= coeffs_t[6]-1./3.*coeffs_t[4]-coeffs_t[5]; 

	// complex_t C0t8= coeffs_t[7]+coeffs_t[4]; 
    // for (int i=0; i<8; i++) {
    //     for (int j=0; j<8;j++) {
    //         if (j<6)
	// 	    {
	// 		    coeffs_b[i] += rh.get_matrix().V0[i][j]*coeffs_t[j];
	// 	    }
	// 	    if (j==6)
	// 	    {
	// 	    	coeffs_b[i] += rh.get_matrix().V0[i][j]*C0t7;
	// 	    }
	// 	    if (j==7)
	// 	    {
	// 		    coeffs_b[i] += rh.get_matrix().V0[i][j]*C0t8;
	// 	    }
    //     }
    // }
    // for (int j=0; j<8; j++) {
    //     coeffs_b[8] += 4.*PI/wilson_p("WPARAM_RUN_SM", 1)*(rh.get_matrix().V0[9-1][j]*coeffs_t[j]);
    // }

    // coeffs_b[9] = coeffs_t[9];

    // for (int i=0; i<coeffs.size(); i++) {
    //     this->at(coeffs[i])->set_WilsonCoeffRun("LO", coeffs_b[i]);
    // }
    // this->base["LO"] = 2;
}

void BCoefficientGroup::set_base_1_NLO() {
    // std::vector<std::string> coeffs {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"};
    // complex_t C7_eff= this->at("C7")->get_CoefficientMatchingValue("NLO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("NLO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    // -20./3.*this->at("C5")->get_CoefficientMatchingValue("NLO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 
	// complex_t C8_eff= this->at("C8")->get_CoefficientMatchingValue("NLO")+this->at("C3")->get_CoefficientMatchingValue("NLO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    // +20.*this->at("C5")->get_CoefficientMatchingValue("NLO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 

	// complex_t C7_eff_0= this->at("C7")->get_CoefficientMatchingValue("LO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("LO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("LO")
    // -20./3.*this->at("C5")->get_CoefficientMatchingValue("LO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("LO"); 
	// complex_t C8_eff_0= this->at("C8")->get_CoefficientMatchingValue("LO")+this->at("C3")->get_CoefficientMatchingValue("LO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("LO")
    // +20.*this->at("C5")->get_CoefficientMatchingValue("LO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("LO"); 

    // auto calculateC1b = [&](int ie, int je, BCoefficientGroup::iterator& iterator) {
    //     complex_t u0_term = (rh.get_matrix().U0)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") : (je == 6 ? C7_eff : C8_eff));
    //     complex_t u1_term = (rh.get_matrix().U1)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff_0 : C8_eff_0));
    //     return wilson_p("WPARAM_RUN_SM", 2) * (u0_term + u1_term);
    // };


    // BCoefficientGroup::iterator it = this->begin();
	// for (int ie = 0; ie < 8; ie++) {
    //     BCoefficientGroup::iterator ite = this->begin();
    //     complex_t _{};
    //     for (int je = 0; je < 8; je++) {
    //         _+= calculateC1b(ie, je,  ite);
    //     }
    //     this->find(coeffs[ie])->second->set_WilsonCoeffRun("NLO", _);
    // }

	// double fourPiOverAlphasMu = 4.0 * PI / wilson_p("WPARAM_RUN_SM", 1);

    // auto updateC1b = [&](int je, BCoefficientGroup::iterator& iterator) {
    //     return wilson_p("WPARAM_RUN_SM", 2) * ((rh.get_matrix().U0)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") + (rh.get_matrix().U1)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO"));
    // };

    // BCoefficientGroup::iterator iterator = this->begin();
    // complex_t _{};
    // for (int je = 0; je < 8; je++) {

    //     _ += fourPiOverAlphasMu * updateC1b(je, iterator);
    //     iterator++;
    // }

    // _ += fourPiOverAlphasMu * wilson_p("WPARAM_RUN_SM", 2) * (rh.get_matrix().U0)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("LO");
    // this->find("C9")->second->set_WilsonCoeffRun("NLO", _);
    // this->find("C10")->second->set_WilsonCoeffRun("NLO", wilson_p("WPARAM_RUN_SM", 2) * this->find("C10")->second->get_CoefficientMatchingValue("NLO"));
    // this->base["NLO"] = 1;
}

void BCoefficientGroup::set_base_2_NLO() {

    // std::vector<std::string> coeffs {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    // std::vector<complex_t> coeffs_t {};
    // std::vector<complex_t> coeffs_t_0 {};
    // // std::vector<complex_t> coeffs_b_0 {};
    // std::vector<complex_t> coeffs_b(10);

    // // std::cout << "WOOOOWOOOO" << std::endl;
    // coeffs_t_0.push_back(this->at("C1")->get_CoefficientMatchingValue("LO")/2.);
	// coeffs_t.push_back((-5.*this->at("C1")->get_CoefficientMatchingValue("LO")+3.*(-4.*this->at("C2")->get_CoefficientMatchingValue("LO")+this->at("C1")->get_CoefficientMatchingValue("NLO")))/6.);
	// coeffs_t_0.push_back(-this->at("C1")->get_CoefficientMatchingValue("LO")/6.+this->at("C2")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(-((11.*this->at("C1")->get_CoefficientMatchingValue("LO"))/18.)+(2.*this->at("C2")->get_CoefficientMatchingValue("LO"))/3.
    // -this->at("C1")->get_CoefficientMatchingValue("NLO")/6.+this->at("C2")->get_CoefficientMatchingValue("NLO"));
	// coeffs_t_0.push_back((54.*this->at("C3")->get_CoefficientMatchingValue("LO")-9.*this->at("C4")->get_CoefficientMatchingValue("LO")+864.*this->at("C5")->get_CoefficientMatchingValue("LO")
    // -144.*this->at("C6")->get_CoefficientMatchingValue("LO"))/54.); 
	// coeffs_t.push_back((36.*this->at("C3")->get_CoefficientMatchingValue("LO")-33.*this->at("C4")->get_CoefficientMatchingValue("LO")+4416.*this->at("C5")->get_CoefficientMatchingValue("LO")
    // -700.*this->at("C6")->get_CoefficientMatchingValue("LO")+54.*this->at("C3")->get_CoefficientMatchingValue("NLO")-9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    // +864.*this->at("C5")->get_CoefficientMatchingValue("NLO")-144.*this->at("C6")->get_CoefficientMatchingValue("NLO"))/54.);
	// coeffs_t_0.push_back(this->at("C4")->get_CoefficientMatchingValue("LO")/2.+8.*this->at("C6")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back((-36.*this->at("C3")->get_CoefficientMatchingValue("LO")-15.*this->at("C4")->get_CoefficientMatchingValue("LO")-960.*this->at("C5")->get_CoefficientMatchingValue("LO")
    // -644.*this->at("C6")->get_CoefficientMatchingValue("LO")+9.*this->at("C4")->get_CoefficientMatchingValue("NLO")+144.*this->at("C6")->get_CoefficientMatchingValue("NLO"))/18.);
	// coeffs_t_0.push_back((54.*this->at("C3")->get_CoefficientMatchingValue("LO")-9.*this->at("C4")->get_CoefficientMatchingValue("LO")+216.*this->at("C5")->get_CoefficientMatchingValue("LO")
    // -36.*this->at("C6")->get_CoefficientMatchingValue("LO"))/54.);
   	// coeffs_t.push_back((-36.*this->at("C3")->get_CoefficientMatchingValue("LO")+33.*this->at("C4")->get_CoefficientMatchingValue("LO")-4080.*this->at("C5")->get_CoefficientMatchingValue("LO")
    // -40.*this->at("C6")->get_CoefficientMatchingValue("LO")+54.*this->at("C3")->get_CoefficientMatchingValue("NLO")-9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    // +216.*this->at("C5")->get_CoefficientMatchingValue("NLO")-36.*this->at("C6")->get_CoefficientMatchingValue("NLO"))/54.);
 	// coeffs_t_0.push_back(this->at("C4")->get_CoefficientMatchingValue("LO")/2.+2.*this->at("C6")->get_CoefficientMatchingValue("LO"));
 	// coeffs_t.push_back((36.*this->at("C3")->get_CoefficientMatchingValue("LO")+15.*this->at("C4")->get_CoefficientMatchingValue("LO")+624.*this->at("C5")->get_CoefficientMatchingValue("LO")
    // +808.*this->at("C6")->get_CoefficientMatchingValue("LO")+9.*this->at("C4")->get_CoefficientMatchingValue("NLO")+36.*this->at("C6")->get_CoefficientMatchingValue("NLO"))/18.);
	// coeffs_t_0.push_back(this->at("C7")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(this->at("C7")->get_CoefficientMatchingValue("NLO"));
	// coeffs_t_0.push_back(this->at("C8")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(this->at("C8")->get_CoefficientMatchingValue("NLO"));
	// coeffs_t_0.push_back(this->at("C9")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(this->at("C9")->get_CoefficientMatchingValue("NLO"));
	// coeffs_t_0.push_back(this->at("C10")->get_CoefficientMatchingValue("LO"));
	// coeffs_t.push_back(this->at("C10")->get_CoefficientMatchingValue("NLO"));

    // complex_t C0t7= coeffs_t_0[6]-1./3.*coeffs_t_0[4]-coeffs_t_0[5]; 
	// complex_t C0t8= coeffs_t_0[7]+coeffs_t_0[4];

	// complex_t C1t7= coeffs_t[6]-1./3.*coeffs_t[4]-coeffs_t[5]; 
	// complex_t C1t8= coeffs_t[7]+coeffs_t[4];

    // for (int i=0; i<8; i++) {
    //     for (int j=0; j<8;j++) {
    //         if (j<6)
	// 	    {
    //             coeffs_b[i] += wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[i][j]*coeffs_t[j]+rh.get_matrix().V1[i][j]*coeffs_t_0[j]);
	// 	    }
	// 	    if (j==6)
	// 	    {
    //             coeffs_b[i] += wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[i][j]*C1t7+rh.get_matrix().V1[i][j]*C0t7);
	// 	    }
	// 	    if (j==7)
	// 	    {
    //             coeffs_b[i] += wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[i][j]*C1t8+rh.get_matrix().V1[i][j]*C0t8);
	// 	    }
    //     }
    // }

    // for (int j=0; j<8; j++) {
    //     coeffs_b[8] += 4.*PI/wilson_p("WPARAM_RUN_SM", 1)*(wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[9-1][j]*coeffs_t[j]+rh.get_matrix().V1[9-1][j]*coeffs_t_0[j]));
    // }

	// coeffs_b[8] += 4.*PI/wilson_p("WPARAM_RUN_SM", 1)*(wilson_p("WPARAM_RUN_SM", 2)*(rh.get_matrix().V0[9-1][9-1]*coeffs_t_0[8]));

	// coeffs_b[9]=wilson_p("WPARAM_RUN_SM", 2)*coeffs_t[9];

    // for (int i=0; i<coeffs.size(); i++) {
    //     this->at(coeffs[i])->set_WilsonCoeffRun("NLO", coeffs_b[i]);
    // }
    
    // this->base["NLO"] = 2;
}

void BCoefficientGroup::set_base_1_NNLO() {
    // std::vector<std::string> coeffs {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"};
    // complex_t C7_eff= this->at("C7")->get_CoefficientMatchingValue("NNLO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("NNLO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("NNLO")
    // -20./3.*this->at("C5")->get_CoefficientMatchingValue("NNLO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("NNLO"); 
	// complex_t C8_eff= this->at("C8")->get_CoefficientMatchingValue("NNLO")+this->at("C3")->get_CoefficientMatchingValue("NNLO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("NNLO")
    // +20.*this->at("C5")->get_CoefficientMatchingValue("NNLO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("NNLO"); 

	// complex_t C7_eff_0= this->at("C7")->get_CoefficientMatchingValue("LO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("LO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("LO")
    // -20./3.*this->at("C5")->get_CoefficientMatchingValue("LO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("LO"); 
	// complex_t C8_eff_0= this->at("C8")->get_CoefficientMatchingValue("LO")+this->at("C3")->get_CoefficientMatchingValue("LO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("LO")
    // +20.*this->at("C5")->get_CoefficientMatchingValue("LO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("LO"); 

	// complex_t C7_eff_1= this->at("C7")->get_CoefficientMatchingValue("NLO")-1./3.*this->at("C3")->get_CoefficientMatchingValue("NLO")-4./9.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    // -20./3.*this->at("C5")->get_CoefficientMatchingValue("NLO")-80./9.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 
	// complex_t C8_eff_1= this->at("C8")->get_CoefficientMatchingValue("NLO")+this->at("C3")->get_CoefficientMatchingValue("NLO")-1./6.*this->at("C4")->get_CoefficientMatchingValue("NLO")
    // +20.*this->at("C5")->get_CoefficientMatchingValue("NLO")-10./3.*this->at("C6")->get_CoefficientMatchingValue("NLO"); 

    // auto calculateC2b = [&](int ie, int je, BCoefficientGroup::iterator& iterator) {
    //     complex_t u0_term = (rh.get_matrix().U0)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NNLO") : (je == 6 ? C7_eff : C8_eff));
    //     complex_t u1_term = (rh.get_matrix().U1)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") : (je == 6 ? C7_eff_1 : C8_eff_1));
    //     complex_t u2_term = (rh.get_matrix().U2)[ie][je] * (je < 6 ? this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO") : (je == 6 ? C7_eff_0 : C8_eff_0));
    //     return wilson_p("WPARAM_RUN_SM", 2) * wilson_p("WPARAM_RUN_SM", 2) * (u0_term + u1_term + u2_term);
    // };

    // BCoefficientGroup::iterator it = this->begin();
	// for (int ie = 0; ie < 8; ie++) {
    //     BCoefficientGroup::iterator iterator = this->begin();
    //     complex_t _{};
    //     for (int je = 0; je < 8; je++) {
    //         _ += calculateC2b(ie, je, iterator);
    //     }
    //     this->find(coeffs[ie])->second->set_WilsonCoeffRun("NNLO", _);
    // }

	// double fourPiOverAlphasMu = 4.0 * PI / wilson_p("WPARAM_RUN_SM", 1);

    // auto updateC2b = [&](int je, BCoefficientGroup::iterator& iterator) {
    //     return wilson_p("WPARAM_RUN_SM", 2) * wilson_p("WPARAM_RUN_SM", 2) * ((rh.get_matrix().U0)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NNLO") 
    //     + (rh.get_matrix().U1)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("NLO") + (rh.get_matrix().U2)[8][je] * this->find(coeffs[je])->second->get_CoefficientMatchingValue("LO"));
    // };

    // BCoefficientGroup::iterator iterator = this->begin();
    // complex_t _{};
    // for (int je = 0; je < 8; je++) {
    //     _ += fourPiOverAlphasMu * updateC2b(je, iterator);
    //     iterator++;
    // }

    // _ += fourPiOverAlphasMu * wilson_p("WPARAM_RUN_SM", 2) * wilson_p("WPARAM_RUN_SM", 2) * ((rh.get_matrix().U0)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("NLO") + (rh.get_matrix().U1)[8][8] * this->find("C9")->second->get_CoefficientMatchingValue("LO"));
    // this->find("C9")->second->set_WilsonCoeffRun("NNLO", _);
    // this->find("C10")->second->set_WilsonCoeffRun("NNLO",wilson_p("WPARAM_RUN_SM", 2) * wilson_p("WPARAM_RUN_SM", 2) * this->find("C10")->second->get_CoefficientMatchingValue("NNLO"));
    // this->base["NNLO"] = 1;
}

void BCoefficientGroup::set_base_2_NNLO() {
    
}

void BCoefficientGroup::init_running_blocks() {
    WilsonParamComposer composer;

    LOG_DEBUG("Init running matrices blocks of B Coefficient group");
	std::unordered_map<ParameterType, std::vector<std::string>> eta_powers_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM"}}};

    auto eta_powers_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();

		for (int i = 0; i < BWilsonRunningParameters::array_size; ++i) {
            dep_block->store_or_assign(LhaID(1, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS", LhaID(1, i)}, std::pow(eta, (BWilsonRunningParameters::ai)[i]), 0., 0.));
            dep_block->store_or_assign(LhaID(2, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS", LhaID(2, i)}, std::pow(eta, (BWilsonRunningParameters::ai2)[i]), 0., 0.));
		}
    };

    std::unordered_map<ParameterType, std::vector<std::string>> mtx_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM", "ETA_POWS"}}};

    auto U_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        auto pid = [] (int n, int k, int l) {
            return ParamId{ParameterType::WILSON, "U_MATRIX", LhaID(n, k, l)};
        };

        double U0, U1, U2;
        double eta_ai;
        using BRP = BWilsonRunningParameters;
		for (int ke = 0; ke < BRP::array_size; ++ke) {
            for (int le = 0; le < BRP::array_size; ++le) {
                U0 = U1 = U2 = 0;
                for (int ie = 0; ie < BRP::array_size; ++ie) {
                    eta_ai = src.at("ETA_POWS")->retrieve(LhaID(1, ie))->get_val();
                    U0 += BRP::m00[ke][le][ie] * eta_ai;    
                    U1 += (BRP::m10[ke][le][ie] + BRP::m11[ke][le][ie] / eta) * eta_ai;
                    U2 += (BRP::m20[ke][le][ie] + BRP::m21[ke][le][ie] / eta + BRP::m22[ke][le][ie] / (eta * eta)) * eta_ai;
                }
                dep_block->store_or_assign(LhaID(0, ke, le), std::make_shared<Parameter>(pid(0, ke, le), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, ke, le), std::make_shared<Parameter>(pid(1, ke, le), U1, 0., 0.));
                dep_block->store_or_assign(LhaID(2, ke, le), std::make_shared<Parameter>(pid(2, ke, le), U2, 0., 0.));
            }
        }
    };

    auto V_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        auto pid = [] (int n, int k, int l) {
            return ParamId{ParameterType::WILSON, "V_MATRIX", LhaID(n, k, l)};
        };

        double V0, V1;
        double eta_ai;
        using BRP = BWilsonRunningParameters;
		for (int ke = 0; ke < BRP::array_size; ++ke) {
            for (int le = 0; le < BRP::array_size; ++le) {
                V0 = V1 = 0;
                for (int ie = 0; ie < BRP::array_size; ++ie) {
                    eta_ai = src.at("ETA_POWS")->retrieve(LhaID(1, ie))->get_val();
                    V0 += BRP::l00[ke][le][ie] * eta_ai;    
                    V1 += (BRP::l10[ke][le][ie] + BRP::l11[ke][le][ie] / eta) * eta_ai;
                }
                dep_block->store_or_assign(LhaID(0, ke, le), std::make_shared<Parameter>(pid(0, ke, le), V0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, ke, le), std::make_shared<Parameter>(pid(1, ke, le), V1, 0., 0.));
            }
        }
    };

    composer.compose_block("ETA_POWS", eta_powers_src, eta_powers_func);
    composer.compose_block("U_MATRIX", mtx_src, U_func);
    composer.compose_block("V_MATRIX", mtx_src, V_func);

    LOG_INFO("Running matrices updated");
}

void BScalarCoefficientGroup::set_base_1_LO() {
    // complex_t coeff_temp= this->at("CQ1")->get_CoefficientMatchingValue("LO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5));
    // this->at("CQ1")->set_WilsonCoeffRun("LO", coeff_temp);
    // complex_t coeff_temp2= this->at("CQ2")->get_CoefficientMatchingValue("LO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5));
    // this->at("CQ2")->set_WilsonCoeffRun("LO", coeff_temp2);

}

void BScalarCoefficientGroup::set_base_1_NLO() {
    // complex_t coeff_temp= this->at("CQ1")->get_CoefficientMatchingValue("NLO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5))*wilson_p("WPARAM_RUN_SM", 2);
    // this->at("CQ1")->set_WilsonCoeffRun("NLO", coeff_temp);
    // complex_t coeff_temp2= this->at("CQ2")->get_CoefficientMatchingValue("NLO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5))*wilson_p("WPARAM_RUN_SM", 2);
    // this->at("CQ2")->set_WilsonCoeffRun("NLO", coeff_temp2);
    // // std::cout << coeff_temp2 << std::endl;
}

void BPrimeCoefficientGroup::set_base_1_LO() {
    // complex_t coeff_temp= this->at("CP7")->get_CoefficientMatchingValue("LO")* std::pow(wilson_p("WPARAM_RUN_SM", 2), 16. / 23.) ;
    // this->at("CP7")->set_WilsonCoeffRun("LO", coeff_temp);
    // complex_t coeff_temp2= this->at("CP8")->get_CoefficientMatchingValue("LO")* std::pow(wilson_p("WPARAM_RUN_SM", 2), 14. / 23.);
    // this->at("CP8")->set_WilsonCoeffRun("LO", coeff_temp2);

    // complex_t coeff_temp5= this->at("CP9")->get_CoefficientMatchingValue("LO");
    // this->at("CP9")->set_WilsonCoeffRun("LO", coeff_temp5);
    // complex_t coeff_temp6= this->at("CP10")->get_CoefficientMatchingValue("LO");
    // this->at("CP10")->set_WilsonCoeffRun("LO", coeff_temp6);

    // complex_t coeff_temp3= this->at("CPQ1")->get_CoefficientMatchingValue("LO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5));
    // this->at("CPQ1")->set_WilsonCoeffRun("LO", coeff_temp3);
    // complex_t coeff_temp4= this->at("CPQ2")->get_CoefficientMatchingValue("LO")* pow(wilson_p("WPARAM_RUN_SM", 2),-4./wilson_p("WPARAM_SI_SM", 5));
    // this->at("CPQ2")->set_WilsonCoeffRun("LO", coeff_temp4);

}

void CoefficientGroup::claim_coefficients() {
    for (auto& [_, coeff]: *this) {
        coeff->set_owned(true);
    }
}
