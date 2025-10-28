#include "WilsonParametersHelper.h"

void WilsonParameterHelper::init(int gen, WGroupId grp) {
	if (initialized) {
		std::cout << "wilson_param_helper already done " << std::endl;
		return;
	}
	std::cout << "Initializing WilsonParameterHelper" << std::endl;
	LOG_DEBUG("Initializing WilsonParameterHelper");
	init_scale_independent_block(gen);
	init_matching_block();
	init_running_block(grp);
	initialized = true;
}

void WilsonParameterHelper::init_scale_independent_block(int gen) {
	LOG_DEBUG("Init scale-independent wparam block");
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"SMINPUTS", "MASS"}}};

    auto func = [gen] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double xh = pow(src.at("MASS")->retrieve(25)->get_val() / src.at("MASS")->retrieve(24)->get_val(), 2);
        printf("25 in the SM (LO) : %.8lf\n", src.at("MASS")->retrieve(25)->get_val().real());
        printf("24 in the SM (LO) : %.8lf\n", src.at("MASS")->retrieve(24)->get_val().real());
		int nf = 5;
		int id = 1;
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, xh, 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, gen, 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, src.at("MASS")->retrieve(9 + 2 * gen)->get_val(), 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, src.at("SMINPUTS")->retrieve({7, 1})->get_val(), 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, 11.-2./3.*nf, 0., 0.)); //TODO, beta0
    };
    iblock_c->compose_block("WPARAM_SI_SM", src, func);
}

void WilsonParameterHelper::init_matching_block() {
	LOG_DEBUG("Init matching scale dependent wparam block");
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS", "QCD"}}, {ParameterType::WILSON, {"EW_SCALE"}}};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double mu_W = src.at("EW_SCALE")->retrieve(1)->get_val();
        double alphas_muW = QCDHelper::alpha_s(mu_W);
		double mass_top_muW = QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR);
		double mass_b_muW_mbrun = QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR);
		double mass_b_muW_mbpole = QCDHelper::msbar_mass(5, mu_W, MassType::POLE);
		double mass_c_muW = QCDHelper::msbar_mass(4, mu_W, MassType::POLE);
		printf("mass_b_muW : %.14lf\n",mass_b_muW_mbrun);
		double m_W = src.at("MASS")->retrieve(24)->get_val();
		double xt = pow(mass_top_muW / m_W, 2);
		double L = log(std::pow(mu_W / m_W, 2));
		double xtW = pow(QCDHelper::msbar_mass(6, m_W) / m_W, 2); // mass top at pole for mtot param
		double xtt = pow(src.at("QCD")->retrieve(6)->get_val() / m_W, 2.); // 24 -> W

		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 1}, alphas_muW, 0., 0.));
		dep_block->store_or_assign(LhaID(2, 1), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)}, xt, 0., 0.));
		dep_block->store_or_assign(LhaID(2, 2), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 2)}, std::pow(xt, 2), 0., 0.));
		dep_block->store_or_assign(LhaID(2, 3), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 3)}, std::pow(xt, 3), 0., 0.));
		dep_block->store_or_assign(LhaID(2, 4), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 4)}, std::pow(xt, 4), 0., 0.));
		dep_block->store_or_assign(3, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 3}, L, 0., 0.));
		dep_block->store_or_assign(4, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 4}, mass_c_muW, 0., 0.));
		dep_block->store_or_assign(LhaID(5, 1), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, mass_b_muW_mbrun, 0., 0.));
		dep_block->store_or_assign(LhaID(5, 2), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 2)}, mass_b_muW_mbpole, 0., 0.));
		dep_block->store_or_assign(6, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 6}, mass_top_muW, 0., 0.));
		dep_block->store_or_assign(7, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 7}, xtW, 0., 0.));
		dep_block->store_or_assign(8, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 8}, xtt, 0., 0.));

		LOG_DEBUG("Update matching block");
    };

	

    iblock_c->compose_block("WPARAM_MATCH_SM", src, func);


}

void WilsonParameterHelper::init_running_parameter_blocks_B() {

    LOG_DEBUG("Init running matrices blocks of B Coefficient group");
	std::unordered_map<ParameterType, std::vector<std::string>> eta_powers_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM"}}};

    auto eta_powers_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();

		for (int i = 0; i < BRP::array_size; ++i) {
            dep_block->store_or_assign(LhaID(1, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS", LhaID(1, i)}, std::pow(eta, (BRP::ai)[i]), 0., 0.));
            dep_block->store_or_assign(LhaID(2, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS", LhaID(2, i)}, std::pow(eta, (BRP::ai2)[i]), 0., 0.));
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
        using BRP = BRP;
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
        using BRP = BRP;
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

    iblock_c->compose_block("ETA_POWS", eta_powers_src, eta_powers_func);
    iblock_c->compose_block("U_MATRIX", mtx_src, U_func);
    iblock_c->compose_block("V_MATRIX", mtx_src, V_func);

    LOG_VERBOSE("Running matrices updated");
}

void WilsonParameterHelper::init_running_parameter_blocks_MM() {
    // WilsonParamComposer composer;

    LOG_DEBUG("Init running matrices blocks of Meson Mixing Coefficient group");
	std::unordered_map<ParameterType, std::vector<std::string>> eta_powers_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM"}}};

    auto eta_powers_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta5 = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        double eta4 = src.at("WPARAM_RUN_SM")->retrieve(3)->get_val();

		for (int i = 0; i < MMRP::n_pows; ++i) {
            dep_block->store_or_assign(LhaID(1, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS_MIXING", LhaID(1, i)}, std::pow(eta5, (MMRP::ai)[i]), 0., 0.));
            dep_block->store_or_assign(LhaID(2, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS_MIXING", LhaID(2, i)}, std::pow(eta4, (MMRP::bi)[i]), 0., 0.));
		}
    };

    std::unordered_map<ParameterType, std::vector<std::string>> mtx_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM", "ETA_POWS_MIXING"}}};

    auto U_5_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta_5 = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        auto pid = [] (int n, int k, int l) {
            return ParamId{ParameterType::WILSON, "UM_MATRIX_5", LhaID(n, k, l)};
        };

        using MMRP = MesonMixingRunningParameters;

        double U0, U1;

        // V
        double eta_V = src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 0))->get_val();
        U0 = MMRP::a0_V_5 * eta_V;
        U1 = (MMRP::a1_V_5 + MMRP::b_V_5 * eta_5) * eta_V;
        dep_block->store_or_assign(LhaID(0, 0, 0), std::make_shared<Parameter>(pid(0, 0, 0), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(0, 5, 5), std::make_shared<Parameter>(pid(0, 5, 5), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 0, 0), std::make_shared<Parameter>(pid(1, 0, 0), U1, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 5, 5), std::make_shared<Parameter>(pid(1, 5, 5), U1, 0., 0.));

        // LR
        std::array<double, 2> eta_LR {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 1))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 2))->get_val(),
        };

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                U0 = U1 = 0;
                for (int k = 0; k < 2; ++k) {
                    U0 += MMRP::a0_LR_5[i][j][k] * eta_LR[k];    
                    U1 += (MMRP::a1_LR_5[i][j][k] + MMRP::b_LR_5[i][j][k] * eta_5) * eta_LR[k];
                }
                dep_block->store_or_assign(LhaID(0, 1+i, 1+j), std::make_shared<Parameter>(pid(0, 1+i, 1+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 1+i, 1+j), std::make_shared<Parameter>(pid(1, 1+i, 1+j), U1, 0., 0.));
            }
        }

        // S
		std::array<double, 2> eta_S {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 3))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 4))->get_val(),
        };

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                U0 = U1 = 0;
                for (int k = 0; k < 2; ++k) {
                    U0 += MMRP::a0_S_5[i][j][k] * eta_S[k];    
                    U1 += (MMRP::a1_S_5[i][j][k] + MMRP::b_S_5[i][j][k] * eta_5) * eta_S[k];
                }
                dep_block->store_or_assign(LhaID(0, 3+i, 3+j), std::make_shared<Parameter>(pid(0, 3+i, 3+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 3+i, 3+j), std::make_shared<Parameter>(pid(1, 3+i, 3+j), U1, 0., 0.));
                dep_block->store_or_assign(LhaID(0, 6+i, 6+j), std::make_shared<Parameter>(pid(0, 6+i, 6+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 6+i, 6+j), std::make_shared<Parameter>(pid(1, 6+i, 6+j), U1, 0., 0.));
            }
        }
    };

    auto U_4_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta_5 = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        double eta_4 = src.at("WPARAM_RUN_SM")->retrieve(3)->get_val();
        auto pid = [] (int n, int k, int l) {
            return ParamId{ParameterType::WILSON, "UM_MATRIX_4", LhaID(n, k, l)};
        };

        using MMRP = MesonMixingRunningParameters;

        double U0, U1;

        // V
        double eta_5_V = src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 0))->get_val();
        double eta_4_V = src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 0))->get_val();
        U0 = MMRP::a0_V_4 * eta_4_V * eta_5_V;
        U1 = (MMRP::a1_V_4 + MMRP::b_V_4 * eta_4 + MMRP::c_V_4 * eta_4 * eta_5) * eta_4_V * eta_5_V;
        dep_block->store_or_assign(LhaID(0, 0, 0), std::make_shared<Parameter>(pid(0, 0, 0), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(0, 5, 5), std::make_shared<Parameter>(pid(0, 5, 5), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 0, 0), std::make_shared<Parameter>(pid(1, 0, 0), U1, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 5, 5), std::make_shared<Parameter>(pid(1, 5, 5), U1, 0., 0.));

        // LR
        std::array<double, 2> eta_5_LR {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 1))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 2))->get_val(),
        };

        std::array<double, 2> eta_4_LR {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 1))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 2))->get_val(),
        };

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                U0 = U1 = 0;
                for (int k = 0; k < 2; ++k) for (int l = 0; l < 2; ++l) {
                    U0 += MMRP::a0_LR_4[i][j][k][l] * eta_4_LR[k] * eta_5_LR[l];    
                    U1 += (MMRP::a1_LR_4[i][j][k][l] + MMRP::b_LR_4[i][j][k][l] * eta_4 + MMRP::c_LR_4[i][j][k][l] * eta_4 * eta_5) * eta_4_LR[k] * eta_5_LR[l];
                }
                dep_block->store_or_assign(LhaID(0, 1+i, 1+j), std::make_shared<Parameter>(pid(0, 1+i, 1+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 1+i, 1+j), std::make_shared<Parameter>(pid(1, 1+i, 1+j), U1, 0., 0.));
            }
        }

        // S
		std::array<double, 2> eta_5_S {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 3))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 4))->get_val(),
        };

        std::array<double, 2> eta_4_S {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 3))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 4))->get_val(),
        };

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                U0 = U1 = 0;
                for (int k = 0; k < 2; ++k) for (int l = 0; l < 2; ++l) {
                    U0 += MMRP::a0_S_4[i][j][k][l] * eta_4_S[k] * eta_5_S[l];    
                    U1 += (MMRP::a1_S_4[i][j][k][l] + MMRP::b_S_4[i][j][k][l] * eta_4 + MMRP::c_S_4[i][j][k][l] * eta_4 * eta_5) * eta_4_S[k] * eta_5_S[l];
                }
                dep_block->store_or_assign(LhaID(0, 3+i, 3+j), std::make_shared<Parameter>(pid(0, 3+i, 3+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 3+i, 3+j), std::make_shared<Parameter>(pid(1, 3+i, 3+j), U1, 0., 0.));
                dep_block->store_or_assign(LhaID(0, 6+i, 6+j), std::make_shared<Parameter>(pid(0, 6+i, 6+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 6+i, 6+j), std::make_shared<Parameter>(pid(1, 6+i, 6+j), U1, 0., 0.));
            }
        }
    };

    iblock_c->compose_block("ETA_POWS_MIXING", eta_powers_src, eta_powers_func);
    iblock_c->compose_block("UM_MATRIX_5", mtx_src, U_5_func);
    iblock_c->compose_block("UM_MATRIX_4", mtx_src, U_4_func);
}

void WilsonParameterHelper::init_running_block(WGroupId grp) {
	LOG_DEBUG("Init running scale dependent wparam block");
	std::unordered_map<ParameterType, std::vector<std::string>> src = {
		{ParameterType::WILSON, {"EW_SCALE", "B_SCALE"}},
		{ParameterType::SM, {"QCD", "MASS"}}
	};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double mu_W = src.at("EW_SCALE")->retrieve(1)->get_val();
		double mu_h = src.at("B_SCALE")->retrieve(1)->get_val();

		int n_f_final = QCDHelper::get_nf(mu_h, MassType::MSBAR);
        double alphas_mu_h = QCDHelper::alpha_s(mu_h);
		std::cout << "n_f_final" << n_f_final << std::endl;
		double eta_5 {0}, eta_4 {0}, eta_3 {0};
		if (n_f_final == 5) {
			eta_5 = QCDHelper::alpha_s(mu_W) / alphas_mu_h;
			double mu_b = src.at("QCD")->retrieve({5, 1})->get_val(); //TODO : check this part
			double alpha_s_mu_b = QCDHelper::alpha_s(mu_b, MassType::MSBAR);
			eta_4 = alpha_s_mu_b / alphas_mu_h;
		} else if (n_f_final == 4) {
			double mu_b = src.at("QCD")->retrieve({5, 1})->get_val();
			double alpha_s_mu_b = QCDHelper::alpha_s(mu_b, MassType::MSBAR);
			eta_5 = QCDHelper::alpha_s(mu_W) / alpha_s_mu_b;
			eta_4 = alpha_s_mu_b / alphas_mu_h;
		} else if (n_f_final == 3) {
			double mu_b = src.at("QCD")->retrieve({5, 1})->get_val();
			double mu_c = src.at("MASS")->retrieve(4)->get_val();
			double alpha_s_mu_b = QCDHelper::alpha_s(mu_b, MassType::MSBAR);
			double alpha_s_mu_c = QCDHelper::alpha_s(mu_c, MassType::MSBAR);
			eta_5 = QCDHelper::alpha_s(mu_W) / alpha_s_mu_b;
			eta_4 = alpha_s_mu_b / alpha_s_mu_c;
			std::cout << "eta_5" << eta_5 << std::endl;
			std::cout << "eta_4" << eta_4 << std::endl;
			eta_3 = alpha_s_mu_c / alphas_mu_h;
		} else {
			LOG_ERROR("ValueError", "In WilsonParameterHelper::init_running_block() : Hadronic mass is out of allowed range in Hyperiso.");
		}	
		
		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 1}, alphas_mu_h, 0., 0.));
		dep_block->store_or_assign(2, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 2}, eta_5, 0., 0.));
		dep_block->store_or_assign(3, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 3}, eta_4, 0., 0.));
		dep_block->store_or_assign(4, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 4}, eta_3, 0., 0.));
    };

    iblock_c->compose_block("WPARAM_RUN_SM", src, func);

	if (grp == GroupMapper::to_id(WGroup::B)) {
		init_running_parameter_blocks_B();

	} else if (grp == GroupMapper::to_id(WGroup::MESON_MIXING)) {
		init_running_parameter_blocks_MM();
	} else {
		init_running_parameter_blocks_B(); //TODO : better ?
	}
}

void WilsonParameterHelper::cleanup() {
	initialized = false;
}