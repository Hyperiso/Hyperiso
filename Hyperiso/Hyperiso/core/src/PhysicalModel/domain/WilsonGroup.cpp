#include "WilsonGroup.h"

void BCoefficientGroup::init_running_block(QCDOrder order, BWilsonBasis basis) {
    LOG_INFO("In BCoefficientGroup::init_running_block");

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {"B_MATCH", "WPARAM_RUN_SM"}},
    };

    if (basis == BWilsonBasis::STANDARD) {
        src.at(ParameterType::WILSON).push_back("WPARAM_MATCH_SM");
        src.at(ParameterType::WILSON).push_back("B_SCALE");
        src.at(ParameterType::WILSON).push_back("U_MATRIX");
        src.emplace(std::make_pair<ParameterType, std::vector<std::string>>(ParameterType::SM, {"SMINPUTS", "MASS"}));

        auto func = [this, order] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
            switch (order) {
            case QCDOrder::NNLO:
                BCoefficientGroup::base_1_NNLO_calculation(src, dep_block);
            case QCDOrder::NLO:
                BCoefficientGroup::base_1_NLO_calculation(src, dep_block);
            case QCDOrder::LO:
                BCoefficientGroup::base_1_LO_calculation(src, dep_block);
            }
        };
        WilsonParamComposer().compose_block("B_HADRONIC", src, func);
        basis = BWilsonBasis::STANDARD;
    } else {
        src.at(ParameterType::WILSON).push_back("V_MATRIX");

        auto func = [this, order] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
            switch (order) {
            case QCDOrder::NNLO:
                LOG_WARN("NNLO running is undefined in teh traditional basis of B Wilson coefficients.");
            case QCDOrder::NLO:
                BCoefficientGroup::base_2_NLO_calculation(src, dep_block);
            case QCDOrder::LO:
                BCoefficientGroup::base_2_LO_calculation(src, dep_block);
            }
        };
        WilsonParamComposer().compose_block("B_HADRONIC", src, func);
        basis = BWilsonBasis::TRADITIONAL;
    }
}

void BCoefficientGroup::switch_basis() {
    WilsonParamComposer().remove_block("B_HADRONIC");
    init_running_block(current_order, basis == BWilsonBasis::STANDARD ? BWilsonBasis::TRADITIONAL : BWilsonBasis::STANDARD);
}

void BCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block)
{
    auto ensure_coef = [src] (const LhaID& id) -> complex_t {
        return src.at("B_MATCH")->contains(id) ? src.at("B_MATCH")->retrieve(id)->get_val() : complex_t(0);
    };

    auto U0 = [src] (int k, int l) -> double {
        return src.at("U_MATRIX")->retrieve(LhaID(0, k, l))->get_val();
    };

    std::array<complex_t, 10> Ci_match = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 10; k++) {
        Ci_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::LO, false));
    }

    Ci_match[6] = BRP::C7_eff_std(Ci_match); 
    Ci_match[7] = BRP::C8_eff_std(Ci_match);

    std::array<complex_t, 10> Ci_run {};

    // C1 - C9
    for (size_t k = 0; k < 9; k++) {
        for (size_t l = 0; l < 9; l++) {
            Ci_run[k] += U0(k, l) * Ci_match[l];
        }
        LOG_INFO("C_run_", k + 1, "=", Ci_run[k]);
    }

    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] *= fact;

    // C10
    // Ask Nazila : Pourquoi les N et NNLO dans le running de C10 à LO et pas séparé en plusieurs ordres ?
    double alpha_ew = 1 / src.at("SMINPUTS")->retrieve(1)->get_val();
    double m_h = src.at("MASS")->retrieve(25)->get_val();
    double m_t_muW = src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val();
    double sw2OS = src.at("SMINPUTS")->retrieve(LhaID(7, 2))->get_val();
    double mu_h = src.at("B_SCALE")->retrieve(1)->get_val();
    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();

    complex_t C1_NLO = ensure_coef(WCoefMapper::flha_full(WCoef::C1, QCDOrder::NLO, false));
    complex_t C4_NLO = ensure_coef(WCoefMapper::flha_full(WCoef::C4, QCDOrder::NLO, false));

    complex_t C10_02=0.;
    for(size_t ie = 0; ie < 8; ie++) {
        C10_02 += BRP::b[ie] * pow(eta, BRP::a[ie]) * Ci_match[1];
    } 
    
    complex_t C10_12 = -0.11060 * std::log(eta) / eta * Ci_match[1]  + (1 / eta - 1) * (0.26087 * Ci_match[8] + 1.15942 * Ci_match[9]);
    for(int ie=0;ie<=7;ie++) {
        C10_12 += pow(eta, BRP::a[ie] + 1.) * (
            (BRP::d_2a[ie] / eta + BRP::d_2b[ie]) * Ci_match[1] 
            + BRP::d_1[ie] * C1_NLO
            + BRP::d_4[ie] * C4_NLO);
    }

    double Delta_alpha = 0.06; 
    double Delta_rhosw2 = -0.03; 
    double Delta_rem = 0.01;
    double Deltar = Delta_alpha + Delta_rhosw2 + Delta_rem; 

    double Gmu1_Gmu0 = 4 * PI / alpha_ew * Deltar;
    double L = std::log(mu_h * mu_h);
    complex_t C1022 = (46.9287715663914 - 3.102350691200236 * L + 0.0992974073578769 * L * L + 0.175877 * (m_t_muW - 163.5) + 0.0173725 * (m_h - 125.9)) / sw2OS;
    C1022 += -Ci_match[9] * Gmu1_Gmu0;

    complex_t C9_NLO = ensure_coef(WCoefMapper::flha_full(WCoef::C9, QCDOrder::NLO, false));
    complex_t C10_NLO = ensure_coef(WCoefMapper::flha_full(WCoef::C10, QCDOrder::NLO, false));
    
    complex_t C10_22 = (0.27924 * C1_NLO + 0.33157 * C4_NLO + 2.35917 * Ci_match[8] + 3.29679 * Ci_match[9]) * log(eta) + (1 - eta) * (0.26087 * C9_NLO + 1.15942 * C10_NLO) + C1022;
    
    complex_t C1_NNLO = ensure_coef(WCoefMapper::flha_full(WCoef::C1, QCDOrder::NNLO, false));
    complex_t C2_NNLO = ensure_coef(WCoefMapper::flha_full(WCoef::C2, QCDOrder::NNLO, false));
    complex_t C3_NNLO = ensure_coef(WCoefMapper::flha_full(WCoef::C3, QCDOrder::NNLO, false));
    complex_t C4_NNLO = ensure_coef(WCoefMapper::flha_full(WCoef::C4, QCDOrder::NNLO, false));
    complex_t C5_NNLO = ensure_coef(WCoefMapper::flha_full(WCoef::C5, QCDOrder::NNLO, false));
    complex_t C6_NNLO = ensure_coef(WCoefMapper::flha_full(WCoef::C6, QCDOrder::NNLO, false));

    for(int ie = 0; ie <= 7; ie++) {
        C10_22 += pow(eta, BRP::a[ie] + 2) * (
                  (BRP::e_1a[ie] / eta + BRP::e_1b[ie]) * C1_NLO
                + (BRP::e_4a[ie] / eta + BRP::e_4b[ie]) * C4_NLO
                + BRP::e_1[ie] * C1_NNLO
                + BRP::e_2[ie] * C2_NNLO
                + BRP::e_3[ie] * C3_NNLO
                + BRP::e_4[ie] * C4_NNLO
                + BRP::e_5[ie] * C5_NNLO
                + BRP::e_6[ie] * C6_NNLO);
    }

    double alpha_s_mu_h = src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[9] = Ci_match[9] + alpha_ew / alpha_s_mu_h * (fact * C10_02 + C10_12) + alpha_ew / (4 * PI) * C10_22;
    
    // Store
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {ParameterType::WILSON, "B_HADRONIC", WCoefMapper::flha_full(ids[k], QCDOrder::LO, false)};
        LOG_INFO("Storing coefficient", WCoefMapper::flha_full(ids[k], QCDOrder::LO, false));
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, Ci_run[k], 0., 0.));
    }
}

void BCoefficientGroup::base_2_LO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block)
{
    auto ensure_coef = [src] (const LhaID& id) -> complex_t {
        return src.at("B_MATCH")->contains(id) ? src.at("B_MATCH")->retrieve(id)->get_val() : complex_t(0);
    };

    auto V0 = [src] (int k, int l) -> double {
        return src.at("V_MATRIX")->retrieve(LhaID(0, k, l))->get_val();
    };

    std::array<complex_t, 10> Ci_match_trad = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 8; k++) {
        for (size_t j = 0; j < 8; j++) {
            Ci_match_trad[k] = BRP::std_to_trad_LO[k][j] * ensure_coef(WCoefMapper::flha_full(ids[j], QCDOrder::LO, false));
        }
    }
    Ci_match_trad[8] = ensure_coef(WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, false));

    std::array<complex_t, 10> Ci_run {};
    for (size_t k = 0; k < 9; k++) {
        for (size_t l = 0; l < 9; l++) {
            Ci_run[k] += V0(k, l) * Ci_match_trad[l];
        }
        LOG_INFO("C_run_", k + 1, "=", Ci_run[k]);
    }

    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] *= fact;
    Ci_run[9] = ensure_coef(WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, false));

    // Store
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {ParameterType::WILSON, "B_HADRONIC", WCoefMapper::flha_full(ids[k], QCDOrder::LO, false)};
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, Ci_run[k], 0., 0.));;
    }
}

void BCoefficientGroup::base_1_NLO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block)
{
    auto ensure_coef = [src] (const LhaID& id) -> complex_t {
        return src.at("B_MATCH")->contains(id) ? src.at("B_MATCH")->retrieve(id)->get_val() : complex_t(0);
    };

    auto U0 = [src] (int k, int l) -> double {
        return src.at("U_MATRIX")->retrieve(LhaID(0, k, l))->get_val();
    };

    auto U1 = [src] (int k, int l) -> double {
        return src.at("U_MATRIX")->retrieve(LhaID(1, k, l))->get_val();
    };

    auto C_eff = [] (const std::array<complex_t, 10>& Ci, const double* vec) -> complex_t {
        complex_t C {0};
        for (size_t k = 0; k < 8; k++) {
            C += Ci[k] * vec[k];
        }
        return C;
    };

    std::array<complex_t, 10> Ci_0_match = {};
    std::array<complex_t, 10> Ci_1_match = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 10; k++) {
        Ci_0_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::LO, false));
        Ci_1_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::NLO, false));
    }

    Ci_0_match[6] = C_eff(Ci_0_match, BRP::y_std);
    Ci_0_match[7] = C_eff(Ci_0_match, BRP::z_std);
    Ci_1_match[6] = C_eff(Ci_1_match, BRP::y_std); 
    Ci_1_match[7] = C_eff(Ci_1_match, BRP::z_std);

    std::array<complex_t, 10> Ci_run {};

    // C1 - C9
    for (size_t k = 0; k < 9; k++) {
        for (size_t l = 0; l < 9; l++) {
            Ci_run[k] += U0(k, l) * Ci_1_match[l] + U1(k, l) * Ci_0_match[l];
        }
        LOG_INFO("C_run_", k + 1, "=", Ci_run[k]);
    }

    // C9 special treatment
    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] = fact * (Ci_run[8] + U0(8, 8) * Ci_0_match[8]);

    // C10
    Ci_run[9] = ensure_coef(WCoefMapper::flha_full(WCoef::C10, QCDOrder::NLO, false));
    
    // Store
    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {ParameterType::WILSON, "B_HADRONIC", WCoefMapper::flha_full(ids[k], QCDOrder::NLO, false)};
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, eta * Ci_run[k], 0., 0.));
    }
}

void BCoefficientGroup::base_2_NLO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block)
{
    auto ensure_coef = [src] (const LhaID& id) -> complex_t {
        return src.at("B_MATCH")->contains(id) ? src.at("B_MATCH")->retrieve(id)->get_val() : complex_t(0);
    };

    auto V0 = [src] (int k, int l) -> double {
        return src.at("V_MATRIX")->retrieve(LhaID(0, k, l))->get_val();
    };

    auto V1 = [src] (int k, int l) -> double {
        return src.at("V_MATRIX")->retrieve(LhaID(1, k, l))->get_val();
    };

    std::array<complex_t, 10> Ci_0_match_trad = {};
    std::array<complex_t, 10> Ci_1_match_trad = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 8; k++) {
        for (size_t j = 0; j < 8; j++) {
            complex_t Cj_0_match_std = ensure_coef(WCoefMapper::flha_full(ids[j], QCDOrder::LO, false));
            complex_t Cj_1_match_std = ensure_coef(WCoefMapper::flha_full(ids[j], QCDOrder::NLO, false));
            Ci_0_match_trad[k] = BRP::std_to_trad_LO[k][j] * Cj_0_match_std;
            Ci_1_match_trad[k] = BRP::std_to_trad_LO[k][j] * Cj_1_match_std + BRP::std_to_trad_NLO[k][j] * Cj_0_match_std;
        }
    }
    Ci_0_match_trad[8] = ensure_coef(WCoefMapper::flha_full(WCoef::C9, QCDOrder::LO, false));
    Ci_1_match_trad[8] = ensure_coef(WCoefMapper::flha_full(WCoef::C9, QCDOrder::NLO, false));

    std::array<complex_t, 10> Ci_run {};
    for (size_t k = 0; k < 9; k++) {
        for (size_t l = 0; l < 8; l++) {
            Ci_run[k] += V0(k, l) * Ci_1_match_trad[l] + V1(k, l) * Ci_0_match_trad[l];
        }
        LOG_INFO("C_run_", k + 1, "=", Ci_run[k]);
    }

    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] = fact * (Ci_run[8] + V0(8, 8) * Ci_0_match_trad[8]);
    Ci_run[9] = ensure_coef(WCoefMapper::flha_full(WCoef::C10, QCDOrder::LO, false));

    // Store
    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {ParameterType::WILSON, "B_HADRONIC", WCoefMapper::flha_full(ids[k], QCDOrder::NLO, false)};
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, eta * Ci_run[k], 0., 0.));
    }
}

void BCoefficientGroup::base_1_NNLO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block)
{
    auto ensure_coef = [src] (const LhaID& id) -> complex_t {
        return src.at("B_MATCH")->contains(id) ? src.at("B_MATCH")->retrieve(id)->get_val() : complex_t(0);
    };

    auto U = [src] (int order, int k, int l) -> double {
        return src.at("U_MATRIX")->retrieve(LhaID(order, k, l))->get_val();
    };

    std::array<complex_t, 10> Ci_0_match = {};
    std::array<complex_t, 10> Ci_1_match = {};
    std::array<complex_t, 10> Ci_2_match = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 10; k++) {
        Ci_0_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::LO, false));
        Ci_1_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::NLO, false));
        Ci_2_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::NNLO, false));
    }

    Ci_0_match[6] = BRP::C7_eff_std(Ci_0_match);
    Ci_0_match[7] = BRP::C8_eff_std(Ci_0_match);
    Ci_1_match[6] = BRP::C7_eff_std(Ci_1_match);
    Ci_1_match[7] = BRP::C8_eff_std(Ci_1_match);
    Ci_2_match[6] = BRP::C7_eff_std(Ci_2_match);
    Ci_2_match[7] = BRP::C8_eff_std(Ci_2_match);

    std::array<complex_t, 10> Ci_run {};

    // C1 - C9
    for (size_t k = 0; k < 9; k++) {
        for (size_t l = 0; l < 9; l++) {
            Ci_run[k] += U(2, k, l) * Ci_0_match[l] + U(1, k, l) * Ci_1_match[l] + U(0, k, l) * Ci_2_match[l];
        }
        LOG_INFO("C_run_", k + 1, "=", Ci_run[k]);
    }

    // C9 special treatment
    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] = fact * (Ci_run[8] + U(1, 8, 8) * Ci_0_match[8] + U(0, 8, 8) * Ci_1_match[8]);

    // C10
    Ci_run[9] = ensure_coef(WCoefMapper::flha_full(WCoef::C10, QCDOrder::NLO, false));
    
    // Store
    double eta_sq = pow(src.at("WPARAM_RUN_SM")->retrieve(2)->get_val(), 2);
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {ParameterType::WILSON, "B_HADRONIC", WCoefMapper::flha_full(ids[k], QCDOrder::NLO, false)};
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, eta_sq * Ci_run[k], 0., 0.));
    }
}


void BCoefficientGroup::init_running_parameter_blocks() {
    WilsonParamComposer composer;

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

    composer.compose_block("ETA_POWS", eta_powers_src, eta_powers_func);
    composer.compose_block("U_MATRIX", mtx_src, U_func);
    composer.compose_block("V_MATRIX", mtx_src, V_func);

    LOG_INFO("Running matrices updated");
}

void BScalarCoefficientGroup::set_base_1_LO() {
    LOG_INFO("In BScalarCoefficientGroup::set_base_1_LO");

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {"B_SCALAR_MATCH", "WPARAM_RUN_SM", "WPARAM_SI_SM"}},
    };

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        auto ensure_coef = [src] (const LhaID& id) -> complex_t {
            return src.at("B_SCALAR_MATCH")->contains(id) ? src.at("B_SCALAR_MATCH")->retrieve(id)->get_val() : complex_t(0);
        };

        std::array<complex_t, 10> CQi_match = {};
        auto ids = WCoefMapper::get_group(WGroup::BScalar);
        for (size_t k = 0; k < ids.size(); k++) {
            CQi_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::LO, false));
        }

        double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        double beta_0 = src.at("WPARAM_SI_SM")->retrieve(5)->get_val(); // TODO : change to QCD params
        double fact = pow(eta, -4 / beta_0);
        
        // Store
        for (size_t k = 0; k < ids.size(); k++) {
            ParamId pid {ParameterType::WILSON, "B_SCALAR_HADRONIC", WCoefMapper::flha_full(ids[k], QCDOrder::LO, false)};
            dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, fact * CQi_match[k], 0., 0.));;
        }
    };

    WilsonParamComposer().compose_block("B_SCALAR_HADRONIC", src, func);
}

void BScalarCoefficientGroup::set_base_1_NLO() {
    LOG_INFO("In BScalarCoefficientGroup::set_base_1_NLO");

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {"B_SCALAR_MATCH", "WPARAM_RUN_SM", "WPARAM_SI_SM"}},
    };

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        auto ensure_coef = [src] (const LhaID& id) -> complex_t {
            return src.at("B_SCALAR_MATCH")->contains(id) ? src.at("B_SCALAR_MATCH")->retrieve(id)->get_val() : complex_t(0);
        };

        auto ids = WCoefMapper::get_group(WGroup::BScalar);
        std::array<complex_t, 2> CQi_match = {};
        for (size_t k = 0; k < ids.size(); k++) {
            CQi_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::NLO, false));
        }

        double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        double beta_0 = src.at("WPARAM_SI_SM")->retrieve(5)->get_val(); // TODO : change to QCD params
        double fact = pow(eta, 1 - 4 / beta_0);
        
        // Store
        for (size_t k = 0; k < ids.size(); k++) {
            ParamId pid {ParameterType::WILSON, "B_SCALAR_HADRONIC", WCoefMapper::flha_full(ids[k], QCDOrder::NLO, false)};
            dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, fact * CQi_match[k], 0., 0.));;
        }
    };

    WilsonParamComposer().compose_block("B_SCALAR_HADRONIC", src, func);
}

void BPrimeCoefficientGroup::set_base_1_LO() {
    LOG_INFO("In BPrimeCoefficientGroup::set_base_1_LO");

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {"B_PRIME_MATCH", "WPARAM_RUN_SM", "WPARAM_SI_SM"}},
    };

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        auto ensure_coef = [src] (const LhaID& id) -> complex_t {
            return src.at("B_PRIME_MATCH")->contains(id) ? src.at("B_PRIME_MATCH")->retrieve(id)->get_val() : complex_t(0);
        };

        std::array<complex_t, 12> CPi_match = {};
        auto ids = WCoefMapper::get_group(WGroup::BPrime);
        for (size_t k = 0; k < ids.size(); k++) {
            CPi_match[k] = ensure_coef(WCoefMapper::flha_full(ids[k], QCDOrder::LO, false));
        }

        double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        
        // Store
        for (size_t k = 0; k < ids.size(); k++) {
            ParamId pid {ParameterType::WILSON, "B_PRIME_HADRONIC", WCoefMapper::flha_full(ids[k], QCDOrder::LO, false)};
            complex_t CPk_run = pow(eta, BRP::exp_prime_running[k]) * CPi_match[k]; 
            dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, CPk_run, 0., 0.));;
        }
    };

    WilsonParamComposer().compose_block("B_PRIME_HADRONIC", src, func);
}

void CoefficientGroup::claim_coefficients() {
    for (auto& [_, coeff]: *this) {
        coeff->set_owned(true);
    }
}
