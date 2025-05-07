#include "WilsonGroup.h"

CoefficientGroup::CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs) {
    this->insert(coeffs.begin(), coeffs.end());
    
    QCDOrder max_order = QCDOrder::LO;
    for (const auto& [_, wil] : coeffs) {
        if (wil->get_max_order() > max_order) {
            max_order = wil->get_max_order();
        }

        if (max_order == QCDOrder::NNLO) break; 
    }

    this->init(max_order == QCDOrder::NONE ? QCDOrder::LO : max_order);
}

void CoefficientGroup::init(QCDOrder order) {
    this->claim_coefficients();
    for (auto& coeff : *this) {
        LOG_DEBUG("Initializing Wilson coefficient", coeff.first, "at", OrderMapper::str(order));
        coeff.second->init(order);
    }
    this->current_order = order;
}

complex_t CoefficientGroup::get_matching_coefficient(std::string coeff, std::string order) const { 
    return this->at(coeff)->get_matching_value(order); 
}

complex_t CoefficientGroup::get_running_coefficient(std::string coeff, std::string order) const {
    auto coef = this->at(coeff);
    ParameterProxy wilson_p = ParameterProxy(ParameterType::WILSON);
    // std::cout << "before : " << coef->id(OrderMapper::enum_elt(order)) << std::endl;
    return complex_t(wilson_p(GroupMapper::str(this->id, ScaleType::HADRONIC, false, this->basis.value_or(BWilsonBasis::STANDARD)), coef->id(OrderMapper::enum_elt(order))));
}

QCDOrder CoefficientGroup::get_order(){
    return this->current_order;
}

bool CoefficientGroup::is_double_basis() const {
    return this->basis.has_value();
}

BCoefficientGroup::BCoefficientGroup() {
    LOG_TRACE("In BCoefficientGroup constructor");
    init_running_parameter_blocks();
    this->basis = BWilsonBasis::STANDARD;
    this->id = WGroup::B;

    if (UseMarty().get()) {
        for (auto&& coeff : {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"}) {
            this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
        }
        this->wilson_type = ContributionType::TOTAL;
        return;
    }

    this->insert(std::make_pair("C1", std::make_shared<C1>())); 
    this->insert(std::make_pair("C2", std::make_shared<C2>())); 
    this->insert(std::make_pair("C3", std::make_shared<C3>()));
    this->insert(std::make_pair("C4", std::make_shared<C4>()));  
    this->insert(std::make_pair("C5", std::make_shared<C5>())); 
    this->insert(std::make_pair("C6", std::make_shared<C6>())); 
    this->insert(std::make_pair("C7", std::make_shared<C7>()));  
    this->insert(std::make_pair("C8", std::make_shared<C8>()));  
    this->insert(std::make_pair("C9", std::make_shared<C9>())); 
    this->insert(std::make_pair("C10", std::make_shared<C10>())); 
}

std::shared_ptr<CoefficientGroup> BCoefficientGroup::clone() const {
    return std::make_shared<BCoefficientGroup>(*this);
}

void BCoefficientGroup::init_running_blocks(QCDOrder order) {
    LOG_INFO("In BCoefficientGroup::init_running_blocks");

    std::unordered_map<ParameterType, std::vector<std::string>> src_1 = {
        {ParameterType::WILSON, {GroupMapper::str(this->id, ScaleType::MATCHING), "WPARAM_RUN_SM", "WPARAM_MATCH_SM", "B_SCALE", "U_MATRIX"}},
        {ParameterType::SM, {"SMINPUTS", "MASS"}},
    };

    auto func_1 = [order, this] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        switch (order) {
        case QCDOrder::NNLO:
            BCoefficientGroup::base_1_NNLO_calculation(src, dep_block, this->wilson_type);
        case QCDOrder::NLO:
            BCoefficientGroup::base_1_NLO_calculation(src, dep_block, this->wilson_type);
        case QCDOrder::LO:
            BCoefficientGroup::base_1_LO_calculation(src, dep_block, this->wilson_type);
        }
    };
    
    WilsonParamComposer().compose_block(GroupMapper::str(this->id, ScaleType::HADRONIC, false, BWilsonBasis::STANDARD) + "INTER", src_1, func_1);

    std::unordered_map<ParameterType, std::vector<std::string>> src_2 = {
        {ParameterType::WILSON, {GroupMapper::str(this->id, ScaleType::MATCHING), "WPARAM_RUN_SM", "V_MATRIX"}}
    };

    auto func_2 = [order, this] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        switch (order) {
        case QCDOrder::NNLO:
            LOG_WARN("NNLO running is undefined in the traditional basis of B Wilson coefficients.");
        case QCDOrder::NLO:
            BCoefficientGroup::base_2_NLO_calculation(src, dep_block, this->wilson_type);
        case QCDOrder::LO:
            BCoefficientGroup::base_2_LO_calculation(src, dep_block, this->wilson_type);
        }
    };

    WilsonParamComposer().compose_block(GroupMapper::str(this->id, ScaleType::HADRONIC, false, BWilsonBasis::TRADITIONAL) + "INTER", src_2, func_2);

    init_full_running_block(BWilsonBasis::STANDARD);
    init_full_running_block(BWilsonBasis::TRADITIONAL);
}

void CoefficientGroup::init_full_running_block(BWilsonBasis basis) {
    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {GroupMapper::str(this->id, ScaleType::HADRONIC) + "INTER", "WPARAM_RUN_SM"}}
    };

    auto func = [basis, this] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double alpha_s_mu_h = src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
        double fact = alpha_s_mu_h / 4 * PI;
        
        for (WCoef coef_id : WCoefMapper::get_group(this->id)) {
            complex_t coef_full {0.};
            for (int order=0; order < 2; order++) {
                coef_full += ensure_coef(coef_id, (QCDOrder)(order + 1), this->wilson_type, GroupMapper::str(this->id, ScaleType::HADRONIC) + "INTER");
            }

            auto coef_lha_base = WCoefMapper::flha_base(coef_id);
            ParamId pid {
                ParameterType::WILSON, 
                GroupMapper::str(this->id, ScaleType::HADRONIC, true, basis) + "INTER", 
                LhaID(coef_lha_base.first, coef_lha_base.second, (int)this->wilson_type)
            };
            // LOG_INFO("Storing full coefficient", LhaID(coef_lha_base.first, coef_lha_base.second, (int)this->wilson_type));
            dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, coef_full, 0., 0.));
        }
    };

    WilsonParamComposer().compose_block(GroupMapper::str(this->id, ScaleType::HADRONIC, true, basis) + "INTER", src, func);
}

void CoefficientGroup::switch_basis() {
    if (!basis.has_value()) return;

    basis = (basis.value() == BWilsonBasis::STANDARD ? BWilsonBasis::TRADITIONAL : BWilsonBasis::STANDARD);
}

void BCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block, 
    ContributionType type)
{
    LOG_DEBUG("Init LO running of BCoefficientGroup in standard basis");

    auto U0 = [src] (int k, int l) -> double {
        return src.at("U_MATRIX")->retrieve(LhaID(0, k, l))->get_val();
    };

    std::array<complex_t, 10> Ci_match = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 10; k++) {
        Ci_match[k] = ensure_coef(ids[k], QCDOrder::LO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    }

    Ci_match[6] = BRP::C7_eff_std(Ci_match); 
    Ci_match[7] = BRP::C8_eff_std(Ci_match);

    std::array<complex_t, 10> Ci_run {};

    // C1 - C9
    for (size_t k = 0; k < 9; k++) {
        for (size_t l = 0; l < 9; l++) {
            Ci_run[k] += U0(k, l) * Ci_match[l];
        }
        LOG_VERBOSE("C_match_", k + 1, "=", Ci_match[k]);
        LOG_VERBOSE("C_run_", k + 1, "=", Ci_run[k]);
    }

    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] *= fact;

    // C10
    // Ask Nazila : Pourquoi les N et NNLO dans le running de C10 à LO et pas séparé en plusieurs ordres ?
    // Answer : Ask Siavash, need to understand if pb comes from litt or implementation.
    // Answer : Fck around and find out
    double alpha_ew = 1. / src.at("SMINPUTS")->retrieve(1)->get_val();
    double m_h = src.at("MASS")->retrieve(25)->get_val();
    double m_t_muW = src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val();
    double sw2OS = src.at("SMINPUTS")->retrieve(LhaID(7, 2))->get_val();
    double mu_h = src.at("B_SCALE")->retrieve(1)->get_val();
    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();

    complex_t C1_NLO = ensure_coef(WCoef::C1, QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    complex_t C4_NLO = ensure_coef(WCoef::C4, QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));

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

    complex_t C9_NLO = ensure_coef(WCoef::C9, QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    complex_t C10_NLO = ensure_coef(WCoef::C10, QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    
    complex_t C10_22 = (0.27924 * C1_NLO + 0.33157 * C4_NLO + 2.35917 * Ci_match[8] + 3.29679 * Ci_match[9]) * log(eta) + (1 - eta) * (0.26087 * C9_NLO + 1.15942 * C10_NLO) + C1022;
    
    complex_t C1_NNLO = ensure_coef(WCoef::C1, QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    complex_t C2_NNLO = ensure_coef(WCoef::C1, QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    complex_t C3_NNLO = ensure_coef(WCoef::C1, QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    complex_t C4_NNLO = ensure_coef(WCoef::C1, QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    complex_t C5_NNLO = ensure_coef(WCoef::C1, QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    complex_t C6_NNLO = ensure_coef(WCoef::C1, QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));

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
        ParamId pid {
            ParameterType::WILSON, 
            GroupMapper::str(WGroup::B, ScaleType::HADRONIC, false, BWilsonBasis::STANDARD) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::LO, type)
        };
        LOG_INFO("Storing coefficient", WCoefMapper::flha_full(ids[k], QCDOrder::LO, type));
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, Ci_run[k], 0., 0.));
    }
}

void BCoefficientGroup::base_2_LO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block,
    ContributionType type)
{
    LOG_DEBUG("Init LO running of BCoefficientGroup in traditional basis");

    auto V0 = [src] (int k, int l) -> double {
        return src.at("V_MATRIX")->retrieve(LhaID(0, k, l))->get_val();
    };

    std::array<complex_t, 10> Ci_match_trad = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 8; k++) {
        for (size_t j = 0; j < 8; j++) {
            Ci_match_trad[k] = BRP::std_to_trad_LO[k][j] * ensure_coef(ids[j], QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));;
        }
    }
    Ci_match_trad[8] = ensure_coef(WCoef::C9, QCDOrder::LO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));

    std::array<complex_t, 10> Ci_run {};
    for (size_t k = 0; k < 9; k++) {
        for (size_t l = 0; l < 9; l++) {
            Ci_run[k] += V0(k, l) * Ci_match_trad[l];
        }
        LOG_VERBOSE("C_run_", k + 1, "=", Ci_run[k]);
    }

    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] *= fact;
    Ci_run[9] = ensure_coef(WCoef::C10, QCDOrder::LO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));

    // Store
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {
            ParameterType::WILSON, 
            GroupMapper::str(WGroup::B, ScaleType::HADRONIC, false, BWilsonBasis::TRADITIONAL) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::LO, type)
        };
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, Ci_run[k], 0., 0.));;
    }
}

void BCoefficientGroup::base_1_NLO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block,
    ContributionType type)
{
    LOG_DEBUG("Init NLO running of BCoefficientGroup in standard basis");

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
        Ci_0_match[k] = ensure_coef(ids[k], QCDOrder::LO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
        Ci_1_match[k] = ensure_coef(ids[k], QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
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
        LOG_VERBOSE("C_run_", k + 1, "=", Ci_run[k]);
    }

    // C9 special treatment
    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] = fact * (Ci_run[8] + U0(8, 8) * Ci_0_match[8]);

    // C10
    Ci_run[9] = ensure_coef(WCoef::C10, QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    
    // Store
    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {
            ParameterType::WILSON, 
            GroupMapper::str(WGroup::B, ScaleType::HADRONIC, false, BWilsonBasis::STANDARD) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::NLO, type)
        };
        LOG_INFO("Storing coefficient", WCoefMapper::flha_full(ids[k], QCDOrder::NLO, type));
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, eta * Ci_run[k], 0., 0.));
    }
}

void BCoefficientGroup::base_2_NLO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block,
    ContributionType type)
{
    LOG_DEBUG("Init NLO running of BCoefficientGroup in traditional basis");

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
            complex_t Cj_0_match_std = ensure_coef(ids[j], QCDOrder::LO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
            complex_t Cj_1_match_std = ensure_coef(ids[j], QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
            Ci_0_match_trad[k] = BRP::std_to_trad_LO[k][j] * Cj_0_match_std;
            Ci_1_match_trad[k] = BRP::std_to_trad_LO[k][j] * Cj_1_match_std + BRP::std_to_trad_NLO[k][j] * Cj_0_match_std;
        }
    }
    Ci_0_match_trad[8] = ensure_coef(WCoef::C9, QCDOrder::LO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    Ci_1_match_trad[8] = ensure_coef(WCoef::C9, QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));

    std::array<complex_t, 10> Ci_run {};
    for (size_t k = 0; k < 9; k++) {
        for (size_t l = 0; l < 8; l++) {
            Ci_run[k] += V0(k, l) * Ci_1_match_trad[l] + V1(k, l) * Ci_0_match_trad[l];
        }
        LOG_VERBOSE("C_run_", k + 1, "=", Ci_run[k]);
    }

    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] = fact * (Ci_run[8] + V0(8, 8) * Ci_0_match_trad[8]);
    Ci_run[9] = ensure_coef(WCoef::C10, QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));

    // Store
    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {
            ParameterType::WILSON,
            GroupMapper::str(WGroup::B, ScaleType::HADRONIC, false, BWilsonBasis::TRADITIONAL) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::NLO, type)
        };
        LOG_INFO("Storing coefficient", WCoefMapper::flha_full(ids[k], QCDOrder::NLO, type));
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, eta * Ci_run[k], 0., 0.));
    }
}

void BCoefficientGroup::base_1_NNLO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block,
    ContributionType type)
{
    LOG_DEBUG("Init NNLO running of BCoefficientGroup in standard basis");

    auto U = [src] (int order, int k, int l) -> double {
        return src.at("U_MATRIX")->retrieve(LhaID(order, k, l))->get_val();
    };

    std::array<complex_t, 10> Ci_0_match = {};
    std::array<complex_t, 10> Ci_1_match = {};
    std::array<complex_t, 10> Ci_2_match = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 10; k++) {
        Ci_0_match[k] = ensure_coef(ids[k], QCDOrder::LO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
        Ci_1_match[k] = ensure_coef(ids[k], QCDOrder::NLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
        Ci_2_match[k] = ensure_coef(ids[k], QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
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
        LOG_VERBOSE("C_run_", k + 1, "=", Ci_run[k]);
    }

    // C9 special treatment
    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
    Ci_run[8] = fact * (Ci_run[8] + U(1, 8, 8) * Ci_0_match[8] + U(0, 8, 8) * Ci_1_match[8]);

    // C10
    Ci_run[9] = ensure_coef(WCoef::C10, QCDOrder::NNLO, type, GroupMapper::str(WGroup::B, ScaleType::MATCHING));
    
    // Store
    double eta_sq = pow(src.at("WPARAM_RUN_SM")->retrieve(2)->get_val(), 2);
    for (size_t k = 0; k < 10; k++) {
        ParamId pid {
            ParameterType::WILSON,
            GroupMapper::str(WGroup::B, ScaleType::HADRONIC, false, BWilsonBasis::STANDARD) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::NNLO, type)
        };
        LOG_INFO("Storing coefficient", WCoefMapper::flha_full(ids[k], QCDOrder::NNLO, type));
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

    LOG_VERBOSE("Running matrices updated");
}

std::shared_ptr<CoefficientGroup> BScalarCoefficientGroup::clone() const {
    return std::make_shared<BScalarCoefficientGroup>(*this);
}

BScalarCoefficientGroup::BScalarCoefficientGroup() {
    if (UseMarty().get()) {
        for (auto&& coeff : {"CQ1", "CQ2"}) {
            this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
        }
        return;
    }
    this->insert(std::make_pair("CQ1", std::make_shared<CQ1>()));
    this->insert(std::make_pair("CQ2", std::make_shared<CQ2>()));

    this->id = WGroup::BScalar;
}

void BScalarCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block,
    ContributionType type)
{

    std::array<complex_t, 10> CQi_match = {};
    auto ids = WCoefMapper::get_group(WGroup::BScalar);
    for (size_t k = 0; k < ids.size(); k++) {
        CQi_match[k] = ensure_coef(ids[k], QCDOrder::LO, type, GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING));
    }

    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    double beta_0 = src.at("WPARAM_SI_SM")->retrieve(5)->get_val(); // TODO : change to QCD params
    double fact = pow(eta, -4 / beta_0);
    
    // Store
    for (size_t k = 0; k < ids.size(); k++) {
        ParamId pid {
            ParameterType::WILSON,
            GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::LO, ContributionType::SM)
        };
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, fact * CQi_match[k], 0., 0.));;
    }
}

void BScalarCoefficientGroup::base_1_NLO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block,
    ContributionType type)
{
    auto ids = WCoefMapper::get_group(WGroup::BScalar);
    std::array<complex_t, 2> CQi_match = {};
    for (size_t k = 0; k < ids.size(); k++) {
        CQi_match[k] = ensure_coef(ids[k], QCDOrder::NLO, type, GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING));
    }

    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    double beta_0 = src.at("WPARAM_SI_SM")->retrieve(5)->get_val(); // TODO : change to QCD params
    double fact = pow(eta, 1 - 4 / beta_0);
    
    // Store
    for (size_t k = 0; k < ids.size(); k++) {
        ParamId pid {
            ParameterType::WILSON, 
            GroupMapper::str(WGroup::BScalar, ScaleType::HADRONIC) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::NLO, ContributionType::SM)
        };
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, fact * CQi_match[k], 0., 0.));;
    }
}

void BScalarCoefficientGroup::init_running_blocks(QCDOrder order) {
    LOG_DEBUG("In BScalarCoefficientGroup::init_running_block");

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {GroupMapper::str(WGroup::BScalar, ScaleType::MATCHING), "WPARAM_RUN_SM", "WPARAM_SI_SM"}},
    };

    auto func = [order, this] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        switch (order) {
        case QCDOrder::NNLO:
            LOG_WARN("Scalar B Coefficients are not defined at NNLO.");
        case QCDOrder::NLO:
            BScalarCoefficientGroup::base_1_NLO_calculation(src, dep_block, this->wilson_type);
        case QCDOrder::LO:
            BScalarCoefficientGroup::base_1_LO_calculation(src, dep_block, this->wilson_type);
        }
    };

    WilsonParamComposer().compose_block(GroupMapper::str(this->id, ScaleType::HADRONIC) + "INTER", src, func);

    init_full_running_block(BWilsonBasis::STANDARD);
}

BPrimeCoefficientGroup::BPrimeCoefficientGroup() {
    if (UseMarty().get()) {
        for (auto&& coeff : {"CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"}) {
            this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
        }
        return;
    }

    this->insert(std::make_pair("CP1", std::make_shared<CP1>())); 
    this->insert(std::make_pair("CP2", std::make_shared<CP2>())); 
    this->insert(std::make_pair("CP3", std::make_shared<CP3>()));
    this->insert(std::make_pair("CP4", std::make_shared<CP4>()));  
    this->insert(std::make_pair("CP5", std::make_shared<CP5>())); 
    this->insert(std::make_pair("CP6", std::make_shared<CP6>())); 
    this->insert(std::make_pair("CP7", std::make_shared<CP7>()));  
    this->insert(std::make_pair("CP8", std::make_shared<CP8>()));  
    this->insert(std::make_pair("CP9", std::make_shared<CP9>())); 
    this->insert(std::make_pair("CP10", std::make_shared<CP10>())); 
    this->insert(std::make_pair("CPQ1", std::make_shared<CPQ1>())); 
    this->insert(std::make_pair("CPQ2", std::make_shared<CPQ2>()));

    this->id = WGroup::BPrime;
}

std::shared_ptr<CoefficientGroup> BPrimeCoefficientGroup::clone() const {
    return std::make_shared<BPrimeCoefficientGroup>(*this);
}

void BPrimeCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<std::string, std::shared_ptr<Block>> &src,
    std::shared_ptr<DependentBlock> dep_block,
    ContributionType type)
{
    std::array<complex_t, 12> CPi_match = {};
    auto ids = WCoefMapper::get_group(WGroup::BPrime);
    for (size_t k = 0; k < ids.size(); k++) {
        CPi_match[k] = ensure_coef(ids[k], QCDOrder::LO, type, GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING));
    }

    double eta = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    
    // Store
    for (size_t k = 0; k < ids.size(); k++) {
        ParamId pid {
            ParameterType::WILSON,
            GroupMapper::str(WGroup::BPrime, ScaleType::HADRONIC) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::LO, ContributionType::SM)
        };
        complex_t CPk_run = pow(eta, BRP::exp_prime_running[k]) * CPi_match[k]; 
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, CPk_run, 0., 0.));
    }
}

void BPrimeCoefficientGroup::init_running_blocks(QCDOrder order) {
    LOG_DEBUG("In BPrimeCoefficientGroup::init_running_block");

    if (order > QCDOrder::LO) {
        LOG_WARN("Primed coefficients are only defined at leading order, defaulting to LO.");
    }

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {GroupMapper::str(WGroup::BPrime, ScaleType::MATCHING), "WPARAM_RUN_SM", "WPARAM_SI_SM"}},
    };

    auto func = [this] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        BPrimeCoefficientGroup::base_1_LO_calculation(src, dep_block, this->wilson_type);
    };

    WilsonParamComposer().compose_block(GroupMapper::str(this->id, ScaleType::HADRONIC) + "INTER", src, func);

    init_full_running_block(BWilsonBasis::STANDARD);
}

void CoefficientGroup::claim_coefficients() {
    for (auto& [_, coeff]: *this) {
        coeff->set_owned(true);
        coeff->set_storage_block(GroupMapper::str(this->id, ScaleType::MATCHING));
    }
}

complex_t CoefficientGroup::ensure_coef(WCoef coef, QCDOrder order, ContributionType type, std::string matching_block) {
    ParameterProxy pp {ParameterType::WILSON};
    return pp(matching_block, WCoefMapper::flha_full(coef, order, type));
}

BlnuCoefficientGroup::BlnuCoefficientGroup() {
    if (UseMarty().get()) {
        for (auto&& coeff : {"C_Blnu_A", "C_Blnu_P"}) {
            this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
        }
        return;
    }

    this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A>()));
    this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P>()));

    this->id = WGroup::Blnu;
}

std::shared_ptr<CoefficientGroup> BlnuCoefficientGroup::clone() const {
    return std::make_shared<BlnuCoefficientGroup>(*this);
}

void BlnuCoefficientGroup::init_running_blocks(QCDOrder order) {
    LOG_DEBUG("In BlnuCoefficientGroup::init_running_block");

    if (order > QCDOrder::LO) {
        LOG_WARN("Charged current coefficients are only defined at leading order, defaulting to LO.");
    }

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {GroupMapper::str(WGroup::Blnu, ScaleType::MATCHING)}},
    };

    auto func = [this] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        BlnuCoefficientGroup::base_1_LO_calculation(src, dep_block, this->wilson_type);
    };

    WilsonParamComposer().compose_block(GroupMapper::str(WGroup::Blnu, ScaleType::HADRONIC) + "INTER", src, func);

    init_full_running_block(BWilsonBasis::STANDARD);
}

void BlnuCoefficientGroup::base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>> &src, std::shared_ptr<DependentBlock> dep_block, ContributionType type) {
    std::array<complex_t, 2> C_match = {};
    auto ids = WCoefMapper::get_group(WGroup::Blnu);
    for (size_t k = 0; k < ids.size(); k++) {
        C_match[k] = ensure_coef(ids[k], QCDOrder::LO, type, GroupMapper::str(WGroup::Blnu, ScaleType::MATCHING));
    }
    
    // Store
    for (size_t k = 0; k < ids.size(); k++) {
        ParamId pid {
            ParameterType::WILSON, 
            GroupMapper::str(WGroup::Blnu, ScaleType::HADRONIC) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::LO, ContributionType::SM)
        };
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, C_match[k], 0., 0.));
    }
}

BclnuCoefficientGroup::BclnuCoefficientGroup() {
    if (UseMarty().get()) {
        for (auto&& coeff : {"C_V1", "C_V2", "C_S1", "C_S2", "C_T"}) {
            this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
        }
        return;
    }
    this->insert(std::make_pair("C_V1", std::make_shared<C_V1>()));
    this->insert(std::make_pair("C_V2", std::make_shared<C_V2>()));
    this->insert(std::make_pair("C_S1", std::make_shared<C_S1>()));
    this->insert(std::make_pair("C_S2", std::make_shared<C_S2>()));
    this->insert(std::make_pair("C_T", std::make_shared<C_T>()));

    this->id = WGroup::BCLNU;
}

std::shared_ptr<CoefficientGroup> BclnuCoefficientGroup::clone() const {
    return std::make_shared<BclnuCoefficientGroup>(*this);
}

void BclnuCoefficientGroup::init_running_blocks(QCDOrder order) {
    LOG_TRACE("In BclnuCoefficientGroup::init_running_block");

    if (order > QCDOrder::LO) {
        LOG_WARN("Charged current coefficients are only defined at leading order, defaulting to LO.");
    }

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::WILSON, {GroupMapper::str(WGroup::BCLNU, ScaleType::MATCHING)}},
    };

    auto func = [this] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        BclnuCoefficientGroup::base_1_LO_calculation(src, dep_block, this->wilson_type);
    };

    WilsonParamComposer().compose_block(GroupMapper::str(WGroup::BCLNU, ScaleType::HADRONIC) + "INTER", src, func);

    init_full_running_block(BWilsonBasis::STANDARD);
}

void BclnuCoefficientGroup::base_1_LO_calculation(const std::unordered_map<std::string,std::shared_ptr<Block>>&src, std::shared_ptr<DependentBlock> dep_block, ContributionType type) {
    std::array<complex_t, 5> C_match = {};
    auto ids = WCoefMapper::get_group(WGroup::BCLNU);
    for (size_t k = 0; k < ids.size(); k++) {
        C_match[k] = ensure_coef(ids[k], QCDOrder::LO, type, GroupMapper::str(WGroup::BCLNU, ScaleType::MATCHING));
    }
    
    // Store
    for (size_t k = 0; k < ids.size(); k++) {
        ParamId pid {
            ParameterType::WILSON, 
            GroupMapper::str(WGroup::BCLNU, ScaleType::HADRONIC) + "INTER", 
            WCoefMapper::flha_full(ids[k], QCDOrder::LO, ContributionType::SM)
        };
        dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, C_match[k], 0., 0.));
    }
}

std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs) {
    for(auto& [name, coeff] : coeffs) {
        os << name << " --------------------------------" << std::endl;
        os << "LO at mu_W: "    << coeffs.get_matching_coefficient(name, "LO")      << std::endl;
        os << "LO at mu_h: "    << coeffs.get_running_coefficient(name, "LO")       << std::endl;
        os << "NLO at mu_W: "   << coeffs.get_matching_coefficient(name, "NLO")     << std::endl;
        os << "NLO at mu_h: "   << coeffs.get_running_coefficient(name, "NLO")      << std::endl;
        os << "NNLO at mu_W: "  << coeffs.get_matching_coefficient(name, "NNLO")    << std::endl;
        os << "NNLO at mu_h: "  << coeffs.get_running_coefficient(name, "NNLO")     << std::endl;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs) {
    return os << *coeffs;
} 