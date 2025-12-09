#include "BWilsonGroup.h"

BCoefficientGroup::BCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    LOG_TRACE("In BCoefficientGroup constructor");
    this->id = GroupMapper::to_id(WGroup::B);
}

std::shared_ptr<CoefficientGroup> BCoefficientGroup::clone() const {
    return std::make_shared<BCoefficientGroup>(*this);
}


std::unordered_map<WCoefId, scalar_t> BCoefficientGroup::base_1_LO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching, 
    const BlockSrc& src
)
{
    LOG_DEBUG("Init LO running of BCoefficientGroup in standard basis");

    auto U0 = [src] (int k, int l) -> double {
        return src.get_val("U_MATRIX",LhaID(0, k, l));
    };

    std::array<complex_t, 10> Ci_match = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 10; k++) {
        Ci_match[k] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[k]));
    }

    std::array<complex_t, 10> Ci_run {};

    Ci_match[6] = BRP::C7_eff_std(Ci_match); 
    Ci_match[7] = BRP::C8_eff_std(Ci_match);
    Ci_match[8] *= src.get_val("WPARAM_MATCH_SM", 1)/ (4. * PI);


    for (size_t k = 0; k < 8; k++) {
        for (size_t l = 0; l < 8; l++) {
            Ci_run[k] += U0(k, l) * Ci_match[l];
        }
    }

    for (size_t l = 0; l < 9; l++) {
        Ci_run[8] += 4 * PI / src.get_val("WPARAM_RUN_SM",1) * U0(8, l) * Ci_match[l];
    }

    // C10
    // Ask Nazila : Pourquoi les N et NNLO dans le running de C10 à LO et pas séparé en plusieurs ordres ?
    // Answer : Ask Siavash, need to understand if pb comes from litt or implementation.
    // Answer : Fck around and find out
  

    Ci_run[9] = Ci_match[9];
    
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = Ci_run[k];
    }

    return Ci_run_map;
}

std::unordered_map<WCoefId, scalar_t> BCoefficientGroup::base_2_LO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    LOG_DEBUG("Init LO running of BCoefficientGroup in traditional basis");

    auto V0 = [src] (int k, int l) -> double {
        return src.get_val("V_MATRIX",LhaID(0, k, l));
    };

    std::array<complex_t, 10> Ci_match_trad = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 8; k++) {
        for (size_t j = 0; j < 8; j++) {
            Ci_match_trad[k] += BRP::std_to_trad_LO[k][j] * coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[j]));
            
        }
    }
    Ci_match_trad[8] = src.get_val("WPARAM_MATCH_SM",1)/ (4. * PI) * coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(WCoef::C9));

    std::array<complex_t, 10> Ci_run {};

    for (size_t k = 0; k < 8; k++) {
        for (size_t l = 0; l < 8; l++) {
            Ci_run[k] += V0(k, l) * Ci_match_trad[l];
        }
        LOG_VERBOSE("C_run_", k + 1, "=", Ci_run[k]);
    }

    for (size_t l = 0; l < 9; l++) {
        Ci_run[8] += 4 * PI / src.get_val("WPARAM_RUN_SM",1)* V0(8, l) * Ci_match_trad[l];
    }

    Ci_run[9] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(WCoef::C10));

    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = Ci_run[k];
    }

    return Ci_run_map;
}

std::unordered_map<WCoefId, scalar_t> BCoefficientGroup::base_1_NLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    LOG_DEBUG("Init NLO running of BCoefficientGroup in standard basis");

    auto U0 = [src] (int k, int l) -> double {
        return src.get_val("U_MATRIX",LhaID(0, k, l));
    };

    auto U1 = [src] (int k, int l) -> double {
        return src.get_val("U_MATRIX",LhaID(1, k, l));
    };

    auto C_eff = [] (const std::array<complex_t, 10>& Ci, const double* vec) -> complex_t {
        complex_t C {0};
        for (size_t k = 0; k < 8; k++) {
            C += Ci[k] * vec[k];
        }
        return C;
    };

    double eta = src.get_val("WPARAM_RUN_SM",2);

    std::array<complex_t, 10> Ci_0_match = {};
    std::array<complex_t, 10> Ci_1_match = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 10; k++) {
        Ci_0_match[k] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[k]));
        Ci_1_match[k] = coef_matching.at(QCDOrder::NLO).at(WCoefMapper::to_id(ids[k]));
    }

    std::array<complex_t, 10> Ci_run {};

    Ci_0_match[6] = C_eff(Ci_0_match, BRP::y_std);
    Ci_0_match[7] = C_eff(Ci_0_match, BRP::z_std);
    Ci_1_match[6] = C_eff(Ci_1_match, BRP::y_std); 
    Ci_1_match[7] = C_eff(Ci_1_match, BRP::z_std);
    Ci_0_match[8] *= src.get_val("WPARAM_MATCH_SM",1) / (4 * PI);
    Ci_1_match[8] *= src.get_val("WPARAM_MATCH_SM",1) / (4 * PI);
    
    for (size_t k = 0; k < 8; k++) {
        for (size_t l = 0; l < 8; l++) {
            Ci_run[k] += U0(k, l) * Ci_1_match[l] + U1(k, l) * Ci_0_match[l];
        }
    }
    
    for (size_t l = 0; l < 9; l++) {
        Ci_run[8] += 4 * PI / src.get_val("WPARAM_RUN_SM",1) * (U0(8, l) * Ci_1_match[l] + U1(8, l) * Ci_0_match[l]);
    }

    Ci_run[9] = coef_matching.at(QCDOrder::NLO).at(WCoefMapper::to_id(WCoef::C10));
        
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = eta * Ci_run[k];
    }

    return Ci_run_map;
}

std::unordered_map<WCoefId, scalar_t> BCoefficientGroup::base_2_NLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    LOG_DEBUG("Init NLO running of BCoefficientGroup in traditional basis");

    auto V0 = [src] (int k, int l) -> double {
        return src.get_val("V_MATRIX",LhaID(0, k, l));
    };

    auto V1 = [src] (int k, int l) -> double {
        return src.get_val("V_MATRIX",LhaID(1, k, l));
    };

    std::array<complex_t, 10> Ci_0_match_trad = {};
    std::array<complex_t, 10> Ci_1_match_trad = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 8; k++) {
        for (size_t j = 0; j < 8; j++) {
            complex_t Cj_0_match_std = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[j]));
            complex_t Cj_1_match_std = coef_matching.at(QCDOrder::NLO).at(WCoefMapper::to_id(ids[j]));
            Ci_0_match_trad[k] += BRP::std_to_trad_LO[k][j] * Cj_0_match_std;
            Ci_1_match_trad[k] += BRP::std_to_trad_LO[k][j] * Cj_1_match_std + BRP::std_to_trad_NLO[k][j] * Cj_0_match_std;
        }
    }
    Ci_0_match_trad[8] = src.get_val("WPARAM_MATCH_SM",1) / (4 * PI) * coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(WCoef::C9));
    Ci_1_match_trad[8] = src.get_val("WPARAM_MATCH_SM",1) / (4 * PI) * coef_matching.at(QCDOrder::NLO).at(WCoefMapper::to_id(WCoef::C9));

    std::array<complex_t, 10> Ci_run {};
    for (size_t k = 0; k < 8; k++) {
        for (size_t l = 0; l < 8; l++) {
            Ci_run[k] += V0(k, l) * Ci_1_match_trad[l] + V1(k, l) * Ci_0_match_trad[l];
        }
    }

    for (size_t l = 0; l < 9; l++) {
        Ci_run[8] += 4 * PI / src.get_val("WPARAM_RUN_SM",1) * (V0(8, l) * Ci_1_match_trad[l] + V1(8, l) * Ci_0_match_trad[l]);
    }

    Ci_run[9] = coef_matching.at(QCDOrder::NLO).at(WCoefMapper::to_id(WCoef::C10));

    double eta = src.get_val("WPARAM_RUN_SM",2);
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = eta * Ci_run[k];
    }

    return Ci_run_map;
}

std::unordered_map<WCoefId, scalar_t> BCoefficientGroup::base_1_NNLO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    LOG_DEBUG("Init NNLO running of BCoefficientGroup in standard basis");

    auto U = [src] (int order, int k, int l) -> double {
        return src.get_val("U_MATRIX",LhaID(order, k, l));
    };

    std::array<complex_t, 10> Ci_0_match = {};
    std::array<complex_t, 10> Ci_1_match = {};
    std::array<complex_t, 10> Ci_2_match = {};
    auto ids = WCoefMapper::get_group(WGroup::B);
    for (size_t k = 0; k < 10; k++) {
        Ci_0_match[k] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[k]));
        Ci_1_match[k] = coef_matching.at(QCDOrder::NLO).at(WCoefMapper::to_id(ids[k]));
        Ci_2_match[k] = coef_matching.at(QCDOrder::NNLO).at(WCoefMapper::to_id(ids[k]));
    }

    Ci_0_match[6] = BRP::C7_eff_std(Ci_0_match);
    Ci_0_match[7] = BRP::C8_eff_std(Ci_0_match);
    Ci_1_match[6] = BRP::C7_eff_std(Ci_1_match);
    Ci_1_match[7] = BRP::C8_eff_std(Ci_1_match);
    Ci_2_match[6] = BRP::C7_eff_std(Ci_2_match);
    Ci_2_match[7] = BRP::C8_eff_std(Ci_2_match);
    Ci_0_match[8] *= src.get_val("WPARAM_MATCH_SM",1) / (4 * PI);
    Ci_1_match[8] *= src.get_val("WPARAM_MATCH_SM",1) / (4 * PI);
    Ci_2_match[8] *= src.get_val("WPARAM_MATCH_SM",1) / (4 * PI);

    std::array<complex_t, 10> Ci_run {};

    for (size_t k = 0; k < 8; k++) {
        for (size_t l = 0; l < 8; l++) {
            Ci_run[k] += U(2, k, l) * Ci_0_match[l] + U(1, k, l) * Ci_1_match[l] + U(0, k, l) * Ci_2_match[l];
        }
    }

    for (size_t l = 0; l < 9; l++) {
        Ci_run[8] += 4 * PI / src.get_val("WPARAM_RUN_SM",1) * (U(2, 8, l) * Ci_0_match[l] + U(1, 8, l) * Ci_1_match[l] + U(0, 8, l) * Ci_2_match[l]);
    }

    Ci_run[9] = coef_matching.at(QCDOrder::NNLO).at(WCoefMapper::to_id(WCoef::C10));
    
    double eta_sq = pow(src.get_val("WPARAM_RUN_SM",2), 2);
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = eta_sq * Ci_run[k];
    }

    return Ci_run_map;
}

