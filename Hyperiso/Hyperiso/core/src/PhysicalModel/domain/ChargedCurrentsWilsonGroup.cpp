#include "ChargedCurrentsWilsonGroup.h"

std::unordered_map<WCoef, scalar_t> BclnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src
) {
    std::array<complex_t, 5> C_match = {};
    auto ids = WCoefMapper::get_group(WGroup::BCLNU);
    for (size_t k = 0; k < ids.size(); k++) {
        C_match[k] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }
    
    // Store
    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[ids[k]] = C_match[k];
    }

    return Ci_run_map;
}

std::unordered_map<WCoef, scalar_t> BlnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src
) {
std::array<complex_t, 2> C_match = {};
    auto ids = WCoefMapper::get_group(WGroup::Blnu);
    for (size_t k = 0; k < ids.size(); k++) {
        C_match[k] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }
    
    // Store
    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[ids[k]] = C_match[k];
    }

    return Ci_run_map;
}