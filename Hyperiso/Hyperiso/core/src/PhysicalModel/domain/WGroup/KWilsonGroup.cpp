#include "KWilsonGroup.h"

std::unordered_map<WCoef, scalar_t> KCoefficientGroup::base_1_LO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, 
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src
)
{

    auto ids = WCoefMapper::get_group(WGroup::B);

    std::array<complex_t, 10> Ci_run {};
    
    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[ids[k]] = Ci_run[k];
        // LOG_DEBUG("At hadronic scale:", k+1, "=", Ci_run_map[ids[k]]);
    }

    return Ci_run_map;
}