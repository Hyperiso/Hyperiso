#include "KWilsonGroup.h"

std::unordered_map<WCoefId, scalar_t> KCoefficientGroup::base_1_LO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&, 
    const BlockSrc&
)
{

    auto ids = WCoefMapper::get_group(WGroup::K);

    std::array<complex_t, 10> Ci_run {};
    
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = Ci_run[k];
    }

    return Ci_run_map;
}