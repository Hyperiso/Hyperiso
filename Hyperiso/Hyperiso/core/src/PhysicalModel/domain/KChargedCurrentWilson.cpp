#include "KChargedCurrentWilson.h"

C_V1_su::C_V1_su() : WilsonCoefficient("C_V1_su", GroupMapper::str(WGroup::BCC_su, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        WCoefMapper::flha_full(WCoef::C_V1_su, QCDOrder::LO, this->type)
    };
}

double C_V1_su::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&) {
    return 1.0;
}