#include "PIChargedCurrentWilson.h"

C_V1_du::C_V1_du() : WilsonCoefficient("C_V1_du", GroupMapper::str(WGroup::CC_du, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        WCoefMapper::flha_full(WCoef::C_V1_du, QCDOrder::LO, this->type)
    };
}

double C_V1_du::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&) {
    return 1.0;
}