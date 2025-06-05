#include "ChargedCurrentWilsonSuper.h"

C_V1::C_V1() : WilsonCoefficient("C_V1", GroupMapper::str(WGroup::BCC, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        WCoefMapper::flha_full(WCoef::C_V1, QCDOrder::LO, this->type)
    };
}

double C_V1::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&) {
    return 1.0;
}