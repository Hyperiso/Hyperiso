#include "DChargedCurrentWilson.h"

C_V1_cs::C_V1_cs() : WilsonCoefficient("C_V1_cs", GroupMapper::str(WGroup::CC_cs, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        WCoefMapper::flha_full(WCoef::C_V1_cs, QCDOrder::LO, this->type)
    };
}

double C_V1_cs::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&) {
    return 1.0;
}



C_V1_cd::C_V1_cd() : WilsonCoefficient("C_V1_cd", GroupMapper::str(WGroup::CC_cd, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        WCoefMapper::flha_full(WCoef::C_V1_cd, QCDOrder::LO, this->type)
    };
}

double C_V1_cd::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&) {
    return 1.0;
}