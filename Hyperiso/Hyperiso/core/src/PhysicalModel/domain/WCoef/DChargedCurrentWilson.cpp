#include "DChargedCurrentWilson.h"

C_V1_cs::C_V1_cs() : WilsonCoefficient("C_V1_cs", GroupMapper::str(WGroup::CC_cs, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_V1_cs::compute_LO(const ParamSrc&) {
    return 1.0;
}



C_V1_cd::C_V1_cd() : WilsonCoefficient("C_V1_cd", GroupMapper::str(WGroup::CC_cd, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_V1_cd::compute_LO(const ParamSrc&) {
    return 1.0;
}