#include "ChargedCurrentWilson.h"

C_V1_bc::C_V1_bc() : WilsonCoefficient("C_V1_bc", GroupMapper::str(WGroup::CC_bc, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_V1_bc::compute_LO(const ParamSrc&) {
    return 1.0;
}

C_V1_bu::C_V1_bu() : WilsonCoefficient("C_V1_bu", GroupMapper::str(WGroup::CC_bu, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_V1_bu::compute_LO(const ParamSrc&) {
    return 1.0;
}