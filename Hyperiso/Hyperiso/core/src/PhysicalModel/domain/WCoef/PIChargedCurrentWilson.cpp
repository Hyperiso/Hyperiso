#include "PIChargedCurrentWilson.h"

C_V1_du::C_V1_du() : WilsonCoefficient("C_V1_du", GroupMapper::str(WGroup::CC_du, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_V1_du::compute_LO(const ParamSrc&) {
    return 1.0;
}