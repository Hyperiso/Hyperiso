#include "KChargedCurrentWilson.h"

C_V1_su::C_V1_su() : WilsonCoefficient("C_V1_su", GroupMapper::str(WGroup::CC_su, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {},
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_V1_su::compute_LO(const ParamSrc& src) {
    return 1.0;
}