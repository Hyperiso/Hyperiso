#include "Wilson.h"

void WilsonCoefficient::init(QCDOrder order) {
    if (!is_owned) {
        WilsonParameterHelper::init(2);
    }

    switch (order) {
    case QCDOrder::NNLO:
        NNLO_calculation();
        is_calculated["NNLO"] = true;
    case QCDOrder::NLO:
        NLO_calculation();
        is_calculated["NLO"] = true;
    case QCDOrder::LO:
        LO_calculation();
        is_calculated["LO"] = true;
    default:
        break;
    }
}

void WilsonCoefficient::set_owned(bool owned) {
    if (this->is_owned && owned) {
        LOG_ERROR("LogicError", "WilsonCoefficient is already owned by a WilsonGroup and cannot be shared.");
    }

    this->is_owned = owned;
}

bool WilsonCoefficient::operator==(const WilsonCoefficient& other) const {
    return this->coeffName == other.coeffName
            && this->type == other.type
            && this->is_owned == other.is_owned
            && this->from_lha == other.from_lha;
}

complex_t WilsonCoefficient::get_CoefficientMatchingValue(std::string order) const {
    LhaID code(WCoefMapper::flha_full(WCoefMapper::enum_elt(this->coeffName), OrderMapper::enum_elt(order), type));
    ParameterProvider wilson_p = ParameterProvider(ParameterType::WILSON);
    return complex_t(wilson_p("B_MATCH", code));
}
