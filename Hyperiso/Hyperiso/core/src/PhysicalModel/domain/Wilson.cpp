#include "Wilson.h"

void WilsonCoefficient::init(QCDOrder order) {
    if (!is_owned) {
        WilsonParameterHelper::init(2);
    }

    max_order = order;
    std::cout << OrderMapper::str(order) << " ;;;" << std::endl;
    switch (order) {
    case QCDOrder::NNLO:
        NNLO_calculation();
    case QCDOrder::NLO:
        NLO_calculation();
    case QCDOrder::LO:
        LO_calculation();
        break;
    default:
        LOG_ERROR("logicerror", "QCDOrder cannot be none.");
        break;
    }
}

void WilsonCoefficient::set_owned(bool owned) {
    if (this->is_owned && owned) {
        LOG_ERROR("LogicError", "WilsonCoefficient is already owned by a WilsonGroup and cannot be shared.");
    }

    this->is_owned = owned;
}

LhaID WilsonCoefficient::id(QCDOrder order) const {
    return WCoefMapper::flha_full(WCoefMapper::enum_elt(this->coeffName), order, type);
}

bool WilsonCoefficient::operator==(const WilsonCoefficient &other) const {
    return this->coeffName == other.coeffName
            && this->type == other.type
            && this->is_owned == other.is_owned
            && this->from_lha == other.from_lha;
}

complex_t WilsonCoefficient::get_matching_value(std::string order) const {
    ParameterProxy wilson_p {ParameterType::WILSON};
    return complex_t(wilson_p(storage_block, this->id(OrderMapper::enum_elt(order))));
}
