#include "WilsonSuper.h"


std::function<double(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&)> WilsonCoefficient::get_func(QCDOrder order) {
    return this->matching_info[order].compute;
}


std::unordered_set<ParamId> WilsonCoefficient::get_sources(QCDOrder order) {
    return this->matching_info[order].sources;
}

LhaID WilsonCoefficient::get_lhaid(QCDOrder order) {
    return this->matching_info[order].lhaid;
}
void WilsonCoefficient::set_owned(bool owned) {
    if (this->is_owned && owned) {
        LOG_ERROR("LogicError", "WilsonCoefficient is already owned by a WilsonGroup and cannot be shared.");
    }

    this->is_owned = owned;
}

void WilsonCoefficient::set_storage_block(std::string block_name) {
    this->storage_block = block_name;
}

void WilsonCoefficient::set_contribution_type(ContributionType type) {
    this->type = type;
}

//TODO : disgusting
LhaID WilsonCoefficient::id(QCDOrder order) const {
    std::string name = this->coeffName;
    size_t pos = name.rfind('_');
    std::string trimmedName = (pos != std::string::npos) ? name.substr(0, pos) : name;
    return WCoefMapper::flha_full(WCoefMapper::enum_elt(trimmedName), order, type);
}

bool WilsonCoefficient::operator==(const WilsonCoefficient &other) const {
    return this->coeffName == other.coeffName
            && this->type == other.type
            && this->is_owned == other.is_owned;
}

complex_t WilsonCoefficient::get_matching_value(std::string order) const {
    ParameterProxy wilson_p {ParameterType::WILSON};
    return complex_t(wilson_p(storage_block, this->id(OrderMapper::enum_elt(order))));
}
