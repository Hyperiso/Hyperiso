#include "Wilson.h"


std::function<scalar_t(const ParamSrc&)> WilsonCoefficient::get_func(QCDOrder order) {
    return this->matching_info[order].compute;
}


std::unordered_set<ParamId> WilsonCoefficient::get_sources(QCDOrder order) {
    return this->matching_info[order].sources;
}

LhaID WilsonCoefficient::get_lhaid(QCDOrder order) {
    return this->matching_info[order].lhaid;
}

LhaID WilsonCoefficient::get_lhaid_from_name(QCDOrder order) {
    std::pair<int,int> base_id = WCoefMapper::flha_base(WCoefMapper::id_of(this->get_base_name()));
    int cont = type == ContributionType::SM ? 0 : type == ContributionType::BSM ? 1 : 2;
    int ord = order == QCDOrder::LO ? 0 : order == QCDOrder::NLO ? 1 : 2;
    return LhaID(base_id.first, base_id.second, ord, cont);
}

std::string WilsonCoefficient::get_base_name() const {
    std::string name = this->coeffName;

    if (ends_with(name, "_THDM")) {
        name = name.substr(0, name.size() - 5);
    } else if (ends_with(name, "_SUSY")) {
        name = name.substr(0, name.size() - 5);
    }

    return name;
}

void WilsonCoefficient::set_owned(bool owned)
{
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

LhaID WilsonCoefficient::id(QCDOrder order, ContributionType typ) const {
    return WCoefMapper::flha_full(WCoefMapper::enum_elt(get_base_name()), order, typ);
}

bool WilsonCoefficient::operator==(const WilsonCoefficient &other) const {
    return this->coeffName == other.coeffName
            && this->type == other.type
            && this->is_owned == other.is_owned;
}

complex_t WilsonCoefficient::get_matching_value(std::string order, ContributionType cont_type, std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_p) const {
    return complex_t((*wilson_p)(storage_block, this->id(OrderMapper::enum_elt(order), cont_type)));
}
