#include "Wilson.h"

WilsonCoefficient::WilsonCoefficient(const std::string& name, const std::string& storage_block) : coeffName(name), storage_block(storage_block) {
    if (ends_with(coeffName, "_THDM") || ends_with(coeffName, "_SUSY")) {
        type = ContributionType::BSM;
    }


    for (auto order : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
        matching_info[order] = MatchingInfo(this->id(order, type));
    }
    }

WilsonCoefficient::WilsonCoefficient(const LhaID &name, const std::string& storage_block, ContributionType ct) : coeffName(name.to_string()), storage_block(storage_block) {
    type = ct;


    for (auto order : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
        auto parts = name.get_parts();
        int order_int = order == QCDOrder::LO ? 0 : order == QCDOrder::NLO ? 1 : 2;
        int ct_int = ct == ContributionType::SM ? 0 : ct == ContributionType::BSM ? 1 : 2;
        LhaID full(parts[0], parts[1], order_int, ct_int);
        matching_info[order] = MatchingInfo(full);
    }
}

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
        return;
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
    if (!(*wilson_p).exist(storage_block, this->id(OrderMapper::enum_elt(order), cont_type))) {
            return complex_t(0.,0.);
    }
    return complex_t((*wilson_p)(storage_block, this->id(OrderMapper::enum_elt(order), cont_type)));
}
