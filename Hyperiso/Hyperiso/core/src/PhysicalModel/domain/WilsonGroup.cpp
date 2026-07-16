#include "WilsonGroup.h"

#include <algorithm>

void CoefficientGroup::claim_coefficients() {
    for (auto& [_, coeff]: *this) {
        coeff->set_owned(true);
        coeff->set_storage_block(this->get_matching_storage_block());
        coeff->set_contribution_type(this->wilson_type);
    }
}

CoefficientGroup::CoefficientGroup(const CoefficientGroup& other)
    : std::map<std::string, std::shared_ptr<WilsonCoefficient>> (), adapters(other.adapters) 
{
    for (const auto& [key, ptr] : other) {
        if (ptr) {
            (*this)[key] = ptr->clone();
        } else {
            (*this)[key] = nullptr;
        }
    }
    wilson_type = other.wilson_type;
    current_order = other.current_order;
    id = other.id;
    member_ids = other.member_ids;
    sources = other.sources;
    block_name = other.block_name;
}

CoefficientGroup::CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs, WilsonGroupAdapterConfig adapters) : CoefficientGroup(adapters) {
    this->insert(coeffs.begin(), coeffs.end());
    
    QCDOrder max_order = QCDOrder::NNLO;

    this->init(max_order == QCDOrder::NONE ? QCDOrder::LO : max_order);
}

void CoefficientGroup::init(QCDOrder max_order) {
    this->claim_coefficients();
    for (int order = 1; order <= (int)max_order; order++) {
        for (auto& coeff : *this) {
            auto func_wrapper = [&coeff, order](const ParamSrc& src,
                                        std::shared_ptr<DependentParameter> dep_param) {
                auto func = coeff.second->get_func((QCDOrder)order);
                if (!func) {
                    dep_param->set_expected(0.);
                    return;
                }
                dep_param->set_expected(func(src));
            };
            adapters.iblock_c->compose_parameter(ParamId{coeff.second->get_storage_block(), coeff.second->get_lhaid((QCDOrder)order)}, coeff.second->get_sources((QCDOrder)order), func_wrapper);
        }
    }
    this->current_order = max_order;
}

void CoefficientGroup::set_matching_storage_block(std::string name) {
    block_name = name;
}

void CoefficientGroup::add_member_id(WCoefId id) {
    if (std::find(member_ids.begin(), member_ids.end(), id) == member_ids.end()) {
        member_ids.emplace_back(std::move(id));
    }
}
complex_t CoefficientGroup::get_matching_coefficient(std::string coeff, std::string order, ContributionType cont_type) const {
    auto it = this->find(coeff);
    if (it == this->end()) {
        throw std::out_of_range("Coefficient '" + coeff + "' not found in group '" + GroupMapper::str(this->id) + "'");
    }
    return it->second->get_matching_value(order, cont_type, adapters.wilson_proxy);
}

complex_t CoefficientGroup::get_running_coefficient(std::string coeff, std::string order, ContributionType cont_type, WilsonBasis basis) const {
    auto it = this->find(coeff);
    if (it == this->end()) {
        throw std::out_of_range("Coefficient '" + coeff + "' not found in group '" + GroupMapper::str(this->id) + "'");
    }
    auto& c = it->second;

    if (!(*adapters.wilson_proxy).exist(GroupMapper::str(this->id, ScaleType::HADRONIC, basis),
        c->id(OrderMapper::enum_elt(order), cont_type))) {
            return complex_t(0.,0.);
    }

    return complex_t((*adapters.wilson_proxy)(
        GroupMapper::str(this->id, ScaleType::HADRONIC, basis),
        c->id(OrderMapper::enum_elt(order), cont_type)
    ));
}

QCDOrder CoefficientGroup::get_order(){
    return this->current_order;
}

std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs) {
    for(auto& [name, coeff] : coeffs) {
        os << name << " --------------------------------" << std::endl;
        os << "LO at mu_W: "    << coeffs.get_matching_coefficient(name, "LO", ContributionType::TOTAL)      << std::endl;
        os << "LO at mu_h: "    << coeffs.get_running_coefficient(name, "LO", ContributionType::TOTAL)       << std::endl;
        os << "NLO at mu_W: "   << coeffs.get_matching_coefficient(name, "NLO", ContributionType::TOTAL)     << std::endl;
        os << "NLO at mu_h: "   << coeffs.get_running_coefficient(name, "NLO", ContributionType::TOTAL)      << std::endl;
        os << "NNLO at mu_W: "  << coeffs.get_matching_coefficient(name, "NNLO", ContributionType::TOTAL)    << std::endl;
        os << "NNLO at mu_h: "  << coeffs.get_running_coefficient(name, "NNLO", ContributionType::TOTAL)     << std::endl;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs) {
    return os << *coeffs;
} 
