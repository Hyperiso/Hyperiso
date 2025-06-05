#include "WilsonGroupSuper.h"

void CoefficientGroup::claim_coefficients() {
    for (auto& [_, coeff]: *this) {
        coeff->set_owned(true);
        coeff->set_storage_block(GroupMapper::str(this->id, ScaleType::MATCHING));
        coeff->set_contribution_type(this->wilson_type);
    }
}

CoefficientGroup::CoefficientGroup(const CoefficientGroup& other)
    : std::map<std::string, std::shared_ptr<WilsonCoefficient>>() 
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
}

CoefficientGroup::CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs) {
    this->insert(coeffs.begin(), coeffs.end());
    
    QCDOrder max_order = QCDOrder::LO;
    for (const auto& [_, wil] : coeffs) {
        if (wil->get_max_order() > max_order) {
            max_order = wil->get_max_order();
        }

        if (max_order == QCDOrder::NNLO) break; 
    }

    this->init(max_order == QCDOrder::NONE ? QCDOrder::LO : max_order);
}

void CoefficientGroup::init(QCDOrder max_order) {
    this->claim_coefficients();
    for (int order = 1; order <= (int)max_order; order++) {
        for (auto& coeff : *this) {
            auto func_wrapper = [&coeff, order](const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src,
                                        std::shared_ptr<DependentParameter> dep_param) {
                auto func = coeff.second->get_func((QCDOrder)order);
                if (!func) {
                    dep_param->set_expected(0.);
                    return;
                }
                // std::cout << "value of : " << coeff.second->get_name() << " at " << OrderMapper::str((QCDOrder)order) << " : " << func(src) << std::endl;
                dep_param->set_expected(func(src));
            };
            WilsonParamComposer().compose_parameter(ParamId{coeff.second->get_storage_block(), coeff.second->get_lhaid((QCDOrder)order)}, coeff.second->get_sources((QCDOrder)order), func_wrapper);
        }
    }
    this->current_order = max_order;
}

complex_t CoefficientGroup::get_matching_coefficient(std::string coeff, std::string order, ContributionType cont_type) const { 
    return this->at(coeff)->get_matching_value(order, cont_type); 
}

complex_t CoefficientGroup::get_running_coefficient(std::string coeff, std::string order, ContributionType cont_type, WilsonBasis basis) const {
    auto coef = this->at(coeff);
    ParameterProxy wilson_p = ParameterProxy(ParameterType::WILSON);
    return complex_t(wilson_p(GroupMapper::str(this->id, ScaleType::HADRONIC, basis), coef->id(OrderMapper::enum_elt(order), cont_type)));
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
