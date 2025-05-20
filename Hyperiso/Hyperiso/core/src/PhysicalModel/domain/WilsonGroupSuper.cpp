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
    basis = other.basis;
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
                std::cout << "value of : " << coeff.second->get_name() << " at " << OrderMapper::str((QCDOrder)order) << " : " << func(src) << std::endl;;
                dep_param->set_expected(func(src));
            };
            WilsonParamComposer().compose_parameter(ParamId{coeff.second->get_storage_block(), coeff.second->get_lhaid((QCDOrder)order)}, coeff.second->get_sources((QCDOrder)order), func_wrapper);
        }
    }
    this->current_order = max_order;
}

complex_t CoefficientGroup::get_matching_coefficient(std::string coeff, std::string order) const { 
    return this->at(coeff)->get_matching_value(order); 
}

complex_t CoefficientGroup::get_running_coefficient(std::string coeff, std::string order) const {
    auto coef = this->at(coeff);
    ParameterProxy wilson_p = ParameterProxy(ParameterType::WILSON);
    // std::cout << "before : " << coef->id(OrderMapper::enum_elt(order)) << std::endl;
    // std::cout << GroupMapper::str(this->id, ScaleType::HADRONIC, false, this->basis.value_or(BWilsonBasis::STANDARD)) << " eheh " << coef->id(OrderMapper::enum_elt(order)) << std::endl;
    return complex_t(wilson_p(GroupMapper::str(this->id, ScaleType::HADRONIC, false, this->basis.value_or(BWilsonBasis::STANDARD)), coef->id(OrderMapper::enum_elt(order))));
}

QCDOrder CoefficientGroup::get_order(){
    return this->current_order;
}

bool CoefficientGroup::is_double_basis() const {
    return this->basis.has_value();
}

void CoefficientGroup::init_full_running_block(const std::unordered_map<ParameterType, std::vector<std::string>> &source_names, BWilsonBasis basis, bool inter, std::vector<ContributionType> contribution_type) {
    // std::unordered_map<ParameterType, std::vector<std::string>> src = {
    //     {ParameterType::WILSON, {GroupMapper::str(this->id, ScaleType::HADRONIC) + "INTER", "WPARAM_RUN_SM"}}
    // };

    auto func = [basis, this, inter, contribution_type] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double alpha_s_mu_h = src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();
        double fact = alpha_s_mu_h / 4 * PI;
        for (WCoef coef_id : WCoefMapper::get_group(this->id)) {
            std::map<ContributionType, complex_t> coeff_fulls;
            for (ContributionType type : contribution_type) {
                for (int order=0; order < 2; order++) {
                    coeff_fulls[type] += ensure_coef(coef_id, (QCDOrder)(order + 1), type, GroupMapper::str(this->id, ScaleType::HADRONIC, false, basis) + (inter ? "INTER" : ""));
                }

            }

            auto coef_lha_base = WCoefMapper::flha_base(coef_id);
            // LOG_INFO("Storing full coefficient", LhaID(coef_lha_base.first, coef_lha_base.second, (int)this->wilson_type));
            for (ContributionType type : contribution_type) {
                ParamId pid {
                    ParameterType::WILSON, 
                    GroupMapper::str(this->id, ScaleType::HADRONIC, true, basis) + (inter ? "INTER" : ""), 
                    LhaID(coef_lha_base.first, coef_lha_base.second, (int)type)
                };
                dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, coeff_fulls[type], 0., (int)type));
            }
        }
    };
    WilsonParamComposer().compose_block(GroupMapper::str(this->id, ScaleType::HADRONIC, true, basis) + (inter ? "INTER" : ""), source_names, func);
}

void CoefficientGroup::switch_basis() {
    if (!basis.has_value()) return;

    basis = (basis.value() == BWilsonBasis::STANDARD ? BWilsonBasis::TRADITIONAL : BWilsonBasis::STANDARD);
}

std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs) {
    for(auto& [name, coeff] : coeffs) {
        os << name << " --------------------------------" << std::endl;
        os << "LO at mu_W: "    << coeffs.get_matching_coefficient(name, "LO")      << std::endl;
        os << "LO at mu_h: "    << coeffs.get_running_coefficient(name, "LO")       << std::endl;
        os << "NLO at mu_W: "   << coeffs.get_matching_coefficient(name, "NLO")     << std::endl;
        os << "NLO at mu_h: "   << coeffs.get_running_coefficient(name, "NLO")      << std::endl;
        os << "NNLO at mu_W: "  << coeffs.get_matching_coefficient(name, "NNLO")    << std::endl;
        os << "NNLO at mu_h: "  << coeffs.get_running_coefficient(name, "NNLO")     << std::endl;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs) {
    return os << *coeffs;
} 
