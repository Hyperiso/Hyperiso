#include "BScalarWilsonGroup.h"

std::shared_ptr<CoefficientGroup> BScalarCoefficientGroup::clone() const {
    return std::make_shared<BScalarCoefficientGroup>(*this);
}

BScalarCoefficientGroup::BScalarCoefficientGroup(WilsonGroupAdapterConfig adapters) : CoefficientGroup(adapters) {
    LOG_TRACE("In BScalarCoefficientGroup constructor");
    this->id = GroupMapper::to_id(WGroup::BScalar);
}

std::unordered_map<WCoefId, scalar_t> BScalarCoefficientGroup::base_1_LO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    auto ids = WCoefMapper::get_group(WGroup::BScalar);
    double eta = src.get_val("WPARAM_RUN_SM",2);
    double beta_0 = src.get_val("WPARAM_SI_SM",5);
    double fact = pow(eta, -4 / beta_0);
    
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (auto coef : ids) {
        const auto id = WCoefMapper::to_id(coef);
        Ci_run_map[id] = fact * coef_matching.at(QCDOrder::LO).at(id);
    }

    return Ci_run_map;
}

std::unordered_map<WCoefId, scalar_t> BScalarCoefficientGroup::base_1_NLO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    auto ids = WCoefMapper::get_group(WGroup::BScalar);
    double eta = src.get_val("WPARAM_RUN_SM",2);
    double beta_0 = src.get_val("WPARAM_SI_SM",5);
    double fact = pow(eta, 1 - 4 / beta_0);
    
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (auto coef : ids) {
        const auto id = WCoefMapper::to_id(coef);
        Ci_run_map[id] = fact * coef_matching.at(QCDOrder::NLO).at(id);
    }
    
    return Ci_run_map;
}
