#include "BScalarWilsonGroup.h"

std::shared_ptr<CoefficientGroup> BScalarCoefficientGroup::clone() const {
    return std::make_shared<BScalarCoefficientGroup>(*this);
}

BScalarCoefficientGroup::BScalarCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    LOG_TRACE("In BScalarCoefficientGroup constructor");
    this->id = GroupMapper::to_id(WGroup::BScalar);
}

std::unordered_map<WCoefId, scalar_t> BScalarCoefficientGroup::base_1_LO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    std::array<complex_t, 2> CQi_match = {};
    auto ids = WCoefMapper::get_group(WGroup::BScalar);
    for (size_t k = 0; k < ids.size(); k++) {
        CQi_match[k] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[k]));
    }

    double eta = src.get_val("WPARAM_RUN_SM",2);
    double beta_0 = src.get_val("WPARAM_SI_SM",5); // TODO : change to QCD params
    double fact = pow(eta, -4 / beta_0);
    
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = fact * CQi_match[k];
    }

    return Ci_run_map;
}

std::unordered_map<WCoefId, scalar_t> BScalarCoefficientGroup::base_1_NLO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src
)
{
    auto ids = WCoefMapper::get_group(WGroup::BScalar);
    std::array<complex_t, 2> CQi_match = {};
    for (size_t k = 0; k < ids.size(); k++) {
        CQi_match[k] = coef_matching.at(QCDOrder::NLO).at(WCoefMapper::to_id(ids[k]));
    }

    double eta = src.get_val("WPARAM_RUN_SM",2);
    double beta_0 = src.get_val("WPARAM_SI_SM",5); // TODO : change to QCD params
    double fact = pow(eta, 1 - 4 / beta_0);
    
    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = fact * CQi_match[k];
    }
    
    return Ci_run_map;
}

