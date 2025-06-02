#include "ChargedCurrentsWilsonGroup.h"

std::unordered_map<WCoef, scalar_t> BclnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src
) {
    std::array<complex_t, 5> C_match = {};
    auto ids = WCoefMapper::get_group(WGroup::BCLNU);
    for (size_t k = 0; k < ids.size(); k++) {
        C_match[k] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }
    
    // Store
    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[ids[k]] = C_match[k];
    }

    return Ci_run_map;
}

std::unordered_map<WCoef, scalar_t> BlnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src
) {
std::array<complex_t, 2> C_match = {};
    auto ids = WCoefMapper::get_group(WGroup::Blnu);
    for (size_t k = 0; k < ids.size(); k++) {
        C_match[k] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }
    
    // Store
    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 10; k++) {
        Ci_run_map[ids[k]] = C_match[k];
    }

    return Ci_run_map;
}

BlnuCoefficientGroup::BlnuCoefficientGroup() {
    if (UseMarty().get()) {
        for (auto&& coeff : {"C_Blnu_A", "C_Blnu_P"}) {
            std::string _name = MartyModelNameAPI().get();
            fs::path _path = MartyModelPathAPI().get();
            std::string _block = GroupMapper::str(this->id, ScaleType::MATCHING);
            LhaID _id = WCoefMapper::flha_full(WCoefMapper::enum_elt(coeff), QCDOrder::LO, this->get_type());
            this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(_id, _block, _name, _path)));
        }
        return;
    }

    this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A>()));
    this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P>()));

    this->id = WGroup::Blnu;
}

std::shared_ptr<CoefficientGroup> BlnuCoefficientGroup::clone() const {
    return std::make_shared<BlnuCoefficientGroup>(*this);
}

BclnuCoefficientGroup::BclnuCoefficientGroup() {
    if (UseMarty().get()) {
        for (auto&& coeff : {"C_V1", "C_V2", "C_S1", "C_S2", "C_T"}) {
            std::string _name = MartyModelNameAPI().get();
            fs::path _path = MartyModelPathAPI().get();
            std::string _block = GroupMapper::str(this->id, ScaleType::MATCHING);
            LhaID _id = WCoefMapper::flha_full(WCoefMapper::enum_elt(coeff), QCDOrder::LO, this->get_type());
            this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(_id, _block, _name, _path)));
        }
        return;
    }
    this->insert(std::make_pair("C_V1", std::make_shared<C_V1>()));
    this->insert(std::make_pair("C_V2", std::make_shared<C_V2>()));
    this->insert(std::make_pair("C_S1", std::make_shared<C_S1>()));
    this->insert(std::make_pair("C_S2", std::make_shared<C_S2>()));
    this->insert(std::make_pair("C_T", std::make_shared<C_T>()));

    this->id = WGroup::BCLNU;
}

std::shared_ptr<CoefficientGroup> BclnuCoefficientGroup::clone() const {
    return std::make_shared<BclnuCoefficientGroup>(*this);
}

complex_t CoefficientGroup::ensure_coef(WCoef coef, QCDOrder order, ContributionType type, std::string matching_block) {
    ParameterProxy pp {ParameterType::WILSON};
    return pp(matching_block, WCoefMapper::flha_full(coef, order, type));
}