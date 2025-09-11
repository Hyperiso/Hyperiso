#include "ChargedCurrentsWilsonGroup.h"

std::unordered_map<WCoef, scalar_t> BclnuCoefficientGroup::base_1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src
) {
    auto ids = WCoefMapper::get_group(WGroup::BCC);

    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < ids.size(); k++) {
        Ci_run_map[ids[k]] = coef_matching.at(QCDOrder::LO).at(ids[k]);
    }

    return Ci_run_map;
}

// std::unordered_map<WCoef, scalar_t> BlnuCoefficientGroup::base_1_LO_calculation(
//     const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
//     const std::unordered_map<std::string, std::shared_ptr<Block>>& src
// ) {
//     auto ids = WCoefMapper::get_group(WGroup::Blnu);

//     std::unordered_map<WCoef, scalar_t> Ci_run_map {};
//     for (size_t k = 0; k < ids.size(); k++) {
//         Ci_run_map[ids[k]] = coef_matching.at(QCDOrder::LO).at(ids[k]);
//     }

//     return Ci_run_map;
// }

// BlnuCoefficientGroup::BlnuCoefficientGroup(bool force_sm) {
//     this->id = WGroup::Blnu;
//     init_sources();
//     add_wilson_coefficients(force_sm);
// }

// void BlnuCoefficientGroup::init_sources() {
//     std::map<QCDOrder,CoefficientGroupSources> grp_src;
//     grp_src[QCDOrder::LO].sources = {
//         {ParameterType::WILSON, {this->get_matching_storage_block()}},
//     };
//     grp_src[QCDOrder::LO].func = base_1_LO_calculation;
//     this->sources.insert({WilsonBasis::B_STANDARD, grp_src});
// }

// void BlnuCoefficientGroup::add_wilson_coefficients(bool force_sm) {
//     if (UseMarty().get()) {
//         this->wilson_type = force_sm ? ContributionType::SM : ContributionType::TOTAL;
//         for (auto&& coeff : {"C_Blnu_A", "C_Blnu_P"}) {
//             std::string _name = force_sm ? "SM" : MartyModelNameAPI().get();
//             fs::path _path = force_sm ? fs::path(std::string(project_assets_root.data())+"input_files/marty_model/sm.h") : MartyModelPathAPI().get();
//             std::string _block = GroupMapper::str(this->id, ScaleType::MATCHING);
//             LhaID _id = WCoefMapper::flha_full(WCoefMapper::enum_elt(coeff), QCDOrder::LO, this->get_type());
//             this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(_id, _block, _name, _path)));
//         }
//         return;
//     }

//     this->insert(std::make_pair("C_Blnu_A", std::make_shared<CQ1>()));
//     this->insert(std::make_pair("C_Blnu_P", std::make_shared<CQ2>()));
// }

// std::shared_ptr<CoefficientGroup> BlnuCoefficientGroup::clone() const {
//     return std::make_shared<BlnuCoefficientGroup>(*this);
// }

BclnuCoefficientGroup::BclnuCoefficientGroup(bool force_sm) {
    this->id = WGroup::BCC;
    init_sources();
    add_wilson_coefficients(force_sm);
}

void BclnuCoefficientGroup::init_sources() {
    std::map<QCDOrder,CoefficientGroupSources> grp_src;
    grp_src[QCDOrder::LO].sources = {
        {ParameterType::WILSON, {this->get_matching_storage_block()}},
    };
    grp_src[QCDOrder::LO].func = base_1_LO_calculation;
    this->sources.insert({WilsonBasis::B_STANDARD, grp_src});
}

void BclnuCoefficientGroup::add_wilson_coefficients(bool force_sm) {
    if (UseMarty().get()) {
        this->wilson_type = force_sm ? ContributionType::SM : ContributionType::TOTAL;
        for (auto&& coeff : {"C_V1", "C_V2", "C_S1", "C_S2", "C_T"}) {
            std::string _name = force_sm ? "SM" : MartyModelNameAPI().get();
            fs::path _path = force_sm ? fs::path(std::string(project_assets_root.data())+"input_files/marty_model/sm.h") : MartyModelPathAPI().get();
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
}

std::shared_ptr<CoefficientGroup> BclnuCoefficientGroup::clone() const {
    return std::make_shared<BclnuCoefficientGroup>(*this);
}

complex_t CoefficientGroup::ensure_coef(WCoef coef, QCDOrder order, ContributionType type, std::string matching_block) {
    ParameterProxy pp {ParameterType::WILSON};
    return pp(matching_block, WCoefMapper::flha_full(coef, order, type));
}