// #ifndef WCOEFF_MAPPER_H
// #define WCOEFF_MAPPER_H

// #include "GeneralEnum.h"
// #include "EnumMapper.h"

// #include "General.h"

// class WCoefMapper : public EnumMapperBase<WCoef, WCoefMapper> {
// public:
//     static const std::map<WCoef, std::string>& mapping() {
//         return _mapping();
//     }

//     static const std::map<std::string, WCoef>& inverse_mapping() {
//         static const std::map<std::string, WCoef> inv = invert_map(mapping());
//         return inv;
//     }

//     static std::pair<int, int> flha_base(WCoef coef) {
//         return flha_mapping().at(coef);
//     }

//     static LhaID flha_full(WCoef coef, QCDOrder order, ContributionType type) {
//         auto [a, b] = flha_mapping().at(coef);
//         return LhaID{a, b, static_cast<int>(order) - 1, static_cast<int>(type)};
//     }

//     static WCoef from_flha(int content, int structure) {
//         auto key = std::make_pair(content, structure);
//         if (!inverse_flha_mapping().contains(key)) {
//             LOG_ERROR("General", "Wilson coefficient with ID", content, structure, "is not supported");
//         }
//         return inverse_flha_mapping().at(key);
//     }

//     static std::vector<WCoef> get_group(WGroup group) {
//         switch (group) {
//             case WGroup::B: return B_group();
//             case WGroup::BPrime: return B_prime_group();
//             case WGroup::BScalar: return B_scalar_group();
//             case WGroup::CC_bc: return b_clnu_group();
//             case WGroup::MESON_MIXING: return meson_mixing_group();
//             default:
//                 LOG_ERROR("Invalid WGroup", "get_group function couldn't find your group");
//         }
//     }

//     static size_t n_wilsons() {
//         return mapping().size();
//     }

//     static const std::map<WCoef, std::pair<int, int>>& flha_mapping() {
//         return _flha_mapping();
//     }

//     static const std::map<std::pair<int, int>, WCoef>& inverse_flha_mapping() {
//         static const std::map<std::pair<int, int>, WCoef> inv = invert_map(flha_mapping());
//         return inv;
//     }

//     static const std::vector<WCoef>& B_group() {
//         static const std::vector<WCoef> g = {WCoef::C1, WCoef::C2, WCoef::C3, WCoef::C4, WCoef::C5, WCoef::C6, WCoef::C7, WCoef::C8, WCoef::C9, WCoef::C10};
//         return g;
//     }

//     static const std::vector<WCoef>& B_prime_group() {
//         static const std::vector<WCoef> g = {WCoef::CP1, WCoef::CP2, WCoef::CP3, WCoef::CP4, WCoef::CP5, WCoef::CP6, WCoef::CP7, WCoef::CP8, WCoef::CP9, WCoef::CP10, WCoef::CPQ1, WCoef::CPQ2};
//         return g;
//     }

//     static const std::vector<WCoef>& B_scalar_group() {
//         static const std::vector<WCoef> g = {WCoef::CQ1, WCoef::CQ2};
//         return g;
//     }

//     static const std::vector<WCoef>& b_clnu_group() {
//         static const std::vector<WCoef> g = {WCoef::C_V1, WCoef::C_V2, WCoef::C_S1, WCoef::C_S2, WCoef::C_T};
//         return g;
//     }

//     // Don't change the order in the vectors !
//     static const std::vector<WCoef>& meson_mixing_group() {
//         static const std::vector<WCoef> g = {
//             WCoef::C_BD_1, WCoef::C_BD_2, WCoef::C_BD_3, WCoef::C_BD_4, WCoef::C_BD_5, WCoef::CT_BD_1, WCoef::CT_BD_2, WCoef::CT_BD_3,
//             WCoef::C_BS_1, WCoef::C_BS_2, WCoef::C_BS_3, WCoef::C_BS_4, WCoef::C_BS_5, WCoef::CT_BS_1, WCoef::CT_BS_2, WCoef::CT_BS_3,
//             WCoef::C_SD_1, WCoef::C_SD_2, WCoef::C_SD_3, WCoef::C_SD_4, WCoef::C_SD_5, WCoef::CT_SD_1, WCoef::CT_SD_2, WCoef::CT_SD_3,
//             WCoef::C_CU_1, WCoef::C_CU_2, WCoef::C_CU_3, WCoef::C_CU_4, WCoef::C_CU_5, WCoef::CT_CU_1, WCoef::CT_CU_2, WCoef::CT_CU_3,
//         };
//         return g;
//     }

// private:
//     static const std::map<WCoef, std::string>& _mapping();
//     static const std::map<WCoef, std::pair<int, int>>& _flha_mapping();
// };

// #endif


#pragma once
#include "wcoef_ids.hpp"
