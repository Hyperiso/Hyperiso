#include "MesonMixingWilson.h"

C_mix_bd_1::C_mix_bd_1() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {ParameterType::SM, "SMINPUTS", 2},                   // G_F
            {ParameterType::SM, "MASS", 24},                      // M_W
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

complex_t C_mix_bd_1::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    double G_F = src.at({ParameterType::SM, "SMINPUTS", 2})->get_val();
    double M_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
    complex_t V_td = src.at({ParameterType::SM, "VCKM", {2, 0}})->get_val();
    complex_t V_tb = src.at({ParameterType::SM, "VCKM", {2, 2}})->get_val();
    return pow(G_F * M_W * std::conj(V_td) * V_tb / (2 * PI), 2) * S0(xt);
}

C_mix_bd_1_tilde::C_mix_bd_1_tilde(): WilsonCoefficient("CT_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1050105, 4242, 0, 0));
}

C_mix_bd_2::C_mix_bd_2(): WilsonCoefficient("C_BD_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1050105, 3131, 0, 0));
}

C_mix_bd_2_tilde::C_mix_bd_2_tilde(): WilsonCoefficient("CT_BD_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1050105, 3232, 0, 0));
}

C_mix_bd_3::C_mix_bd_3(): WilsonCoefficient("C_BD_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1050105, 7171, 0, 0));
}

C_mix_bd_3_tilde::C_mix_bd_3_tilde(): WilsonCoefficient("CT_BD_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1050105, 7272, 0, 0));
}

C_mix_bd_4::C_mix_bd_4(): WilsonCoefficient("C_BD_4", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1050105, 3132, 0, 0));
}

C_mix_bd_5::C_mix_bd_5(): WilsonCoefficient("C_BD_5", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1050105, 7172, 0, 0));
}

C_mix_bs_1::C_mix_bs_1() : WilsonCoefficient("C_BS_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {ParameterType::SM, "SMINPUTS", 2},                   // G_F
            {ParameterType::SM, "MASS", 24},                      // M_W
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(3050305, 4141, 0, 0)
    };
}

complex_t C_mix_bs_1::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    double G_F = src.at({ParameterType::SM, "SMINPUTS", 2})->get_val();
    double M_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
    complex_t V_ts = src.at({ParameterType::SM, "VCKM", {2, 1}})->get_val();
    complex_t V_tb = src.at({ParameterType::SM, "VCKM", {2, 2}})->get_val();
    return pow(G_F * M_W * std::conj(V_ts) * V_tb / (2 * PI), 2) * S0(xt);
}

C_mix_bs_1_tilde::C_mix_bs_1_tilde() : WilsonCoefficient("CT_BS_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(3050305, 4242, 0, 0));
}

C_mix_bs_2::C_mix_bs_2() : WilsonCoefficient("C_BS_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(3050305, 3131, 0, 0));
}

C_mix_bs_2_tilde::C_mix_bs_2_tilde() : WilsonCoefficient("CT_BS_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(3050305, 3232, 0, 0));
}

C_mix_bs_3::C_mix_bs_3() : WilsonCoefficient("C_BS_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(3050305, 7171, 0, 0));
}

C_mix_bs_3_tilde::C_mix_bs_3_tilde() : WilsonCoefficient("CT_BS_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(3050305, 7272, 0, 0));
}

C_mix_bs_4::C_mix_bs_4() : WilsonCoefficient("C_BS_4", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(3050305, 3132, 0, 0));
}

C_mix_bs_5::C_mix_bs_5() : WilsonCoefficient("C_BS_5", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(3050305, 7172, 0, 0));
}

C_mix_sd_1::C_mix_sd_1() : WilsonCoefficient("C_SD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {ParameterType::SM, "SMINPUTS", 2},                   // G_F
            {ParameterType::SM, "MASS", 24},                      // M_W
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
        },
        compute_LO,
        LhaID(1030103, 4141, 0, 0)
    };
}

complex_t C_mix_sd_1::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
    double G_F = src.at({ParameterType::SM, "SMINPUTS", 2})->get_val();
    double M_W = src.at({ParameterType::SM, "MASS", 24})->get_val();
    complex_t V_td = src.at({ParameterType::SM, "VCKM", {2, 0}})->get_val();
    complex_t V_ts = src.at({ParameterType::SM, "VCKM", {2, 1}})->get_val();
    return pow(G_F * M_W * abs(std::conj(V_td) * V_ts) / (2 * PI), 2) * S0(xt);
}

C_mix_sd_1_tilde::C_mix_sd_1_tilde() : WilsonCoefficient("CT_SD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1030103, 4242, 0, 0));
}

C_mix_sd_2::C_mix_sd_2() : WilsonCoefficient("C_SD_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1030103, 3131, 0, 0));
}

C_mix_sd_2_tilde::C_mix_sd_2_tilde() : WilsonCoefficient("CT_SD_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1030103, 3232, 0, 0));
}

C_mix_sd_3::C_mix_sd_3() : WilsonCoefficient("C_SD_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1030103, 7171, 0, 0));
}

C_mix_sd_3_tilde::C_mix_sd_3_tilde() : WilsonCoefficient("CT_SD_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1030103, 7272, 0, 0));
}

C_mix_sd_4::C_mix_sd_4() : WilsonCoefficient("C_SD_4", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1030103, 3132, 0, 0));
}

C_mix_sd_5::C_mix_sd_5() : WilsonCoefficient("C_SD_5", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(1030103, 7172, 0, 0));
}

C_mix_cu_1::C_mix_cu_1() : WilsonCoefficient("C_CU_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(2040204, 4141, 0, 0));
}

C_mix_cu_1_tilde::C_mix_cu_1_tilde() : WilsonCoefficient("CT_CU_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(2040204, 4242, 0, 0));
}

C_mix_cu_2::C_mix_cu_2() : WilsonCoefficient("C_CU_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(2040204, 3131, 0, 0));
}

C_mix_cu_2_tilde::C_mix_cu_2_tilde() : WilsonCoefficient("CT_CU_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(2040204, 3232, 0, 0));
}

C_mix_cu_3::C_mix_cu_3() : WilsonCoefficient("C_CU_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(2040204, 7171, 0, 0));
}

C_mix_cu_3_tilde::C_mix_cu_3_tilde() : WilsonCoefficient("CT_CU_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(2040204, 7272, 0, 0));
}

C_mix_cu_4::C_mix_cu_4() : WilsonCoefficient("C_CU_4", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(2040204, 3132, 0, 0));
}

C_mix_cu_5::C_mix_cu_5() : WilsonCoefficient("C_CU_5", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = MatchingInfo(LhaID(2040204, 7172, 0, 0));
}

