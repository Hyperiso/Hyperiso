#include "MesonMixingWilson.h"

C_mix_bd_1::C_mix_bd_1() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {ParameterType::SM, "SMINPUTS", 2},                   // G_F
            {ParameterType::SM, "MASS", 24},                      // M_W
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

complex_t C_mix_bd_1::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});;
    double G_F = src.get_val(ParameterType::SM, "SMINPUTS", 2);
    double M_W = src.get_val(ParameterType::SM, "MASS", 24);
    complex_t V_td = src.get_val(ParameterType::SM, "VCKM", {2, 0});
    complex_t V_tb = src.get_val(ParameterType::SM, "VCKM", {2, 2});
    return pow(G_F * M_W * std::conj(V_td) * V_tb / (2 * PI), 2) * S0(xt);
}

C_mix_bd_1_tilde::C_mix_bd_1_tilde(): WilsonCoefficient("CT_BD_1", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bd_2::C_mix_bd_2(): WilsonCoefficient("C_BD_2", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bd_2_tilde::C_mix_bd_2_tilde(): WilsonCoefficient("CT_BD_2", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bd_3::C_mix_bd_3(): WilsonCoefficient("C_BD_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bd_3_tilde::C_mix_bd_3_tilde(): WilsonCoefficient("CT_BD_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bd_4::C_mix_bd_4(): WilsonCoefficient("C_BD_4", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bd_5::C_mix_bd_5(): WilsonCoefficient("C_BD_5", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bs_1::C_mix_bs_1() : WilsonCoefficient("C_BS_1", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {ParameterType::SM, "SMINPUTS", 2},                   // G_F
            {ParameterType::SM, "MASS", 24},                      // M_W
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

complex_t C_mix_bs_1::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});;
    double G_F = src.get_val(ParameterType::SM, "SMINPUTS", 2);
    double M_W = src.get_val(ParameterType::SM, "MASS", 24);
    complex_t V_ts = src.get_val(ParameterType::SM, "VCKM", {2, 1});
    complex_t V_tb = src.get_val(ParameterType::SM, "VCKM", {2, 2});
    return pow(G_F * M_W * std::conj(V_ts) * V_tb / (2 * PI), 2) * S0(xt);
}

C_mix_bs_1_tilde::C_mix_bs_1_tilde() : WilsonCoefficient("CT_BS_1", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bs_2::C_mix_bs_2() : WilsonCoefficient("C_BS_2", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bs_2_tilde::C_mix_bs_2_tilde() : WilsonCoefficient("CT_BS_2", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bs_3::C_mix_bs_3() : WilsonCoefficient("C_BS_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bs_3_tilde::C_mix_bs_3_tilde() : WilsonCoefficient("CT_BS_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bs_4::C_mix_bs_4() : WilsonCoefficient("C_BS_4", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_bs_5::C_mix_bs_5() : WilsonCoefficient("C_BS_5", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_sd_1::C_mix_sd_1() : WilsonCoefficient("C_SD_1", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", LhaID(2, 1)},  // x_t
            {ParameterType::SM, "SMINPUTS", 2},                   // G_F
            {ParameterType::SM, "MASS", 24},                      // M_W
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

complex_t C_mix_sd_1::compute_LO(const ParamSrc& src) {
    double xt = src.get_val(ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1});;
    double G_F = src.get_val(ParameterType::SM, "SMINPUTS", 2);
    double M_W = src.get_val(ParameterType::SM, "MASS", 24);
    complex_t V_td = src.get_val(ParameterType::SM, "VCKM", {2, 0});
    complex_t V_ts = src.get_val(ParameterType::SM, "VCKM", {2, 1});
    return pow(G_F * M_W * abs(std::conj(V_td) * V_ts) / (2 * PI), 2) * S0(xt);
}

C_mix_sd_1_tilde::C_mix_sd_1_tilde() : WilsonCoefficient("CT_SD_1", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_sd_2::C_mix_sd_2() : WilsonCoefficient("C_SD_2", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_sd_2_tilde::C_mix_sd_2_tilde() : WilsonCoefficient("CT_SD_2", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_sd_3::C_mix_sd_3() : WilsonCoefficient("C_SD_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_sd_3_tilde::C_mix_sd_3_tilde() : WilsonCoefficient("CT_SD_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_sd_4::C_mix_sd_4() : WilsonCoefficient("C_SD_4", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_sd_5::C_mix_sd_5() : WilsonCoefficient("C_SD_5", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_cu_1::C_mix_cu_1() : WilsonCoefficient("C_CU_1", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_cu_1_tilde::C_mix_cu_1_tilde() : WilsonCoefficient("CT_CU_1", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_cu_2::C_mix_cu_2() : WilsonCoefficient("C_CU_2", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_cu_2_tilde::C_mix_cu_2_tilde() : WilsonCoefficient("CT_CU_2", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_cu_3::C_mix_cu_3() : WilsonCoefficient("C_CU_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_cu_3_tilde::C_mix_cu_3_tilde() : WilsonCoefficient("CT_CU_3", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_cu_4::C_mix_cu_4() : WilsonCoefficient("C_CU_4", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

C_mix_cu_5::C_mix_cu_5() : WilsonCoefficient("C_CU_5", GroupMapper::str(WGroup::MESON_MIXING, ScaleType::MATCHING)) {
    matching_info[QCDOrder::LO] = MatchingInfo(get_lhaid_from_name(QCDOrder::LO));
}

