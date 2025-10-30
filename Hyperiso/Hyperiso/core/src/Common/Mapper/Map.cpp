

#include "Map.h"

const std::map<Observables, std::string>& observable_mapping() {
    static const std::map<Observables, std::string> m = {
        {Observables::BR_BS_MUMU, "BR_Bs__mu_mu"},
        {Observables::BR_BS_MUMU_UNTAG, "BRuntag_Bs__mu_mu"},
        {Observables::BR_BD_MUMU, "BR_Bd__mu_mu"},
        {Observables::BR_BU_TAU_NU, "BR_Bu__tau_nu"},
        {Observables::R_TAU_NU, "R_tau_nu"},
        {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, "IA_B__K*_gamma"},
        {Observables::BR_B__KSTAR_GAMMA, "BR_B__K*_gamma"},
        {Observables::BR_B_XS_GAMMA, "BR_B__Xs_gamma"},
        {Observables::BR_B__D_TAU_NU, "BR_B__D_tau_nu"},
        {Observables::A_FB_B__D_TAU_NU, "A_FB_B__D_tau_nu"},
        {Observables::P_TAU_B__D_TAU_NU, "P_tau_B__D_tau_nu"},
        {Observables::R_D, "R_D"},
        {Observables::BR_B__DSTAR_TAU_NU, "BR_B__D*_tau_nu"},
        {Observables::A_FB_B__DSTAR_TAU_NU, "A_FB_B__D*_tau_nu"},
        {Observables::P_TAU_B__DSTAR_TAU_NU, "P_tau_B__D*_tau_nu"},
        {Observables::P_D_B__DSTAR_TAU_NU, "P_D*_B__D*_tau_nu"},
        {Observables::R_DSTAR, "R_D*"},
        {Observables::BR_B__Xs_mu_mu__LOW_Q2, "BR_B__Xs_mu_mu__[1-6]"},
        {Observables::BR_B__Xs_mu_mu__HIGH_Q2, "BR_B__Xs_mu_mu__[>14.4]"},
        {Observables::BR_B__Xs_tau_tau__HIGH_Q2, "BR_B__Xs_tau_tau__[>14.4]"},
        {Observables::DGAMMA_DQ2_B__KSTAR_L_L, "dGamma_dq2_B__K*_l_l"},
        {Observables::DGAMMA_BAR_DQ2_B__KSTAR_L_L, "dGamma_bar_dq2_B__K*_l_l"},
        {Observables::A_FB_B__KSTAR_L_L, "A_FB_B__K*_l_l"},
        {Observables::Q0_A_FB_B__KSTAR_L_L, "q0_A_FB_B__K*_l_l"},
        {Observables::A_CP_B__KSTAR_L_L, "A_CP_B__K*_l_l"},
        {Observables::F_L_B__KSTAR_L_L, "F_L_B__K*_l_l"},
        {Observables::F_T_B__KSTAR_L_L, "F_T_B__K*_l_l"},
        {Observables::A_T_1_B__KSTAR_L_L, "A_T_1_B__K*_l_l"},
        {Observables::A_T_2_B__KSTAR_L_L, "A_T_2_B__K*_l_l"},
        {Observables::A_T_3_B__KSTAR_L_L, "A_T_3_B__K*_l_l"},
        {Observables::A_T_4_B__KSTAR_L_L, "A_T_4_B__K*_l_l"},
        {Observables::A_T_5_B__KSTAR_L_L, "A_T_5_B__K*_l_l"},
        {Observables::A_T_RE_B__KSTAR_L_L, "A_T_RE_B__K*_l_l"},
        {Observables::A_T_RE_CPV_B__KSTAR_L_L, "A_T_RE_CPV_B__K*_l_l"},
        {Observables::A_IM_B__KSTAR_L_L, "A_IM_B__K*_l_l"},
        {Observables::ALPHA_K_B__KSTAR_L_L, "alpha_K_B__K*_l_l"},
        {Observables::H_T_1_B__KSTAR_L_L, "H_T_1_B__K*_l_l"},
        {Observables::H_T_2_B__KSTAR_L_L, "H_T_2_B__K*_l_l"},
        {Observables::H_T_3_B__KSTAR_L_L, "H_T_3_B__K*_l_l"},
        {Observables::P_2_B__KSTAR_L_L, "P_2_B__K*_l_l"},
        {Observables::P_3_B__KSTAR_L_L, "P_3_B__K*_l_l"},
        {Observables::P_6_B__KSTAR_L_L, "P_6_B__K*_l_l"},
        {Observables::P_8_B__KSTAR_L_L, "P_8_B__K*_l_l"},
        {Observables::P_PRIME_4_B__KSTAR_L_L, "P'_4_B__K*_l_l"},
        {Observables::P_PRIME_5_B__KSTAR_L_L, "P'_5_B__K*_l_l"},
        {Observables::P_PRIME_6_B__KSTAR_L_L, "P'_6_B__K*_l_l"},
        {Observables::P_PRIME_8_B__KSTAR_L_L, "P'_8_B__K*_l_l"},
        {Observables::S_3_B__KSTAR_L_L, "S_3_B__K*_l_l"},
        {Observables::S_4_B__KSTAR_L_L, "S_4_B__K*_l_l"},
        {Observables::S_5_B__KSTAR_L_L, "S_5_B__K*_l_l"},
        {Observables::S_6C_B__KSTAR_L_L, "S_6c_B__K*_l_l"},
        {Observables::S_7_B__KSTAR_L_L, "S_7_B__K*_l_l"},
        {Observables::S_8_B__KSTAR_L_L, "S_8_B__K*_l_l"},
        {Observables::S_9_B__KSTAR_L_L, "S_9_B__K*_l_l"},
        {Observables::A_3_B__KSTAR_L_L, "A_3_B__K*_l_l"},
        {Observables::A_4_B__KSTAR_L_L, "A_4_B__K*_l_l"},
        {Observables::A_5_B__KSTAR_L_L, "A_5_B__K*_l_l"},
        {Observables::A_6S_B__KSTAR_L_L, "A_6s_B__K*_l_l"},
        {Observables::A_7_B__KSTAR_L_L, "A_7_B__K*_l_l"},
        {Observables::A_8_B__KSTAR_L_L, "A_8_B__K*_l_l"},
        {Observables::A_9_B__KSTAR_L_L, "A_9_B__K*_l_l"},
        {Observables::P_1_CPV_B__KSTAR_L_L, "P_1_CPV_B__K*_l_l"},
        {Observables::P_2_CPV_B__KSTAR_L_L, "P_2_CPV_B__K*_l_l"},
        {Observables::P_3_CPV_B__KSTAR_L_L, "P_3_CPV_B__K*_l_l"},
        {Observables::P_PRIME_4_CPV_B__KSTAR_L_L, "P'_4_CPV_B__K*_l_l"},
        {Observables::P_PRIME_5_CPV_B__KSTAR_L_L, "P'_5_CPV_B__K*_l_l"},
        {Observables::P_PRIME_6_CPV_B__KSTAR_L_L, "P'_6_CPV_B__K*_l_l"},
        {Observables::P_PRIME_8_CPV_B__KSTAR_L_L, "P'_8_CPV_B__K*_l_l"},
        {Observables::DGAMMA_DQ2_BS__PHI_L_L, "dGamma_dq2_Bs__phi_l_l"},
        {Observables::DGAMMA_BAR_DQ2_BS__PHI_L_L, "dGamma_bar_dq2_Bs__phi_l_l"},
        {Observables::A_FB_CPV_BS__PHI_L_L, "A_FB_CPV_Bs__phi_l_l"},
        {Observables::F_L_BS_PHI_L_L, "F_L_Bs__phi_l_l"},
        {Observables::A_T_2_BS_PHI_L_L, "A_T_2_Bs__phi_l_l"},
        {Observables::A_T_RE_CPV_BS_PHI_L_L, "A_T_Re_CPV_Bs__phi_l_l"},
        {Observables::A_T_IM_CPV_BS_PHI_L_L, "A_T_Im_CPV_Bs__phi_l_l"},
        {Observables::P_PRIME_4_BS_PHI_L_L, "P'_4_Bs__phi_l_l"},
        {Observables::P_PRIME_6_BS_PHI_L_L, "P'_6_Bs__phi_l_l"},
        {Observables::S_2S_BS_PHI_L_L, "S_2s_Bs__phi_l_l"},
        {Observables::S_3_BS_PHI_L_L, "S_3_Bs__phi_l_l"},
        {Observables::S_4_BS_PHI_L_L, "S_4_Bs__phi_l_l"},
        {Observables::S_7_BS_PHI_L_L, "S_7_Bs__phi_l_l"},
        {Observables::A_5_BS_PHI_L_L, "A_5_Bs__phi_l_l"},
        {Observables::A_6C_BS_PHI_L_L, "A_6c_Bs__phi_l_l"},
        {Observables::A_8_BS_PHI_L_L, "A_8_Bs__phi_l_l"},
        {Observables::A_9_BS_PHI_L_L, "A_9_Bs__phi_l_l"},
        {Observables::P_2_CPV_BS_PHI_L_L, "P_2_CPV_Bs__phi_l_l"},
        {Observables::P_3_CPV_BS_PHI_L_L, "P_3_CPV_Bs__phi_l_l"},
        {Observables::P_PRIME_5_CPV_BS_PHI_L_L, "P'_5_CPV_Bs__phi_l_l"},
        {Observables::P_PRIME_8_CPV_BS_PHI_L_L, "P'_8_CPV_Bs__phi_l_l"},
        {Observables::Q_8M_BS_PHI_L_L, "Q_8m_Bs__phi_l_l"},
        {Observables::Q_8P_BS_PHI_L_L, "Q_8p_CPV_Bs__phi_l_l"},
        {Observables::Q_9_BS_PHI_L_L, "Q_9_CPV_Bs__phi_l_l"},
        {Observables::DGAMMA_DQ2_B__K_L_L, "dGamma_dq2_B__K_l_l"},
        {Observables::A_FB_B__K_L_L, "A_FB_B__K_l_l"},
        {Observables::F_H_B__K_L_L, "F_H_B__K_l_l"},
        {Observables::PHI_D, "phi_d"},
        {Observables::DELTA_M_BD, "Delta_M_Bd"},
        {Observables::PHI_S, "phi_s"},
        {Observables::DELTA_M_BS, "Delta_M_Bs"},
        {Observables::A_FS, "a_fs"},
        {Observables::DELTA_M_K, "Delta_M_K"},
        {Observables::ABS_EPSILON_K, "|epsilon_K|"},
        {Observables::X_D, "x_D"},
        {Observables::BR_KL__MU_MU, "BR_K_L__mu_mu"},
        {Observables::BR_KS__MU_MU, "BR_K_S__mu_mu"},
        {Observables::TEST, "test"},
    };
    return m;
}

const std::map<Observables, LhaID>& observable_flha_mapping() {
    static const std::map<Observables, LhaID> m = {
        {Observables::BR_BS_MUMU,        LhaID(531, 1, 2, 13, -13)},
        {Observables::BR_BS_MUMU_UNTAG,  LhaID(531, 15, 2, 13, -13)},
        {Observables::BR_BD_MUMU,        LhaID(511, 1, 2, 13, -13)},
        {Observables::BR_BU_TAU_NU,      LhaID(521, 1, 2, -15, 16)},
        {Observables::R_TAU_NU,          LhaID(521, 2, 2, -15, 16)},
        {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, LhaID(521, 4, 2, 313, 22)},
        {Observables::BR_B_XS_GAMMA,     LhaID(5,   1, 2, 3, 22)},
        {Observables::BR_B__D_TAU_NU,    LhaID(521, 1, 3, 421, -15, 16)},
        {Observables::A_FB_B__D_TAU_NU,  LhaID(521, 5, 3, 421, -15, 16)},
        {Observables::P_TAU_B__D_TAU_NU, LhaID(521, 82,3, 421, -15, 16)},
        {Observables::R_D,               LhaID(521, 11,3, 421, -15, 16)},
        {Observables::BR_B__DSTAR_TAU_NU,LhaID(521, 1, 3, 423, -15, 16)},
        {Observables::A_FB_B__DSTAR_TAU_NU, LhaID(521, 5, 3, 423, -15, 16)},
        {Observables::P_TAU_B__DSTAR_TAU_NU, LhaID(521, 82,3, 423, -15, 16)},
        {Observables::P_D_B__DSTAR_TAU_NU,   LhaID(521, 81,3, 423, -15, 16)},
        {Observables::R_DSTAR,           LhaID(521, 11,3, 423, -15, 16)},
    };
    return m;
}

const std::map<QCDOrder, std::string>& order_mapping() {
    static const std::map<QCDOrder, std::string> m = {
        {QCDOrder::NONE, "None"},
        {QCDOrder::LO, "LO"},
        {QCDOrder::NLO, "NLO"},
        {QCDOrder::NNLO, "NNLO"},
    };
    return m;
}

const std::map<WGroup, std::string>& group_mapping() {
    static const std::map<WGroup, std::string> m = {
        {WGroup::B, "BCoefficients"},
        {WGroup::BPrime, "BPrimeCoefficients"},
        {WGroup::BScalar, "BScalarCoefficients"},
        {WGroup::CC_bc, "BcChargedCurrentCoefficients"},
        {WGroup::CC_bu, "BuChargedCurrentCoefficients"},
        {WGroup::CC_cs, "DsChargedCurrentCoefficients"},
        {WGroup::CC_cd, "DdChargedCurrentCoefficients"},
        {WGroup::CC_su, "KuChargedCurrentCoefficients"},
        {WGroup::MESON_MIXING, "MesonMixing"},
        {WGroup::CUSTOM_GROUP, "Custom_Group"},
        {WGroup::K, "K"},
    };
    return m;
}



const std::map<WCoef, std::string>& wcoef_mapping() {
    static const std::map<WCoef, std::string> m = {
        {WCoef::C1, "C1"}, {WCoef::C2, "C2"}, {WCoef::C3, "C3"}, {WCoef::C4, "C4"}, {WCoef::C5, "C5"}, {WCoef::C6, "C6"},
        {WCoef::C7, "C7"}, {WCoef::C8, "C8"}, {WCoef::C9, "C9"}, {WCoef::C10, "C10"},
        {WCoef::CP1, "CP1"}, {WCoef::CP2, "CP2"}, {WCoef::CP3, "CP3"}, {WCoef::CP4, "CP4"},
        {WCoef::CP5, "CP5"}, {WCoef::CP6, "CP6"}, {WCoef::CP7, "CP7"}, {WCoef::CP8, "CP8"},
        {WCoef::CP9, "CP9"}, {WCoef::CP10, "CP10"}, {WCoef::CQ1, "CQ1"}, {WCoef::CQ2, "CQ2"},
        {WCoef::CPQ1, "CPQ1"}, {WCoef::CPQ2, "CPQ2"},
        {WCoef::C_V1_bc, "C_V1_bc"}, {WCoef::C_V2_bc, "C_V2_bc"}, {WCoef::C_S1_bc, "C_S1_bc"}, {WCoef::C_S2_bc, "C_S2_bc"}, {WCoef::C_T_bc, "C_T_bc"},
        {WCoef::C_V1_bu, "C_V1_bu"}, {WCoef::C_V2_bu, "C_V2_bu"}, {WCoef::C_S1_bu, "C_S1_bu"}, {WCoef::C_S2_bu, "C_S2_bu"}, {WCoef::C_T_bu, "C_T_bu"},
        {WCoef::C_V1_cs, "C_V1_cs"}, {WCoef::C_V2_cs, "C_V2_cs"}, {WCoef::C_S1_cs, "C_S1_cs"}, {WCoef::C_S2_cs, "C_S2_cs"}, {WCoef::C_T_cs, "C_T_cs"},
        {WCoef::C_V1_cd, "C_V1_cd"}, {WCoef::C_V2_cd, "C_V2_cd"}, {WCoef::C_S1_cd, "C_S1_cd"}, {WCoef::C_S2_cd, "C_S2_cd"}, {WCoef::C_T_cd, "C_T_cd"},
        {WCoef::C_V1_su, "C_V1_su"}, {WCoef::C_V2_su, "C_V2_su"}, {WCoef::C_S1_su, "C_S1_su"}, {WCoef::C_S2_su, "C_S2_su"}, {WCoef::C_T_su, "C_T_su"},
        {WCoef::C_BD_1, "C_BD_1"}, {WCoef::CT_BD_1, "CT_BD_1"}, {WCoef::C_BD_2, "C_BD_2"}, {WCoef::CT_BD_2, "CT_BD_2"}, {WCoef::C_BD_3, "C_BD_3"}, {WCoef::CT_BD_3, "CT_BD_3"}, {WCoef::C_BD_4, "C_BD_4"}, {WCoef::C_BD_5, "C_BD_5"},
        {WCoef::C_BS_1, "C_BS_1"}, {WCoef::CT_BS_1, "CT_BS_1"}, {WCoef::C_BS_2, "C_BS_2"}, {WCoef::CT_BS_2, "CT_BS_2"}, {WCoef::C_BS_3, "C_BS_3"}, {WCoef::CT_BS_3, "CT_BS_3"}, {WCoef::C_BS_4, "C_BS_4"}, {WCoef::C_BS_5, "C_BS_5"},
        {WCoef::C_SD_1, "C_SD_1"}, {WCoef::CT_SD_1, "CT_SD_1"}, {WCoef::C_SD_2, "C_SD_2"}, {WCoef::CT_SD_2, "CT_SD_2"}, {WCoef::C_SD_3, "C_SD_3"}, {WCoef::CT_SD_3, "CT_SD_3"}, {WCoef::C_SD_4, "C_SD_4"}, {WCoef::C_SD_5, "C_SD_5"},
        {WCoef::C_CU_1, "C_CU_1"}, {WCoef::CT_CU_1, "CT_CU_1"}, {WCoef::C_CU_2, "C_CU_2"}, {WCoef::CT_CU_2, "CT_CU_2"}, {WCoef::C_CU_3, "C_CU_3"}, {WCoef::CT_CU_3, "CT_CU_3"}, {WCoef::C_CU_4, "C_CU_4"}, {WCoef::C_CU_5, "C_CU_5"},
        {WCoef::CK9, "CK9"}, {WCoef::CPK9, "CPK9"}, {WCoef::CK10, "CK10"}, {WCoef::CPK10, "CPK10"}, {WCoef::CKQ1, "CKQ1"}, {WCoef::CKQ2, "CKQ2"}, {WCoef::CPKQ1, "CPKQ1"}, {WCoef::CPKQ2, "CPKQ2"}, {WCoef::CK_L, "CK_L"}
    };
    return m;
}


const std::map<WCoef, std::pair<int, int>>& wcoef_flha_mapping() {
    static const std::map<WCoef, std::pair<int, int>> m = {
        {WCoef::C1, {3040405, 6161}}, {WCoef::C2, {3040405, 4141}}, {WCoef::C3, {3050707, 4133}},
        {WCoef::C4, {3050707, 6153}}, {WCoef::C5, {3050707, 4536}}, {WCoef::C6, {3050707, 6556}},
        {WCoef::C7, {305, 4422}}, {WCoef::C8, {305, 6421}}, {WCoef::C9, {3051313, 4133}}, {WCoef::C10, {3051313, 4137}},
        {WCoef::CP1, {3040405, 6262}}, {WCoef::CP2, {3040405, 4242}}, {WCoef::CP3, {3050707, 4233}},
        {WCoef::CP4, {3050707, 6253}}, {WCoef::CP5, {3050707, 4636}}, {WCoef::CP6, {3050707, 6656}},
        {WCoef::CP7, {305, 4322}}, {WCoef::CP8, {305, 4321}}, {WCoef::CP9, {3051313, 4233}},
        {WCoef::CP10, {3051313, 4234}}, {WCoef::CQ1, {3051313, 3230}}, {WCoef::CQ2, {3051313, 3233}},
        {WCoef::CPQ1, {3051313, 3130}}, {WCoef::CPQ2, {3051313, 3133}},
        {WCoef::C_V1_bc, {4051516, 4141}}, {WCoef::C_V2_bc, {4051516, 4241}},
        {WCoef::C_S1_bc, {4051516, 3231}}, {WCoef::C_S2_bc, {4051516, 3131}}, {WCoef::C_T_bc, {4051516, 4343}},
        {WCoef::C_V1_bu, {4051516, 1}}, {WCoef::C_V2_bu, {4051516, 2}},
        {WCoef::C_S1_bu, {4051516, 3}}, {WCoef::C_S2_bu, {4051516, 4}}, {WCoef::C_T_bu, {4051516, 5}},
        {WCoef::C_V1_cs, {4051516, 11}}, {WCoef::C_V2_cs, {4051516, 22}},
        {WCoef::C_S1_cs, {4051516, 33}}, {WCoef::C_S2_cs, {4051516, 44}}, {WCoef::C_T_cs, {4051516, 55}},
        {WCoef::C_V1_cd, {4051516, 111}}, {WCoef::C_V2_cd, {4051516, 222}},
        {WCoef::C_S1_cd, {4051516, 333}}, {WCoef::C_S2_cd, {4051516, 444}}, {WCoef::C_T_cd, {4051516, 555}},
        {WCoef::C_V1_su, {4051516, 1111}}, {WCoef::C_V2_su, {4051516, 2222}},
        {WCoef::C_S1_su, {4051516, 3333}}, {WCoef::C_S2_su, {4051516, 4444}}, {WCoef::C_T_su, {4051516, 5555}},
        {WCoef::C_BD_1, {1050105, 4141}}, {WCoef::CT_BD_1, {1050105, 4242}}, {WCoef::C_BD_2, {1050105, 3131}}, {WCoef::CT_BD_2, {1050105, 3232}},
        {WCoef::C_BD_3, {1050105, 7171}}, {WCoef::CT_BD_3, {1050105, 7272}}, {WCoef::C_BD_4, {1050105, 3132}}, {WCoef::C_BD_5, {1050105, 7172}},
        {WCoef::C_BS_1, {3050305, 4141}}, {WCoef::CT_BS_1, {3050305, 4242}}, {WCoef::C_BS_2, {3050305, 3131}}, {WCoef::CT_BS_2, {3050305, 3232}},
        {WCoef::C_BS_3, {3050305, 7171}}, {WCoef::CT_BS_3, {3050305, 7272}}, {WCoef::C_BS_4, {3050305, 3132}}, {WCoef::C_BS_5, {3050305, 7172}},
        {WCoef::C_SD_1, {1030103, 4141}}, {WCoef::CT_SD_1, {1030103, 4242}}, {WCoef::C_SD_2, {1030103, 3131}}, {WCoef::CT_SD_2, {1030103, 3232}},
        {WCoef::C_SD_3, {1030103, 7171}}, {WCoef::CT_SD_3, {1030103, 7272}}, {WCoef::C_SD_4, {1030103, 3132}}, {WCoef::C_SD_5, {1030103, 7172}},
        {WCoef::C_CU_1, {2040204, 4141}}, {WCoef::CT_CU_1, {2040204, 4242}}, {WCoef::C_CU_2, {2040204, 3131}}, {WCoef::CT_CU_2, {2040204, 3232}},
        {WCoef::C_CU_3, {2040204, 7171}}, {WCoef::CT_CU_3, {2040204, 7272}}, {WCoef::C_CU_4, {2040204, 3132}}, {WCoef::C_CU_5, {2040204, 7172}},

        {WCoef::CK9, {0, 1}}, {WCoef::CPK9, {0, 2}}, {WCoef::CK10, {0, 3}}, {WCoef::CPK10, {0, 4}}, //TODO : real coefs
        {WCoef::CKQ1, {0, 5}}, {WCoef::CKQ2, {0, 6}}, {WCoef::CPKQ1, {0, 7}}, {WCoef::CPKQ2, {0, 8}},
        {WCoef::CK_L, {0, 9}},
    };
    return m;
}


const std::map<ParameterType, std::string>& parametertype_mapping() {
    static const std::map<ParameterType, std::string> m = {
        {ParameterType::SM, "SM"},
        {ParameterType::BSM, "BSM"},
        {ParameterType::FLAVOR, "FLAVOR"},
        {ParameterType::WILSON, "WILSON"},
        {ParameterType::DECAY, "DECAY"},
        {ParameterType::PASSTHROUGH, "PASSTHROUGH"},
        {ParameterType::OBSERVABLE, "OBSERVABLE"},
    };
    return m;
}

const std::map<Model, std::string>& model_mapping() {
    static const std::map<Model, std::string> m = {
        {Model::SM, "SM"},
        {Model::SUSY, "SUSY"},
        {Model::THDM, "THDM"},
        {Model::CUSTOM, "CUSTOM"},
    };
    return m;
}

const std::map<WilsonBasis, std::string>& wilsonbasis_mapping() {
    static const std::map<WilsonBasis, std::string> m = {
        {WilsonBasis::B_STANDARD, "STANDARD"},
        {WilsonBasis::B_TRADITIONAL, "TRADITIONAL"},
    };
    return m;
}

const std::map<ContributionType, std::string>& contributiontype_mapping() {
    static const std::map<ContributionType, std::string> m = {
        {ContributionType::SM, "SM"},
        {ContributionType::BSM, "BSM"},
        {ContributionType::TOTAL, "TOTAL"},
    };
    return m;
}

const std::map<Decays, std::string>& decays_mapping() {
    static const std::map<Decays, std::string> m = {
        {Decays::B__D_l_nu, "B__D_l_nu"},
        {Decays::B__Dstar_l_nu, "B__Dstar_l_nu"},
        {Decays::B__Kstar_gamma, "B__Kstar_gamma"},
        {Decays::B__l_l, "B__l_l"},
        {Decays::B__l_nu, "B__l_nu"},
        {Decays::B__Xs, "B__Xs"},
        {Decays::B__Xs_l_l, "B__Xs_ll"},
        {Decays::B__Kstar_l_l, "B__K*_l_l"},
        {Decays::B__K_l_l, "B__K_l_l"},
        {Decays::Bs__phi_l_l, "Bs__phi_l_l"},
        {Decays::M0_Mix, "M0_Mix"},
        {Decays::K__l_l, "K__l_l"},
    };
    return m;
}

const std::map<Decays, std::vector<Observables>>& decay_observable_mapping() {
    static const std::map<Decays, std::vector<Observables>> m = {
        {Decays::B__D_l_nu,     {Observables::BR_B__D_TAU_NU, Observables::A_FB_B__D_TAU_NU, Observables::P_TAU_B__D_TAU_NU, Observables::R_D}},
        {Decays::B__Dstar_l_nu, {Observables::BR_B__DSTAR_TAU_NU, Observables::A_FB_B__DSTAR_TAU_NU, Observables::P_TAU_B__DSTAR_TAU_NU, Observables::P_D_B__DSTAR_TAU_NU, Observables::R_DSTAR}},
        {Decays::B__Kstar_gamma,{Observables::BR_B__KSTAR_GAMMA, Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA}},
        {Decays::B__l_l,        {Observables::BR_BS_MUMU, Observables::BR_BD_MUMU, Observables::BR_BS_MUMU_UNTAG}},
        {Decays::B__l_nu,       {Observables::BR_BU_TAU_NU, Observables::R_TAU_NU}},
        {Decays::B__Xs,         {Observables::BR_B_XS_GAMMA}},
        {Decays::B__Xs_l_l,     {Observables::BR_B__Xs_mu_mu__LOW_Q2, Observables::BR_B__Xs_mu_mu__HIGH_Q2, Observables::BR_B__Xs_tau_tau__HIGH_Q2}},
        {Decays::B__Kstar_l_l,  {Observables::DGAMMA_DQ2_B__KSTAR_L_L, Observables::A_FB_B__KSTAR_L_L, Observables::Q0_A_FB_B__KSTAR_L_L, Observables::A_CP_B__KSTAR_L_L, Observables::F_L_B__KSTAR_L_L, Observables::F_T_B__KSTAR_L_L, Observables::A_T_1_B__KSTAR_L_L, Observables::A_T_2_B__KSTAR_L_L, Observables::A_T_3_B__KSTAR_L_L, Observables::A_T_4_B__KSTAR_L_L, Observables::A_T_5_B__KSTAR_L_L, Observables::A_T_RE_B__KSTAR_L_L, Observables::A_T_RE_CPV_B__KSTAR_L_L, Observables::A_IM_B__KSTAR_L_L, Observables::ALPHA_K_B__KSTAR_L_L, Observables::H_T_1_B__KSTAR_L_L, Observables::H_T_2_B__KSTAR_L_L, Observables::H_T_3_B__KSTAR_L_L, Observables::P_2_B__KSTAR_L_L, Observables::P_3_B__KSTAR_L_L, Observables::P_6_B__KSTAR_L_L, Observables::P_8_B__KSTAR_L_L, Observables::P_PRIME_4_B__KSTAR_L_L, Observables::P_PRIME_5_B__KSTAR_L_L, Observables::P_PRIME_6_B__KSTAR_L_L, Observables::P_PRIME_8_B__KSTAR_L_L, Observables::S_3_B__KSTAR_L_L, Observables::S_4_B__KSTAR_L_L, Observables::S_5_B__KSTAR_L_L, Observables::S_6C_B__KSTAR_L_L, Observables::S_7_B__KSTAR_L_L, Observables::S_8_B__KSTAR_L_L, Observables::S_9_B__KSTAR_L_L, Observables::A_3_B__KSTAR_L_L, Observables::A_4_B__KSTAR_L_L, Observables::A_5_B__KSTAR_L_L, Observables::A_6S_B__KSTAR_L_L, Observables::A_7_B__KSTAR_L_L, Observables::A_8_B__KSTAR_L_L, Observables::A_9_B__KSTAR_L_L, Observables::P_1_CPV_B__KSTAR_L_L, Observables::P_2_CPV_B__KSTAR_L_L, Observables::P_3_CPV_B__KSTAR_L_L, Observables::P_PRIME_4_CPV_B__KSTAR_L_L, Observables::P_PRIME_5_CPV_B__KSTAR_L_L, Observables::P_PRIME_6_CPV_B__KSTAR_L_L, Observables::P_PRIME_8_CPV_B__KSTAR_L_L}},
        {Decays::B__K_l_l,      {Observables::DGAMMA_DQ2_B__K_L_L, Observables::A_FB_B__K_L_L, Observables::F_H_B__K_L_L}},
        {Decays::Bs__phi_l_l,   {Observables::DGAMMA_DQ2_BS__PHI_L_L, Observables::DGAMMA_BAR_DQ2_BS__PHI_L_L, Observables::A_FB_CPV_BS__PHI_L_L, Observables::F_L_BS_PHI_L_L, Observables::A_T_2_BS_PHI_L_L, Observables::A_T_RE_CPV_BS_PHI_L_L, Observables::A_T_IM_CPV_BS_PHI_L_L, Observables::P_PRIME_4_BS_PHI_L_L, Observables::P_PRIME_6_BS_PHI_L_L, Observables::S_2S_BS_PHI_L_L, Observables::S_3_BS_PHI_L_L, Observables::S_4_BS_PHI_L_L, Observables::S_7_BS_PHI_L_L, Observables::A_5_BS_PHI_L_L, Observables::A_6C_BS_PHI_L_L, Observables::A_8_BS_PHI_L_L, Observables::A_9_BS_PHI_L_L, Observables::P_2_CPV_BS_PHI_L_L, Observables::P_3_CPV_BS_PHI_L_L, Observables::P_PRIME_5_CPV_BS_PHI_L_L, Observables::P_PRIME_8_CPV_BS_PHI_L_L, Observables::Q_8M_BS_PHI_L_L, Observables::Q_8P_BS_PHI_L_L, Observables::Q_9_BS_PHI_L_L}},
        {Decays::M0_Mix,        {Observables::PHI_D, Observables::DELTA_M_BD, Observables::PHI_S, Observables::DELTA_M_BS, Observables::A_FS, Observables::DELTA_M_K, Observables::ABS_EPSILON_K, Observables::X_D}},
        {Decays::K__l_l,        {Observables::BR_KL__MU_MU, Observables::BR_KS__MU_MU}},
    };
    return m;
}

const std::map<MassType, std::string>& masstype_mapping() {
    static const std::map<MassType, std::string> m = {
        {MassType::POLE,  "POLE"},
        {MassType::MSBAR, "MSBAR"},
    };
    return m;
}

const std::map<ScaleType, std::string>& scaletype_mapping() {
    static const std::map<ScaleType, std::string> m = {
        {ScaleType::MATCHING, "EW_SCALE"},
        {ScaleType::HADRONIC, "B_SCALE"},
    };
    return m;
}

