

#include "Map.h"

using enum Observables;
using enum Decays;
using enum WCoef;
using enum WGroup;

const std::map<Observables, std::string>& observable_mapping() {
    static const std::map<Observables, std::string> m = {
        {BR_BS_MUMU, "BR_Bs__mu_mu"},
        {BR_BS_MUMU_UNTAG, "BRuntag_Bs__mu_mu"},
        {BR_BD_MUMU, "BR_Bd__mu_mu"},
        {BR_BU_TAU_NU, "BR_Bu__tau_nu"},
        {R_TAU_NU, "R_tau_nu"},
        {ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, "IA_B__K*_gamma"},
        {BR_B__KSTAR_GAMMA, "BR_B__K*_gamma"},
        {BR_B_XS_GAMMA, "BR_B__Xs_gamma"},
        {BR_B__D_TAU_NU, "BR_B__D_tau_nu"},
        {A_FB_B__D_TAU_NU, "A_FB_B__D_tau_nu"},
        {P_TAU_B__D_TAU_NU, "P_tau_B__D_tau_nu"},
        {R_D, "R_D"},
        {BR_B__DSTAR_TAU_NU, "BR_B__D*_tau_nu"},
        {A_FB_B__DSTAR_TAU_NU, "A_FB_B__D*_tau_nu"},
        {P_TAU_B__DSTAR_TAU_NU, "P_tau_B__D*_tau_nu"},
        {P_D_B__DSTAR_TAU_NU, "P_D*_B__D*_tau_nu"},
        {R_DSTAR, "R_D*"},
        {BR_B__Xs_mu_mu__LOW_Q2, "BR_B__Xs_mu_mu__[1-6]"},
        {BR_B__Xs_mu_mu__HIGH_Q2, "BR_B__Xs_mu_mu__[>14.4]"},
        {BR_B__Xs_tau_tau__HIGH_Q2, "BR_B__Xs_tau_tau__[>14.4]"},
        {DGAMMA_DQ2_B__KSTAR_L_L, "dGamma_dq2_B__K*_l_l"},
        {DGAMMA_BAR_DQ2_B__KSTAR_L_L, "dGamma_bar_dq2_B__K*_l_l"},
        {A_FB_B__KSTAR_L_L, "A_FB_B__K*_l_l"},
        {Q0_A_FB_B__KSTAR_L_L, "q0_A_FB_B__K*_l_l"},
        {A_CP_B__KSTAR_L_L, "A_CP_B__K*_l_l"},
        {F_L_B__KSTAR_L_L, "F_L_B__K*_l_l"},
        {F_T_B__KSTAR_L_L, "F_T_B__K*_l_l"},
        {A_T_1_B__KSTAR_L_L, "A_T_1_B__K*_l_l"},
        {A_T_2_B__KSTAR_L_L, "A_T_2_B__K*_l_l"},
        {A_T_3_B__KSTAR_L_L, "A_T_3_B__K*_l_l"},
        {A_T_4_B__KSTAR_L_L, "A_T_4_B__K*_l_l"},
        {A_T_5_B__KSTAR_L_L, "A_T_5_B__K*_l_l"},
        {A_T_RE_B__KSTAR_L_L, "A_T_RE_B__K*_l_l"},
        {A_T_RE_CPV_B__KSTAR_L_L, "A_T_RE_CPV_B__K*_l_l"},
        {A_IM_B__KSTAR_L_L, "A_IM_B__K*_l_l"},
        {ALPHA_K_B__KSTAR_L_L, "alpha_K_B__K*_l_l"},
        {H_T_1_B__KSTAR_L_L, "H_T_1_B__K*_l_l"},
        {H_T_2_B__KSTAR_L_L, "H_T_2_B__K*_l_l"},
        {H_T_3_B__KSTAR_L_L, "H_T_3_B__K*_l_l"},
        {P_2_B__KSTAR_L_L, "P_2_B__K*_l_l"},
        {P_3_B__KSTAR_L_L, "P_3_B__K*_l_l"},
        {P_6_B__KSTAR_L_L, "P_6_B__K*_l_l"},
        {P_8_B__KSTAR_L_L, "P_8_B__K*_l_l"},
        {P_PRIME_4_B__KSTAR_L_L, "P'_4_B__K*_l_l"},
        {P_PRIME_5_B__KSTAR_L_L, "P'_5_B__K*_l_l"},
        {P_PRIME_6_B__KSTAR_L_L, "P'_6_B__K*_l_l"},
        {P_PRIME_8_B__KSTAR_L_L, "P'_8_B__K*_l_l"},
        {S_3_B__KSTAR_L_L, "S_3_B__K*_l_l"},
        {S_4_B__KSTAR_L_L, "S_4_B__K*_l_l"},
        {S_5_B__KSTAR_L_L, "S_5_B__K*_l_l"},
        {S_6C_B__KSTAR_L_L, "S_6c_B__K*_l_l"},
        {S_7_B__KSTAR_L_L, "S_7_B__K*_l_l"},
        {S_8_B__KSTAR_L_L, "S_8_B__K*_l_l"},
        {S_9_B__KSTAR_L_L, "S_9_B__K*_l_l"},
        {A_3_B__KSTAR_L_L, "A_3_B__K*_l_l"},
        {A_4_B__KSTAR_L_L, "A_4_B__K*_l_l"},
        {A_5_B__KSTAR_L_L, "A_5_B__K*_l_l"},
        {A_6S_B__KSTAR_L_L, "A_6s_B__K*_l_l"},
        {A_7_B__KSTAR_L_L, "A_7_B__K*_l_l"},
        {A_8_B__KSTAR_L_L, "A_8_B__K*_l_l"},
        {A_9_B__KSTAR_L_L, "A_9_B__K*_l_l"},
        {P_1_CPV_B__KSTAR_L_L, "P_1_CPV_B__K*_l_l"},
        {P_2_CPV_B__KSTAR_L_L, "P_2_CPV_B__K*_l_l"},
        {P_3_CPV_B__KSTAR_L_L, "P_3_CPV_B__K*_l_l"},
        {P_PRIME_4_CPV_B__KSTAR_L_L, "P'_4_CPV_B__K*_l_l"},
        {P_PRIME_5_CPV_B__KSTAR_L_L, "P'_5_CPV_B__K*_l_l"},
        {P_PRIME_6_CPV_B__KSTAR_L_L, "P'_6_CPV_B__K*_l_l"},
        {P_PRIME_8_CPV_B__KSTAR_L_L, "P'_8_CPV_B__K*_l_l"},
        {DGAMMA_DQ2_BS__PHI_L_L, "dGamma_dq2_Bs__phi_l_l"},
        {DGAMMA_BAR_DQ2_BS__PHI_L_L, "dGamma_bar_dq2_Bs__phi_l_l"},
        {A_FB_CPV_BS__PHI_L_L, "A_FB_CPV_Bs__phi_l_l"},
        {F_L_BS_PHI_L_L, "F_L_Bs__phi_l_l"},
        {A_T_2_BS_PHI_L_L, "A_T_2_Bs__phi_l_l"},
        {A_T_RE_CPV_BS_PHI_L_L, "A_T_Re_CPV_Bs__phi_l_l"},
        {A_T_IM_CPV_BS_PHI_L_L, "A_T_Im_CPV_Bs__phi_l_l"},
        {P_PRIME_4_BS_PHI_L_L, "P'_4_Bs__phi_l_l"},
        {P_PRIME_6_BS_PHI_L_L, "P'_6_Bs__phi_l_l"},
        {S_2S_BS_PHI_L_L, "S_2s_Bs__phi_l_l"},
        {S_3_BS_PHI_L_L, "S_3_Bs__phi_l_l"},
        {S_4_BS_PHI_L_L, "S_4_Bs__phi_l_l"},
        {S_7_BS_PHI_L_L, "S_7_Bs__phi_l_l"},
        {A_5_BS_PHI_L_L, "A_5_Bs__phi_l_l"},
        {A_6C_BS_PHI_L_L, "A_6c_Bs__phi_l_l"},
        {A_8_BS_PHI_L_L, "A_8_Bs__phi_l_l"},
        {A_9_BS_PHI_L_L, "A_9_Bs__phi_l_l"},
        {P_2_CPV_BS_PHI_L_L, "P_2_CPV_Bs__phi_l_l"},
        {P_3_CPV_BS_PHI_L_L, "P_3_CPV_Bs__phi_l_l"},
        {P_PRIME_5_CPV_BS_PHI_L_L, "P'_5_CPV_Bs__phi_l_l"},
        {P_PRIME_8_CPV_BS_PHI_L_L, "P'_8_CPV_Bs__phi_l_l"},
        {Q_8M_BS_PHI_L_L, "Q_8m_Bs__phi_l_l"},
        {Q_8P_BS_PHI_L_L, "Q_8p_CPV_Bs__phi_l_l"},
        {Q_9_BS_PHI_L_L, "Q_9_CPV_Bs__phi_l_l"},
        {DGAMMA_DQ2_B__K_L_L, "dGamma_dq2_B__K_l_l"},
        {A_FB_B__K_L_L, "A_FB_B__K_l_l"},
        {F_H_B__K_L_L, "F_H_B__K_l_l"},
        {DGAMMA_DQ2_CP_AVG_LAMBDA_B__LAMBDA_L_L, "dG_dq2_CPA_Lambda_b__Lambda_l_l"},
        {A_FB_L_LAMBDA_B__LAMBDA_L_L, "A_FB_l_Lambda_b__Lambda_l_l"},
        {A_FB_H_LAMBDA_B__LAMBDA_L_L, "A_FB_h_Lambda_b__Lambda_l_l"},
        {A_FB_LH_LAMBDA_B__LAMBDA_L_L, "A_FB_lh_Lambda_b__Lambda_l_l"},
        {F_L_LAMBDA_B__LAMBDA_L_L, "F_L_Lambda_b__Lambda_l_l"},
        {F_T_LAMBDA_B__LAMBDA_L_L, "F_T_Lambda_b__Lambda_l_l"},
        {PHI_D, "phi_d"},
        {DELTA_M_BD, "Delta_M_Bd"},
        {PHI_S, "phi_s"},
        {DELTA_M_BS, "Delta_M_Bs"},
        {A_FS, "a_fs"},
        {DELTA_M_K, "Delta_M_K"},
        {ABS_EPSILON_K, "|epsilon_K|"},
        {X_D, "x_D"},
        {BR_KL__MU_MU, "BR_K_L__mu_mu"},
        {BR_KS__MU_MU, "BR_K_S__mu_mu"},
        {BR_K__PI_NU_NU, "BR_K__pi_nu_nu"},
        {BR_KL__PI0_NU_NU, "BR_KL__pi0_nu_nu"},
        {BR_K__MU_NU__BR_PI__MU_NU, "BR_K__mu_nu / BR_pi__mu_nu"},
        {R_MU23, "R_mu23"},
        {BR_D__MU_NU, "BR_D__mu_nu"},
        {BR_DS__MU_NU, "BR_Ds__mu_nu"},
        {BR_DS__TAU_NU, "BR_Ds__tau_nu"},
        {TEST, "test"},
    };
    return m;
}

const std::map<Observables, LhaID>& observable_flha_mapping() {
    static const std::map<Observables, LhaID> m = {
        {BR_BS_MUMU,        LhaID(531, 1, 2, 13, -13)},
        {BR_BS_MUMU_UNTAG,  LhaID(531, 15, 2, 13, -13)},
        {BR_BD_MUMU,        LhaID(511, 1, 2, 13, -13)},
        {BR_BU_TAU_NU,      LhaID(521, 1, 2, -15, 16)},
        {R_TAU_NU,          LhaID(521, 2, 2, -15, 16)},
        {ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, LhaID(521, 4, 2, 313, 22)},
        {BR_B_XS_GAMMA,     LhaID(5,   1, 2, 3, 22)},
        {BR_B__D_TAU_NU,    LhaID(521, 1, 3, 421, -15, 16)},
        {A_FB_B__D_TAU_NU,  LhaID(521, 5, 3, 421, -15, 16)},
        {P_TAU_B__D_TAU_NU, LhaID(521, 82,3, 421, -15, 16)},
        {R_D,               LhaID(521, 11,3, 421, -15, 16)},
        {BR_B__DSTAR_TAU_NU,LhaID(521, 1, 3, 423, -15, 16)},
        {A_FB_B__DSTAR_TAU_NU, LhaID(521, 5, 3, 423, -15, 16)},
        {P_TAU_B__DSTAR_TAU_NU, LhaID(521, 82,3, 423, -15, 16)},
        {P_D_B__DSTAR_TAU_NU,   LhaID(521, 81,3, 423, -15, 16)},
        {R_DSTAR,           LhaID(521, 11,3, 423, -15, 16)},
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
        {B, "BCoefficients"},
        {BPrime, "BPrimeCoefficients"},
        {BScalar, "BScalarCoefficients"},
        {CC_bc, "BcChargedCurrentCoefficients"},
        {CC_bu, "BuChargedCurrentCoefficients"},
        {CC_cs, "DsChargedCurrentCoefficients"},
        {CC_cd, "DdChargedCurrentCoefficients"},
        {CC_su, "KuChargedCurrentCoefficients"},
        {CC_du, "PIuChargedCurrentCoefficients"},
        {MESON_MIXING, "MesonMixing"},
        {CUSTOM_GROUP, "Custom_Group"},
        {K, "K"},
    };
    return m;
}



const std::map<WCoef, std::string>& wcoef_mapping() {
    static const std::map<WCoef, std::string> m = {
        {C1, "C1"}, {C2, "C2"}, {C3, "C3"}, {C4, "C4"}, {C5, "C5"}, {C6, "C6"},
        {C7, "C7"}, {C8, "C8"}, {C9, "C9"}, {C10, "C10"},
        {CP1, "CP1"}, {CP2, "CP2"}, {CP3, "CP3"}, {CP4, "CP4"},
        {CP5, "CP5"}, {CP6, "CP6"}, {CP7, "CP7"}, {CP8, "CP8"},
        {CP9, "CP9"}, {CP10, "CP10"}, {CQ1, "CQ1"}, {CQ2, "CQ2"},
        {CPQ1, "CPQ1"}, {CPQ2, "CPQ2"},
        {C_V1_bc, "C_V1_bc"}, {C_V2_bc, "C_V2_bc"}, {C_S1_bc, "C_S1_bc"}, {C_S2_bc, "C_S2_bc"}, {C_T_bc, "C_T_bc"},
        {C_V1_bu, "C_V1_bu"}, {C_V2_bu, "C_V2_bu"}, {C_S1_bu, "C_S1_bu"}, {C_S2_bu, "C_S2_bu"}, {C_T_bu, "C_T_bu"},
        {C_V1_cs, "C_V1_cs"}, {C_V2_cs, "C_V2_cs"}, {C_S1_cs, "C_S1_cs"}, {C_S2_cs, "C_S2_cs"}, {C_T_cs, "C_T_cs"},
        {C_V1_cd, "C_V1_cd"}, {C_V2_cd, "C_V2_cd"}, {C_S1_cd, "C_S1_cd"}, {C_S2_cd, "C_S2_cd"}, {C_T_cd, "C_T_cd"},
        {C_V1_su, "C_V1_su"}, {C_V2_su, "C_V2_su"}, {C_S1_su, "C_S1_su"}, {C_S2_su, "C_S2_su"}, {C_T_su, "C_T_su"},
        {C_V1_du, "C_V1_du"}, {C_V2_du, "C_V2_du"}, {C_S1_du, "C_S1_du"}, {C_S2_du, "C_S2_du"}, {C_T_du, "C_T_du"},
        {C_BD_1, "C_BD_1"}, {CT_BD_1, "CT_BD_1"}, {C_BD_2, "C_BD_2"}, {CT_BD_2, "CT_BD_2"}, {C_BD_3, "C_BD_3"}, {CT_BD_3, "CT_BD_3"}, {C_BD_4, "C_BD_4"}, {C_BD_5, "C_BD_5"},
        {C_BS_1, "C_BS_1"}, {CT_BS_1, "CT_BS_1"}, {C_BS_2, "C_BS_2"}, {CT_BS_2, "CT_BS_2"}, {C_BS_3, "C_BS_3"}, {CT_BS_3, "CT_BS_3"}, {C_BS_4, "C_BS_4"}, {C_BS_5, "C_BS_5"},
        {C_SD_1, "C_SD_1"}, {CT_SD_1, "CT_SD_1"}, {C_SD_2, "C_SD_2"}, {CT_SD_2, "CT_SD_2"}, {C_SD_3, "C_SD_3"}, {CT_SD_3, "CT_SD_3"}, {C_SD_4, "C_SD_4"}, {C_SD_5, "C_SD_5"},
        {C_CU_1, "C_CU_1"}, {CT_CU_1, "CT_CU_1"}, {C_CU_2, "C_CU_2"}, {CT_CU_2, "CT_CU_2"}, {C_CU_3, "C_CU_3"}, {CT_CU_3, "CT_CU_3"}, {C_CU_4, "C_CU_4"}, {C_CU_5, "C_CU_5"},
        {CK9, "CK9"}, {CPK9, "CPK9"}, {CK10, "CK10"}, {CPK10, "CPK10"}, {CKQ1, "CKQ1"}, {CKQ2, "CKQ2"}, {CPKQ1, "CPKQ1"}, {CPKQ2, "CPKQ2"}, {CK_L, "CK_L"}
    };
    return m;
}


const std::map<WCoef, std::pair<int, int>>& wcoef_flha_mapping() {
    static const std::map<WCoef, std::pair<int, int>> m = {
        {C1, {3040405, 6161}}, {C2, {3040405, 4141}}, {C3, {3050707, 4133}},
        {C4, {3050707, 6153}}, {C5, {3050707, 4536}}, {C6, {3050707, 6556}},
        {C7, {305, 4422}}, {C8, {305, 6421}}, {C9, {3051313, 4133}}, {C10, {3051313, 4137}},
        {CP1, {3040405, 6262}}, {CP2, {3040405, 4242}}, {CP3, {3050707, 4233}},
        {CP4, {3050707, 6253}}, {CP5, {3050707, 4636}}, {CP6, {3050707, 6656}},
        {CP7, {305, 4322}}, {CP8, {305, 4321}}, {CP9, {3051313, 4233}},
        {CP10, {3051313, 4234}}, {CQ1, {3051313, 3230}}, {CQ2, {3051313, 3233}},
        {CPQ1, {3051313, 3130}}, {CPQ2, {3051313, 3133}},
        {C_V1_bc, {4051516, 4141}}, {C_V2_bc, {4051516, 4241}},
        {C_S1_bc, {4051516, 3231}}, {C_S2_bc, {4051516, 3131}}, {C_T_bc, {4051516, 4343}},
        {C_V1_bu, {4051516, 1}}, {C_V2_bu, {4051516, 2}},
        {C_S1_bu, {4051516, 3}}, {C_S2_bu, {4051516, 4}}, {C_T_bu, {4051516, 5}},
        {C_V1_cs, {4051516, 11}}, {C_V2_cs, {4051516, 22}},
        {C_S1_cs, {4051516, 33}}, {C_S2_cs, {4051516, 44}}, {C_T_cs, {4051516, 55}},
        {C_V1_cd, {4051516, 111}}, {C_V2_cd, {4051516, 222}},
        {C_S1_cd, {4051516, 333}}, {C_S2_cd, {4051516, 444}}, {C_T_cd, {4051516, 555}},
        {C_V1_su, {4051516, 1111}}, {C_V2_su, {4051516, 2222}},
        {C_S1_su, {4051516, 3333}}, {C_S2_su, {4051516, 4444}}, {C_T_su, {4051516, 5555}},
        {C_V1_du, {4051516, 11111}}, {C_V2_du, {4051516, 22222}},
        {C_S1_du, {4051516, 33333}}, {C_S2_du, {4051516, 44444}}, {C_T_du, {4051516, 55555}},
        {C_BD_1, {1050105, 4141}}, {CT_BD_1, {1050105, 4242}}, {C_BD_2, {1050105, 3131}}, {CT_BD_2, {1050105, 3232}},
        {C_BD_3, {1050105, 7171}}, {CT_BD_3, {1050105, 7272}}, {C_BD_4, {1050105, 3132}}, {C_BD_5, {1050105, 7172}},
        {C_BS_1, {3050305, 4141}}, {CT_BS_1, {3050305, 4242}}, {C_BS_2, {3050305, 3131}}, {CT_BS_2, {3050305, 3232}},
        {C_BS_3, {3050305, 7171}}, {CT_BS_3, {3050305, 7272}}, {C_BS_4, {3050305, 3132}}, {C_BS_5, {3050305, 7172}},
        {C_SD_1, {1030103, 4141}}, {CT_SD_1, {1030103, 4242}}, {C_SD_2, {1030103, 3131}}, {CT_SD_2, {1030103, 3232}},
        {C_SD_3, {1030103, 7171}}, {CT_SD_3, {1030103, 7272}}, {C_SD_4, {1030103, 3132}}, {C_SD_5, {1030103, 7172}},
        {C_CU_1, {2040204, 4141}}, {CT_CU_1, {2040204, 4242}}, {C_CU_2, {2040204, 3131}}, {CT_CU_2, {2040204, 3232}},
        {C_CU_3, {2040204, 7171}}, {CT_CU_3, {2040204, 7272}}, {C_CU_4, {2040204, 3132}}, {C_CU_5, {2040204, 7172}},

        {CK9, {0, 1}}, {CPK9, {0, 2}}, {CK10, {0, 3}}, {CPK10, {0, 4}}, //TODO : real coefs
        {CKQ1, {0, 5}}, {CKQ2, {0, 6}}, {CPKQ1, {0, 7}}, {CPKQ2, {0, 8}},
        {CK_L, {0, 9}},
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
        {Model::MARTY, "MARTY"},
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
        {B__D_l_nu, "B__D_l_nu"},
        {B__Dstar_l_nu, "B__Dstar_l_nu"},
        {B__Kstar_gamma, "B__Kstar_gamma"},
        {B__l_l, "B__l_l"},
        {B__l_nu, "B__l_nu"},
        {B__Xs, "B__Xs"},
        {B__Xs_l_l, "B__Xs_ll"},
        {B__Kstar_l_l, "B__K*_l_l"},
        {B__K_l_l, "B__K_l_l"},
        {Bs__phi_l_l, "Bs__phi_l_l"},
        {Lambda_b__Lambda_l_l, "Lambda_b__Lambda_l_l"},
        {M0_Mix, "M0_Mix"},
        {K__l_l, "K__l_l"},
        {K__pi_nu_nu, "K__pi_nu_nu"},
        {K__l_nu, "K__l_nu"},
        {D__l_nu, "D__l_nu"},
        {Ds__l_nu, "Ds__l_nu"},
    };
    return m;
}

const std::map<Decays, std::vector<Observables>>& decay_observable_mapping() {
    static const std::map<Decays, std::vector<Observables>> m = {
        {B__D_l_nu,             {BR_B__D_TAU_NU, A_FB_B__D_TAU_NU, P_TAU_B__D_TAU_NU, R_D}},
        {B__Dstar_l_nu,         {BR_B__DSTAR_TAU_NU, A_FB_B__DSTAR_TAU_NU, P_TAU_B__DSTAR_TAU_NU, P_D_B__DSTAR_TAU_NU, R_DSTAR}},
        {B__Kstar_gamma,        {BR_B__KSTAR_GAMMA, ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA}},
        {B__l_l,                {BR_BS_MUMU, BR_BD_MUMU, BR_BS_MUMU_UNTAG}},
        {B__l_nu,               {BR_BU_TAU_NU, R_TAU_NU}},
        {B__Xs,                 {BR_B_XS_GAMMA}},
        {B__Xs_l_l,             {BR_B__Xs_mu_mu__LOW_Q2, BR_B__Xs_mu_mu__HIGH_Q2, BR_B__Xs_tau_tau__HIGH_Q2}},
        {B__Kstar_l_l,          {DGAMMA_DQ2_B__KSTAR_L_L, A_FB_B__KSTAR_L_L, Q0_A_FB_B__KSTAR_L_L, A_CP_B__KSTAR_L_L, F_L_B__KSTAR_L_L, F_T_B__KSTAR_L_L, A_T_1_B__KSTAR_L_L, A_T_2_B__KSTAR_L_L, A_T_3_B__KSTAR_L_L, A_T_4_B__KSTAR_L_L, A_T_5_B__KSTAR_L_L, A_T_RE_B__KSTAR_L_L, A_T_RE_CPV_B__KSTAR_L_L, A_IM_B__KSTAR_L_L, ALPHA_K_B__KSTAR_L_L, H_T_1_B__KSTAR_L_L, H_T_2_B__KSTAR_L_L, H_T_3_B__KSTAR_L_L, P_2_B__KSTAR_L_L, P_3_B__KSTAR_L_L, P_6_B__KSTAR_L_L, P_8_B__KSTAR_L_L, P_PRIME_4_B__KSTAR_L_L, P_PRIME_5_B__KSTAR_L_L, P_PRIME_6_B__KSTAR_L_L, P_PRIME_8_B__KSTAR_L_L, S_3_B__KSTAR_L_L, S_4_B__KSTAR_L_L, S_5_B__KSTAR_L_L, S_6C_B__KSTAR_L_L, S_7_B__KSTAR_L_L, S_8_B__KSTAR_L_L, S_9_B__KSTAR_L_L, A_3_B__KSTAR_L_L, A_4_B__KSTAR_L_L, A_5_B__KSTAR_L_L, A_6S_B__KSTAR_L_L, A_7_B__KSTAR_L_L, A_8_B__KSTAR_L_L, A_9_B__KSTAR_L_L, P_1_CPV_B__KSTAR_L_L, P_2_CPV_B__KSTAR_L_L, P_3_CPV_B__KSTAR_L_L, P_PRIME_4_CPV_B__KSTAR_L_L, P_PRIME_5_CPV_B__KSTAR_L_L, P_PRIME_6_CPV_B__KSTAR_L_L, P_PRIME_8_CPV_B__KSTAR_L_L}},
        {B__K_l_l,              {DGAMMA_DQ2_B__K_L_L, A_FB_B__K_L_L, F_H_B__K_L_L}},
        {Bs__phi_l_l,           {DGAMMA_DQ2_BS__PHI_L_L, DGAMMA_BAR_DQ2_BS__PHI_L_L, A_FB_CPV_BS__PHI_L_L, F_L_BS_PHI_L_L, A_T_2_BS_PHI_L_L, A_T_RE_CPV_BS_PHI_L_L, A_T_IM_CPV_BS_PHI_L_L, P_PRIME_4_BS_PHI_L_L, P_PRIME_6_BS_PHI_L_L, S_2S_BS_PHI_L_L, S_3_BS_PHI_L_L, S_4_BS_PHI_L_L, S_7_BS_PHI_L_L, A_5_BS_PHI_L_L, A_6C_BS_PHI_L_L, A_8_BS_PHI_L_L, A_9_BS_PHI_L_L, P_2_CPV_BS_PHI_L_L, P_3_CPV_BS_PHI_L_L, P_PRIME_5_CPV_BS_PHI_L_L, P_PRIME_8_CPV_BS_PHI_L_L, Q_8M_BS_PHI_L_L, Q_8P_BS_PHI_L_L, Q_9_BS_PHI_L_L}},
        {Lambda_b__Lambda_l_l,  {DGAMMA_DQ2_CP_AVG_LAMBDA_B__LAMBDA_L_L, A_FB_L_LAMBDA_B__LAMBDA_L_L, A_FB_H_LAMBDA_B__LAMBDA_L_L, A_FB_LH_LAMBDA_B__LAMBDA_L_L, F_L_LAMBDA_B__LAMBDA_L_L, F_T_LAMBDA_B__LAMBDA_L_L}},
        {M0_Mix,                {PHI_D, DELTA_M_BD, PHI_S, DELTA_M_BS, A_FS, DELTA_M_K, ABS_EPSILON_K, X_D}},
        {K__l_l,                {BR_KL__MU_MU, BR_KS__MU_MU}},
        {K__pi_nu_nu,           {BR_K__PI_NU_NU, BR_KL__PI0_NU_NU}},
        {K__l_nu,               {BR_K__MU_NU__BR_PI__MU_NU, R_MU23}},
        {D__l_nu,               {BR_D__MU_NU}},
        {Ds__l_nu,              {BR_DS__MU_NU, BR_DS__TAU_NU}},
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

