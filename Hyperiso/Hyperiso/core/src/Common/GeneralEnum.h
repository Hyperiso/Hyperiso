#ifndef GENERAL_ENUM_H
#define GENERAL_ENUM_H

enum class Observables {
    /* For testing purposes */
    TEST,
    /* Rare B decays */
    // B > l l
    BR_BS_MUMU,
    BR_BS_MUMU_UNTAG,
    BR_BD_MUMU,
    // B > l nu
    R_TAU_NU,
    BR_BU_TAU_NU,
    // B > K* gamma
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA,
    BR_B__KSTAR_GAMMA,
    // b > s gamma
    BR_B_XS_GAMMA,
    // B > D l nu
    BR_B__D_TAU_NU,
    A_FB_B__D_TAU_NU,
    P_TAU_B__D_TAU_NU,
    R_D,
    // B > D* l nu
    BR_B__DSTAR_TAU_NU,
    A_FB_B__DSTAR_TAU_NU,
    P_TAU_B__DSTAR_TAU_NU,
    P_D_B__DSTAR_TAU_NU,
    R_DSTAR,
    // b > s l l 
    BR_B__Xs_mu_mu__LOW_Q2,
    BR_B__Xs_mu_mu__HIGH_Q2,
    BR_B__Xs_tau_tau__HIGH_Q2,
    // B > K* l l 
    DGAMMA_DQ2_B__KSTAR_L_L,
    DGAMMA_BAR_DQ2_B__KSTAR_L_L,
    A_FB_B__KSTAR_L_L,
    Q0_A_FB_B__KSTAR_L_L,
    A_CP_B__KSTAR_L_L,
    F_L_B__KSTAR_L_L,
    F_T_B__KSTAR_L_L,
    A_T_1_B__KSTAR_L_L,
    A_T_2_B__KSTAR_L_L,
    A_T_3_B__KSTAR_L_L,
    A_T_4_B__KSTAR_L_L,
    A_T_5_B__KSTAR_L_L,
    A_T_RE_B__KSTAR_L_L,
    A_T_RE_CPV_B__KSTAR_L_L,
    A_IM_B__KSTAR_L_L,
    ALPHA_K_B__KSTAR_L_L,
    H_T_1_B__KSTAR_L_L,
    H_T_2_B__KSTAR_L_L,
    H_T_3_B__KSTAR_L_L,
    P_2_B__KSTAR_L_L,
    P_3_B__KSTAR_L_L,
    P_6_B__KSTAR_L_L,
    P_8_B__KSTAR_L_L,
    P_PRIME_4_B__KSTAR_L_L,
    P_PRIME_5_B__KSTAR_L_L,
    P_PRIME_6_B__KSTAR_L_L,
    P_PRIME_8_B__KSTAR_L_L,
    S_3_B__KSTAR_L_L,
    S_4_B__KSTAR_L_L,
    S_5_B__KSTAR_L_L,
    S_6C_B__KSTAR_L_L,
    S_7_B__KSTAR_L_L,
    S_8_B__KSTAR_L_L,
    S_9_B__KSTAR_L_L,
    A_3_B__KSTAR_L_L,
    A_4_B__KSTAR_L_L,
    A_5_B__KSTAR_L_L,
    A_6S_B__KSTAR_L_L,
    A_7_B__KSTAR_L_L,
    A_8_B__KSTAR_L_L,
    A_9_B__KSTAR_L_L,
    P_1_CPV_B__KSTAR_L_L,
    P_2_CPV_B__KSTAR_L_L,
    P_3_CPV_B__KSTAR_L_L,
    P_PRIME_4_CPV_B__KSTAR_L_L,
    P_PRIME_5_CPV_B__KSTAR_L_L,
    P_PRIME_6_CPV_B__KSTAR_L_L,
    P_PRIME_8_CPV_B__KSTAR_L_L,
    // Bs > phi l l
    DGAMMA_DQ2_BS__PHI_L_L,
    DGAMMA_BAR_DQ2_BS__PHI_L_L,
    A_FB_CPV_BS__PHI_L_L,
    F_L_BS_PHI_L_L,
    A_T_2_BS_PHI_L_L,
    A_T_RE_CPV_BS_PHI_L_L,
    A_T_IM_CPV_BS_PHI_L_L,
    P_PRIME_4_BS_PHI_L_L,
    P_PRIME_6_BS_PHI_L_L,
    S_2S_BS_PHI_L_L,
    S_3_BS_PHI_L_L,
    S_4_BS_PHI_L_L,
    S_7_BS_PHI_L_L,
    A_5_BS_PHI_L_L,
    A_6C_BS_PHI_L_L,
    A_8_BS_PHI_L_L,
    A_9_BS_PHI_L_L,
    P_2_CPV_BS_PHI_L_L,
    P_3_CPV_BS_PHI_L_L,
    P_PRIME_5_CPV_BS_PHI_L_L,
    P_PRIME_8_CPV_BS_PHI_L_L,
    Q_8M_BS_PHI_L_L,
    Q_8P_BS_PHI_L_L,
    Q_9_BS_PHI_L_L,
    // B > K l l
    DGAMMA_DQ2_B__K_L_L,
    A_FB_B__K_L_L,
    F_H_B__K_L_L,
    // Lambda_b > Lambda l l
    DGAMMA_DQ2_CP_AVG_LAMBDA_B__LAMBDA_L_L,
    A_FB_L_LAMBDA_B__LAMBDA_L_L,
    A_FB_H_LAMBDA_B__LAMBDA_L_L,
    A_FB_LH_LAMBDA_B__LAMBDA_L_L,
    F_L_LAMBDA_B__LAMBDA_L_L,
    F_T_LAMBDA_B__LAMBDA_L_L,
    /* Neutral Meson Mixing */
    // Bd
    PHI_D,
    DELTA_M_BD,
    // Bs
    PHI_S,
    DELTA_M_BS,
    A_FS,
    // K
    DELTA_M_K,
    ABS_EPSILON_K,
    // D
    X_D,
    /* Rare K decays */
    // K_L,S > mu mu
    BR_KL__MU_MU,
    BR_KS__MU_MU,
    // K > mu nu
    BR_K__MU_NU__BR_PI__MU_NU,
    R_MU23,
    // K > pi nu nu
    BR_K__PI_NU_NU,
    BR_KL__PI0_NU_NU,
    /* Rare D decays */
    BR_D__MU_NU,
    BR_DS__MU_NU,
    BR_DS__TAU_NU,
};

enum class Decays {
    B__D_l_nu,
    B__Dstar_l_nu,
    B__Kstar_gamma,
    B__l_l,
    B__l_nu,
    B__Xs,
    B__Xs_l_l,
    B__K_l_l,
    B__Kstar_l_l,
    Bs__phi_l_l,
    Lambda_b__Lambda_l_l,
    M0_Mix,
    K__l_l,
    K__pi_nu_nu,
    K__l_nu,
    D__l_nu,
    Ds__l_nu,
};

enum class QCDOrder {
    NONE,
    LO,
    NLO,
    NNLO
};

enum class WCoef {
    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10,                                // b > s l l
    CQ1, CQ2,                                                               // b > s l l
    CP1, CP2, CP3, CP4, CP5, CP6, CP7, CP8, CP9, CP10, CPQ1, CPQ2,          // b > s l l
    C_V1_bc, C_V2_bc, C_S1_bc, C_S2_bc, C_T_bc,                             // b > c l nu
    C_V1_bu, C_V2_bu, C_S1_bu, C_S2_bu, C_T_bu,                             // b > u l nu
    C_V1_cs, C_V2_cs, C_S1_cs, C_S2_cs, C_T_cs,                             // c > s l nu
    C_V1_cd, C_V2_cd, C_S1_cd, C_S2_cd, C_T_cd,                             // c > d l nu
    C_V1_su, C_V2_su, C_S1_su, C_S2_su, C_T_su,                             // s > u l nu
    C_V1_du, C_V2_du, C_S1_du, C_S2_du, C_T_du,                             // d > u l nu
    C_BD_1, CT_BD_1, C_BD_2, CT_BD_2, C_BD_3, CT_BD_3, C_BD_4, C_BD_5,      // Bd0 mixing
    C_BS_1, CT_BS_1, C_BS_2, CT_BS_2, C_BS_3, CT_BS_3, C_BS_4, C_BS_5,      // Bs0 mixing
    C_SD_1, CT_SD_1, C_SD_2, CT_SD_2, C_SD_3, CT_SD_3, C_SD_4, C_SD_5,      // K0 mixing
    C_CU_1, CT_CU_1, C_CU_2, CT_CU_2, C_CU_3, CT_CU_3, C_CU_4, C_CU_5,      // D0 mixing
    CK9, CPK9, CK10, CPK10, CKQ1, CKQ2, CPKQ1, CPKQ2, CK_L                        // s > d l l
};

enum class WGroup {
    B, 
    BPrime, 
    BScalar,
    CC_bc,
    CC_bu,
    CC_cs,
    CC_cd,
    CC_su,
    CC_du,
    MESON_MIXING,
    K
};


enum class ParameterType {
    SM,
    BSM,
    FLAVOR,
    WILSON,
    DECAY,
    PASSTHROUGH,
    OBSERVABLE,
};

enum class Model {
    SM,
    SUSY,
    THDM,
    MARTY
};

enum class WilsonBasis {
    B_STANDARD, 
    B_TRADITIONAL
};

enum class ContributionType {
    SM, 
    BSM,
    TOTAL
};

enum class MassType {
    POLE,
    MSBAR
};

enum class ScaleType {
    MATCHING,
    HADRONIC
};

enum class DataType { 
    VALUE,
    STD_STAT,
    STD_SYST,
    STD_COMBINED
};

enum class UncertaintyType {
    STAT,
    SYST,
    COMBINED
};


#endif