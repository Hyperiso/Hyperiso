#ifndef GENERAL_ENUM_H
#define GENERAL_ENUM_H

enum class Observables {
    BR_BS_MUMU,
    BR_BS_MUMU_UNTAG,
    BR_BD_MUMU,
    R_TAU_NU,
    BR_BU_TAU_NU,
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA,
    BR_B_XS_GAMMA,
    BR_B__D_TAU_NU,
    A_FB_B__D_TAU_NU,
    P_TAU_B__D_TAU_NU,
    R_D,
    BR_B__DSTAR_TAU_NU,
    A_FB_B__DSTAR_TAU_NU,
    P_TAU_B__DSTAR_TAU_NU,
    P_D_B__DSTAR_TAU_NU,
    R_DSTAR,
};

enum class Decays {
    B__D_l_nu,
    B__Dstar_l_nu,
    B__Kstar,
    B__l_l,
    B__l_nu,
    B__Xs,
};

enum class QCDOrder {
    NONE,
    LO,
    NLO,
    NNLO
};

enum class WCoef {
    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, CQ1, CQ2, CP1, CP2, CP3, CP4, CP5, CP6, CP7, CP8, CP9, CP10, CPQ1, CPQ2, C_V1, C_V2, C_S1, C_S2, C_T
};

enum class WGroup {
    B, 
    BPrime, 
    BScalar,
    BCC,
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
    CUSTOM
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