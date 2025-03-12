#include "General.h"

const std::map<Decays, std::vector<Observables>> DecayMapper::obs_mapping = {
    {Decays::B__D_l_nu,     {Observables::BR_B__D_TAU_NU, Observables::A_FB_B__D_TAU_NU, Observables::P_TAU_B__D_TAU_NU, Observables::R_D}},
    {Decays::B__Dstar_l_nu, {Observables::BR_B__DSTAR_TAU_NU, Observables::A_FB_B__DSTAR_TAU_NU, Observables::P_TAU_B__DSTAR_TAU_NU, Observables::P_D_B__DSTAR_TAU_NU, Observables::R_DSTAR}},
    {Decays::B__Kstar,      {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA}},
    {Decays::B__l_l,        {Observables::BR_BS_MUMU, Observables::BR_BD_MUMU, Observables::BR_BS_MUMU_UNTAG}},
    {Decays::B__l_nu,       {Observables::BR_BU_TAU_NU, Observables::R_TAU_NU}},
    {Decays::B__Xs,         {Observables::BR_B_XS_GAMMA}},
};

const std::map<QCDOrder, std::string> OrderMapper::mapping = {
    {QCDOrder::NONE, "None"}, 
    {QCDOrder::LO, "LO"}, 
    {QCDOrder::NLO, "NLO"}, 
    {QCDOrder::NNLO, "NNLO"}
};

const std::map<std::string, QCDOrder> OrderMapper::inverse_mapping = {
    {"None", QCDOrder::NONE}, 
    {"LO", QCDOrder::LO}, 
    {"NLO", QCDOrder::NLO},
    {"NNLO", QCDOrder::NNLO}
};

const std::map<WCoef, std::string> WCoefMapper::mapping = {
    {WCoef::C1, "C1"},
    {WCoef::C2, "C2"},
    {WCoef::C3, "C3"},
    {WCoef::C4, "C4"},
    {WCoef::C5, "C5"},
    {WCoef::C6, "C6"},
    {WCoef::C7, "C7"},
    {WCoef::C8, "C8"},
    {WCoef::C9, "C9"},
    {WCoef::C10, "C10"},
    {WCoef::CP1, "CP1"},
    {WCoef::CP2, "CP2"},
    {WCoef::CP3, "CP3"},
    {WCoef::CP4, "CP4"},
    {WCoef::CP5, "CP5"},
    {WCoef::CP6, "CP6"},
    {WCoef::CP7, "CP7"},
    {WCoef::CP8, "CP8"},
    {WCoef::CP9, "CP9"},
    {WCoef::CP10, "CP10"},
    {WCoef::CQ1, "CQ1"},
    {WCoef::CQ2, "CQ2"},
    {WCoef::CPQ1, "CPQ1"},
    {WCoef::CPQ2, "CPQ2"}, 
    {WCoef::CBlnu_A, "C_Blnu_A"}, 
    {WCoef::CBlnu_P, "C_Blnu_P"}, 
    {WCoef::C_V1, "C_V1"}, 
    {WCoef::C_V2, "C_V2"}, 
    {WCoef::C_S1, "C_S1"}, 
    {WCoef::C_S2, "C_S2"}, 
    {WCoef::C_T, "C_T"}, 
};

const std::map<WCoef, std::pair<int, int>> WCoefMapper::flha_mapping = {
    {WCoef::C1, {3040405, 6161}},
    {WCoef::C2, {3040405, 4141}},
    {WCoef::C3, {3050707, 4133}},
    {WCoef::C4, {3050707, 6153}},
    {WCoef::C5, {3050707, 4536}},
    {WCoef::C6, {3050707, 6556}},
    {WCoef::C7, {305, 4422}},
    {WCoef::C8, {305, 6421}},
    {WCoef::C9, {3051313, 4133}},
    {WCoef::C10, {3051313, 4137}},
    {WCoef::CP1, {3040405, 6262}},
    {WCoef::CP2, {3040405, 4242}},
    {WCoef::CP3, {3050707, 4233}},
    {WCoef::CP4, {3050707, 6253}},
    {WCoef::CP5, {3050707, 4636}},
    {WCoef::CP6, {3050707, 6656}},
    {WCoef::CP7, {305, 4322}},
    {WCoef::CP8, {305, 4321}},
    {WCoef::CP9, {3051313, 4233}},
    {WCoef::CP10, {3051313, 4234}},
    {WCoef::CQ1, {3051313, 3230}},
    {WCoef::CQ2, {3051313, 3233}},
    {WCoef::CPQ1, {3051313, 3130}},
    {WCoef::CPQ2, {3051313, 3133}}, 
};

const std::map<std::pair<int, int>, WCoef> WCoefMapper::inverse_flha_mapping = {
    {{3040405, 6161}, WCoef::C1},
    {{3040405, 4141}, WCoef::C2},
    {{3050707, 4133}, WCoef::C3},
    {{3050707, 6153}, WCoef::C4},
    {{3050707, 4536}, WCoef::C5},
    {{3050707, 6556}, WCoef::C6},
    {{305, 4422}, WCoef::C7},
    {{305, 6421}, WCoef::C8},
    {{3051313, 4133}, WCoef::C9},
    {{3051313, 4137}, WCoef::C10},
    {{3040405, 6262}, WCoef::CP1},
    {{3040405, 4242}, WCoef::CP2},
    {{3050707, 4233}, WCoef::CP3},
    {{3050707, 6253}, WCoef::CP4},
    {{3050707, 4636}, WCoef::CP5},
    {{3050707, 6656}, WCoef::CP6},
    {{305, 4322}, WCoef::CP7},
    {{305, 4321}, WCoef::CP8},
    {{3051313, 4233}, WCoef::CP9},
    {{3051313, 4234}, WCoef::CP10},
    {{3051313, 3230}, WCoef::CQ1},
    {{3051313, 3233}, WCoef::CQ2},
    {{3051313, 3130}, WCoef::CPQ1},
    {{3051313, 3133}, WCoef::CPQ2}, 
};

const std::map<std::string, WCoef> WCoefMapper::inverse_mapping = {
    {"C1", WCoef::C1},
    {"C2", WCoef::C2},
    {"C3", WCoef::C3},
    {"C4", WCoef::C4},
    {"C5", WCoef::C5},
    {"C6", WCoef::C6},
    {"C7", WCoef::C7},
    {"C8", WCoef::C8},
    {"C9", WCoef::C9},
    {"C10", WCoef::C10},
    {"CP1", WCoef::CP1},
    {"CP2", WCoef::CP2},
    {"CP3", WCoef::CP3},
    {"CP4", WCoef::CP4},
    {"CP5", WCoef::CP5},
    {"CP6", WCoef::CP6},
    {"CP7", WCoef::CP7},
    {"CP8", WCoef::CP8},
    {"CP9", WCoef::CP9},
    {"CP10", WCoef::CP10},
    {"CQ1", WCoef::CQ1},
    {"CQ2", WCoef::CQ2},
    {"CPQ1", WCoef::CPQ1},
    {"CPQ2", WCoef::CPQ2},
    {"C_Blnu_A", WCoef::CBlnu_A}, 
    {"C_Blnu_P", WCoef::CBlnu_P}, 
    {"C_V1", WCoef::C_V1}, 
    {"C_V2", WCoef::C_V2}, 
    {"C_S1", WCoef::C_S1}, 
    {"C_S2", WCoef::C_S2}, 
    {"C_T", WCoef::C_T}, 
};

const std::vector<WCoef> WCoefMapper::B_group = {
    WCoef::C1, 
    WCoef::C2,
    WCoef::C3,
    WCoef::C4,
    WCoef::C5,
    WCoef::C6,
    WCoef::C7,
    WCoef::C8,
    WCoef::C9,
    WCoef::C10 
};

const std::vector<WCoef> WCoefMapper::B_prime_group = {
    WCoef::CP1, 
    WCoef::CP2,
    WCoef::CP3,
    WCoef::CP4,
    WCoef::CP5,
    WCoef::CP6,
    WCoef::CP7,
    WCoef::CP8,
    WCoef::CP9,
    WCoef::CP10,
    WCoef::CPQ1,
    WCoef::CPQ2, 
};

const std::vector<WCoef> WCoefMapper::B_scalar_group = {
    WCoef::CQ1, 
    WCoef::CQ2 
};

const std::vector<WCoef> WCoefMapper::B_lnu_group = {
    WCoef::CBlnu_A, 
    WCoef::CBlnu_P 
};

const std::vector<WCoef> WCoefMapper::b_clnu_group = {
    WCoef::C_V1, 
    WCoef::C_V2,
    WCoef::C_S1,
    WCoef::C_S2,
    WCoef::C_T
};

const std::map<WGroup, std::string> GroupMapper::mapping = {
    {WGroup::B, "BCoefficients"},
    {WGroup::BPrime, "BPrimeCoefficients"},
    {WGroup::BScalar, "BScalarCoefficients"},
    {WGroup::Blnu, "BlnuCoefficients"},
    {WGroup::BCLNU, "BclnuCoefficients"},
}; 

const std::map<std::string, WGroup> GroupMapper::inverse_mapping = {
    {"BCoefficients", WGroup::B},
    {"BPrimeCoefficients", WGroup::BPrime},
    {"BScalarCoefficients", WGroup::BScalar},
    {"BlnuCoefficients", WGroup::Blnu},
    {"BclnuCoefficients", WGroup::BCLNU},
}; 

const std::map<Observables, std::string> ObservableMapper::mapping = {
    {Observables::BR_BS_MUMU, "531_1_2_13_-13"},
    {Observables::BR_BS_MUMU_UNTAG, "BRuntag_Bs__mu_mu"},
    {Observables::BR_BD_MUMU, "511_1_2_13_-13"},
    {Observables::BR_BU_TAU_NU, "BR_Bu__tau_nu"},
    {Observables::R_TAU_NU, "R_tau_nu"},
    {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, "IA_B__K*_gamma"},
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
};

const std::map<std::string, Observables> ObservableMapper::inverse_mapping = {
    {"531_1_2_13_-13", Observables::BR_BS_MUMU},
    {"BRuntag_Bs__mu_mu", Observables::BR_BS_MUMU_UNTAG},
    {"511_1_2_13_-13", Observables::BR_BD_MUMU},
    {"BR_Bu__tau_nu", Observables::BR_BU_TAU_NU},
    {"R_tau_nu", Observables::R_TAU_NU},
    {"IA_B__K*_gamma", Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA},
    {"BR_B__Xs_gamma", Observables::BR_B_XS_GAMMA},
    {"BR_B__D_tau_nu", Observables::BR_B__D_TAU_NU},
    {"A_FB_B__D_tau_nu", Observables::A_FB_B__D_TAU_NU},
    {"P_tau_B__D_tau_nu", Observables::P_TAU_B__D_TAU_NU},
    {"R_D", Observables::R_D},
    {"BR_B__D*_tau_nu", Observables::BR_B__DSTAR_TAU_NU},
    {"A_FB_B__D*_tau_nu", Observables::A_FB_B__DSTAR_TAU_NU},
    {"P_tau_B__D*_tau_nu", Observables::P_TAU_B__DSTAR_TAU_NU},
    {"P_D_B__D*_tau_nu", Observables::P_D_B__DSTAR_TAU_NU},
    {"R_D*", Observables::R_DSTAR},
};

const std::map<Model, std::string> ModelMapper::mapping = {
    {Model::SM, "SM"},
    {Model::SUSY, "SUSY"},
    {Model::THDM, "THDM"},
    {Model::CUSTOM, "CUSTOM"},
};

const std::map<std::string, Model> ModelMapper::inverse_mapping = {
    {"SM", Model::SM},
    {"SUSY", Model::SUSY},
    {"THDM", Model::THDM},
    {"CUSTOM", Model::CUSTOM},
};


const std::map<ParameterType, std::string> ParameterTypeMapper::mapping = {
    {ParameterType::SM, "SM"},
    {ParameterType::SUSY, "SUSY"},
    {ParameterType::THDM, "THDM"},
    {ParameterType::CUSTOM, "CUSTOM"},
    {ParameterType::FLAVOR, "FLAVOR"},
    {ParameterType::WILSON, "WILSON"},
    {ParameterType::DECAY, "DECAY"},
};

const std::map<std::string, ParameterType> ParameterTypeMapper::inverse_mapping = {
    {"SM", ParameterType::SM},
    {"SUSY", ParameterType::SUSY},
    {"THDM", ParameterType::THDM},
    {"CUSTOM", ParameterType::CUSTOM},
    {"FLAVOR", ParameterType::FLAVOR},
    {"WILSON", ParameterType::WILSON},
    {"DECAY", ParameterType::DECAY},
};

const std::map<Observables, std::vector<ParamId>> DependenciesHelper::dep_lists = {
    {Observables::BR_BS_MUMU, {
        ParamId{ParameterType::SM, "SMINPUTS", 1},
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::SM, "RECKM", 21},
        ParamId{ParameterType::SM, "RECKM", 22},
        ParamId{ParameterType::SM, "IMCKM", 21},
        ParamId{ParameterType::SM, "IMCKM", 22},
        ParamId{ParameterType::FLAVOR, "FMASS", 531},
        ParamId{ParameterType::FLAVOR, "FLIFE", 531},
        ParamId{ParameterType::FLAVOR, "FCONST", 53101}
    }},
    {Observables::BR_BS_MUMU_UNTAG, {
        ParamId{ParameterType::SM, "SMINPUTS", 1},
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::SM, "RECKM", 21},
        ParamId{ParameterType::SM, "RECKM", 22},
        ParamId{ParameterType::SM, "IMCKM", 21},
        ParamId{ParameterType::SM, "IMCKM", 22},
        ParamId{ParameterType::FLAVOR, "FMASS", 531},
        ParamId{ParameterType::FLAVOR, "FLIFE", 531},
        ParamId{ParameterType::FLAVOR, "FCONST", 53101},
        ParamId{ParameterType::DECAY, "B_ll", 1}
    }},
    {Observables::BR_BD_MUMU, {
        ParamId{ParameterType::SM, "SMINPUTS", 1},
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 13},
        ParamId{ParameterType::SM, "RECKM", 20},
        ParamId{ParameterType::SM, "RECKM", 22},
        ParamId{ParameterType::SM, "IMCKM", 20},
        ParamId{ParameterType::SM, "IMCKM", 22},
        ParamId{ParameterType::FLAVOR, "FMASS", 511},
        ParamId{ParameterType::FLAVOR, "FLIFE", 511},
        ParamId{ParameterType::FLAVOR, "FCONST", 51101}
    }},
    {Observables::BR_BU_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::SM, "RECKM", 02},
        // ParamId{ParameterType::SM, "IMCKM", 02},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FLIFE", 521},
        ParamId{ParameterType::FLAVOR, "FCONST", 52101}
    }},
    {Observables::R_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521}
    }},
    {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, {
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "RECKM", 01},
        ParamId{ParameterType::SM, "IMCKM", 01},
        ParamId{ParameterType::SM, "RECKM", 02},
        ParamId{ParameterType::SM, "IMCKM", 02},
        ParamId{ParameterType::SM, "RECKM", 11},
        ParamId{ParameterType::SM, "IMCKM", 11},
        ParamId{ParameterType::SM, "RECKM", 12},
        ParamId{ParameterType::SM, "IMCKM", 12},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 323},
        ParamId{ParameterType::FLAVOR, "FCONST", 52101},
        ParamId{ParameterType::FLAVOR, "FCONST", 32301},
        ParamId{ParameterType::FLAVOR, "FCONST", 32302},
        ParamId{ParameterType::DECAY, "B_Ks", 1},
        ParamId{ParameterType::DECAY, "B_Ks", 2},
        ParamId{ParameterType::DECAY, "B_Ks", 3},
        ParamId{ParameterType::DECAY, "B_Ks", 4},
        ParamId{ParameterType::DECAY, "B_Ks", 5},
        ParamId{ParameterType::DECAY, "B_Ks", 6},
        ParamId{ParameterType::DECAY, "B_Ks", 7},
        ParamId{ParameterType::DECAY, "B_Ks", 8},
        ParamId{ParameterType::DECAY, "B_Ks", 9},
        ParamId{ParameterType::DECAY, "B_Ks", 10},
        ParamId{ParameterType::DECAY, "B_Ks", 11},
        ParamId{ParameterType::DECAY, "B_Ks", 12},
        ParamId{ParameterType::DECAY, "B_Ks", 13}
    }},
    {Observables::BR_B_XS_GAMMA, {
        ParamId{ParameterType::SM, "SMINPUTS", 1},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "RECKM", 12},
        // ParamId{ParameterType::SM, "IMCKM", 12},
        ParamId{ParameterType::SM, "RECKM", 21},
        // ParamId{ParameterType::SM, "IMCKM", 21},
        ParamId{ParameterType::SM, "RECKM", 22},
        // ParamId{ParameterType::SM, "IMCKM", 22},
        ParamId{ParameterType::DECAY, "B_Xs", 1},
        ParamId{ParameterType::DECAY, "B_Xs", 2},
        ParamId{ParameterType::DECAY, "B_Xs", 3},
        ParamId{ParameterType::DECAY, "B_Xs", 4},
        ParamId{ParameterType::DECAY, "B_Xs", 5},
        ParamId{ParameterType::DECAY, "B_Xs", 6},
        ParamId{ParameterType::DECAY, "B_Xs", 7},
        ParamId{ParameterType::DECAY, "B_Xs", 8},
        ParamId{ParameterType::DECAY, "B_Xs", 9},
    }},
    {Observables::BR_B__D_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::SM, "RECKM", 12},
        ParamId{ParameterType::SM, "IMCKM", 12},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 521},
        ParamId{ParameterType::DECAY, "B_Dlnu", 1},
        ParamId{ParameterType::DECAY, "B_Dlnu", 2},
        ParamId{ParameterType::DECAY, "B_Dlnu", 3},
    }},
    {Observables::A_FB_B__D_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 521},
        ParamId{ParameterType::DECAY, "B_Dlnu", 2},
        ParamId{ParameterType::DECAY, "B_Dlnu", 3},
    }},
    {Observables::P_TAU_B__D_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 521},
        ParamId{ParameterType::DECAY, "B_Dlnu", 2},
        ParamId{ParameterType::DECAY, "B_Dlnu", 3},
    }},
    {Observables::R_D, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 521},
        ParamId{ParameterType::DECAY, "B_Dlnu", 2},
        ParamId{ParameterType::DECAY, "B_Dlnu", 3},
    }},
    {Observables::BR_B__DSTAR_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::SM, "RECKM", 12},
        ParamId{ParameterType::SM, "IMCKM", 12},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 523},
        ParamId{ParameterType::DECAY, "B_Dslnu", 1},
        ParamId{ParameterType::DECAY, "B_Dslnu", 2},
        ParamId{ParameterType::DECAY, "B_Dslnu", 3},
        ParamId{ParameterType::DECAY, "B_Dslnu", 4},
    }},
    {Observables::A_FB_B__DSTAR_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 523},
        ParamId{ParameterType::DECAY, "B_Dslnu", 2},
        ParamId{ParameterType::DECAY, "B_Dslnu", 3},
        ParamId{ParameterType::DECAY, "B_Dslnu", 4},
    }},
    {Observables::P_TAU_B__DSTAR_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 523},
        ParamId{ParameterType::DECAY, "B_Dslnu", 2},
        ParamId{ParameterType::DECAY, "B_Dslnu", 3},
        ParamId{ParameterType::DECAY, "B_Dslnu", 4},
    }},
    {Observables::P_D_B__DSTAR_TAU_NU, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 523},
        ParamId{ParameterType::DECAY, "B_Dslnu", 2},
        ParamId{ParameterType::DECAY, "B_Dslnu", 3},
        ParamId{ParameterType::DECAY, "B_Dslnu", 4},
    }},
    {Observables::R_DSTAR, {
        ParamId{ParameterType::SM, "SMINPUTS", 2},
        ParamId{ParameterType::SM, "SMINPUTS", 3},
        ParamId{ParameterType::SM, "SMINPUTS", 4},
        ParamId{ParameterType::SM, "SMINPUTS", 5},
        ParamId{ParameterType::SM, "SMINPUTS", 6},
        ParamId{ParameterType::SM, "MASS", 1},
        ParamId{ParameterType::SM, "MASS", 2},
        ParamId{ParameterType::SM, "MASS", 3},
        ParamId{ParameterType::SM, "MASS", 4},
        ParamId{ParameterType::SM, "MASS", 15},
        ParamId{ParameterType::FLAVOR, "FMASS", 521},
        ParamId{ParameterType::FLAVOR, "FMASS", 421},
        ParamId{ParameterType::FLAVOR, "FLIFE", 523},
        ParamId{ParameterType::DECAY, "B_Dslnu", 2},
        ParamId{ParameterType::DECAY, "B_Dslnu", 3},
        ParamId{ParameterType::DECAY, "B_Dslnu", 4},
    }}
};

std::vector<ParamId> DependenciesHelper::get_allowed_parameters(Observables id) {
    return dep_lists.at(id);
}

bool DependenciesHelper::is_param_allowed(Observables id, ParamId pid) {
    auto allowed = dep_lists.at(id);
    return std::find(allowed.begin(), allowed.end(), pid) != allowed.end();
}

const std::map<std::string, std::vector<std::vector<long>>> LhaParamsHelper::minimal_blocks = {
    {"FMASS", {{211}, {321}, {323}, {411}, {421}, {423}, {431}, {511}, {521}, {531}}},
    {"FLIFE", {{211}, {321}, {323}, {411}, {421}, {431}, {511}, {521}, {531}}},
    {"FCONST", {{511, 1}, {521, 1}, {531, 1}, {323, 1}, {323, 2}}},
};

std::vector<std::vector<long>> LhaParamsHelper::get_minimal_content(const std::string &block_name) {
    if (LhaParamsHelper::minimal_blocks.contains(block_name)) {
        return LhaParamsHelper::minimal_blocks.at(block_name);
    }
    LOG_ERROR("LhaParamsHelper", "Unknown block", block_name);
}

std::ostream &operator<<(std::ostream &os, const LhaID &id) {
    os << id.to_string();
    return os;
}

LhaID::LhaID(const std::string &str_id) {
    for (const auto &num : split(str_id, '_')) {
        parts.emplace_back(std::stol(num));
    }
}

std::string LhaID::to_string() const {
    std::stringstream ss;
    if (!this->parts.empty()) {
        ss << this->parts.at(0);
        for (size_t i = 1; i < this->parts.size(); i++) {
            ss << '_' << this->parts.at(i);
        }
    }
    return ss.str();
}
