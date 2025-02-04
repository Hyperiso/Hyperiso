#include "General.h"

const std::map<Decays, std::vector<Observables>> DecayMapper::obs_mapping = {
    {Decays::B__D_l_nu, {Observables::BR_B__D_TAU_NU, Observables::XI__D_L_NU}},
    {Decays::B__Kstar,  {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA}},
    {Decays::B__l_l,    {Observables::BR_BS_MUMU, Observables::BR_BD_MUMU, Observables::BR_BS_MUMU_UNTAG}},
    {Decays::B__l_nu,   {Observables::BR_BU_TAU_NU, Observables::R_TAU_NU}},
    {Decays::B__Xs,     {Observables::BR_B_XS_GAMMA}},
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
};

const std::map<WCoef, std::string> WCoefMapper::flha_mapping = {
    {WCoef::C1, "03040405|6161"},
    {WCoef::C2, "03040405|4141"},
    {WCoef::C3, "03050707|4133"},
    {WCoef::C4, "03050707|6153"},
    {WCoef::C5, "03050707|4536"},
    {WCoef::C6, "03050707|6556"},
    {WCoef::C7, "00000305|4422"},
    {WCoef::C8, "00000305|6421"},
    {WCoef::C9, "03051313|4133"},
    {WCoef::C10, "03051313|4137"},
    {WCoef::CP1, "03040405|6262"},
    {WCoef::CP2, "03040405|4242"},
    {WCoef::CP3, "03050707|4233"},
    {WCoef::CP4, "03050707|6253"},
    {WCoef::CP5, "03050707|4636"},
    {WCoef::CP6, "03050707|6656"},
    {WCoef::CP7, "00000305|4322"},
    {WCoef::CP8, "00000305|4321"},
    {WCoef::CP9, "03051313|4233"},
    {WCoef::CP10, "03051313|4234"},
    {WCoef::CQ1, "03051313|3230"},
    {WCoef::CQ2, "03051313|3233"},
    {WCoef::CPQ1, "03051313|3130"},
    {WCoef::CPQ2, "03051313|3133"}, 
};

const std::map<std::string, WCoef> WCoefMapper::inverse_flha_mapping = {
    {"03040405|6161", WCoef::C1},
    {"03040405|4141", WCoef::C2},
    {"03050707|4133", WCoef::C3},
    {"03050707|6153", WCoef::C4},
    {"03050707|4536", WCoef::C5},
    {"03050707|6556", WCoef::C6},
    {"00000305|4422", WCoef::C7},
    {"00000305|6421", WCoef::C8},
    {"03051313|4133", WCoef::C9},
    {"03051313|4137", WCoef::C10},
    {"03040405|6262", WCoef::CP1},
    {"03040405|4242", WCoef::CP2},
    {"03050707|4233", WCoef::CP3},
    {"03050707|6253", WCoef::CP4},
    {"03050707|4636", WCoef::CP5},
    {"03050707|6656", WCoef::CP6},
    {"00000305|4322", WCoef::CP7},
    {"00000305|4321", WCoef::CP8},
    {"03051313|4233", WCoef::CP9},
    {"03051313|4234", WCoef::CP10},
    {"03051313|3230", WCoef::CQ1},
    {"03051313|3233", WCoef::CQ2},
    {"03051313|3130", WCoef::CPQ1},
    {"03051313|3133", WCoef::CPQ2}, 
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

const std::map<WGroup, std::string> GroupMapper::mapping = {
    {WGroup::B, "BCoefficients"},
    {WGroup::BPrime, "BPrimeCoefficients"},
    {WGroup::BScalar, "BScalarCoefficients"},
    {WGroup::Blnu, "BlnuCoefficients"},
}; 

const std::map<std::string, WGroup> GroupMapper::inverse_mapping = {
    {"BCoefficients", WGroup::B},
    {"BPrimeCoefficients", WGroup::BPrime},
    {"BScalarCoefficients", WGroup::BScalar},
    {"BlnuCoefficients", WGroup::Blnu},
}; 

const std::map<Observables, std::string> ObservableMapper::mapping = {
    {Observables::BR_BS_MUMU, "BR_Bs__mu_mu"},
    {Observables::BR_BS_MUMU_UNTAG, "BRuntag_Bs__mu_mu"},
    {Observables::BR_BD_MUMU, "BR_Bd__mu_mu"},
    {Observables::BR_BU_TAU_NU, "BR_Bu__tau_nu"},
    {Observables::R_TAU_NU, "R_tau_nu"},
    {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, "IA_B_K*__gamma"},
    {Observables::BR_B_XS_GAMMA, "BR_B__Xs_gamma"},
};

const std::map<std::string, Observables> ObservableMapper::inverse_mapping = {
    {"BR_Bs__mu_mu", Observables::BR_BS_MUMU},
    {"BRuntag_Bs__mu_mu", Observables::BR_BS_MUMU_UNTAG},
    {"BR_Bd__mu_mu", Observables::BR_BD_MUMU},
    {"BR_Bu__tau_nu", Observables::BR_BU_TAU_NU},
    {"R_tau_nu", Observables::R_TAU_NU},
    {"IA_B_K*__gamma", Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA},
    {"BR_B__Xs_gamma", Observables::BR_B_XS_GAMMA}
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
    {ParameterType::FF, "FF"},
};

const std::map<std::string, ParameterType> ParameterTypeMapper::inverse_mapping = {
    {"SM", ParameterType::SM},
    {"SUSY", ParameterType::SUSY},
    {"THDM", ParameterType::THDM},
    {"CUSTOM", ParameterType::CUSTOM},
    {"FLAVOR", ParameterType::FLAVOR},
    {"WILSON", ParameterType::WILSON},
    {"FF", ParameterType::FF},
};