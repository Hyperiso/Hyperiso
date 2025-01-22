#include "General.h"

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

const std::map<BWilsonCoefficients, std::string> WCoefMapper::mapping = {
    {BWilsonCoefficients::C1, "C1"},
    {BWilsonCoefficients::C2, "C2"},
    {BWilsonCoefficients::C3, "C3"},
    {BWilsonCoefficients::C4, "C4"},
    {BWilsonCoefficients::C5, "C5"},
    {BWilsonCoefficients::C6, "C6"},
    {BWilsonCoefficients::C7, "C7"},
    {BWilsonCoefficients::C8, "C8"},
    {BWilsonCoefficients::C9, "C9"},
    {BWilsonCoefficients::C10, "C10"},
    {BWilsonCoefficients::CP1, "CP1"},
    {BWilsonCoefficients::CP2, "CP2"},
    {BWilsonCoefficients::CP3, "CP3"},
    {BWilsonCoefficients::CP4, "CP4"},
    {BWilsonCoefficients::CP5, "CP5"},
    {BWilsonCoefficients::CP6, "CP6"},
    {BWilsonCoefficients::CP7, "CP7"},
    {BWilsonCoefficients::CP8, "CP8"},
    {BWilsonCoefficients::CP9, "CP9"},
    {BWilsonCoefficients::CP10, "CP10"},
    {BWilsonCoefficients::CQ1, "CQ1"},
    {BWilsonCoefficients::CQ2, "CQ2"},
    {BWilsonCoefficients::CPQ1, "CPQ1"},
    {BWilsonCoefficients::CPQ2, "CPQ2"}, 
};

const std::map<BWilsonCoefficients, std::string> WCoefMapper::flha_mapping = {
    {BWilsonCoefficients::C1, "03040405|6161"},
    {BWilsonCoefficients::C2, "03040405|4141"},
    {BWilsonCoefficients::C3, "03050707|4133"},
    {BWilsonCoefficients::C4, "03050707|6153"},
    {BWilsonCoefficients::C5, "03050707|4536"},
    {BWilsonCoefficients::C6, "03050707|6556"},
    {BWilsonCoefficients::C7, "00000305|4422"},
    {BWilsonCoefficients::C8, "00000305|6421"},
    {BWilsonCoefficients::C9, "03051313|4133"},
    {BWilsonCoefficients::C10, "03051313|4137"},
    {BWilsonCoefficients::CP1, "03040405|6262"},
    {BWilsonCoefficients::CP2, "03040405|4242"},
    {BWilsonCoefficients::CP3, "03050707|4233"},
    {BWilsonCoefficients::CP4, "03050707|6253"},
    {BWilsonCoefficients::CP5, "03050707|4636"},
    {BWilsonCoefficients::CP6, "03050707|6656"},
    {BWilsonCoefficients::CP7, "00000305|4322"},
    {BWilsonCoefficients::CP8, "00000305|4321"},
    {BWilsonCoefficients::CP9, "03051313|4233"},
    {BWilsonCoefficients::CP10, "03051313|4234"},
    {BWilsonCoefficients::CQ1, "03051313|3230"},
    {BWilsonCoefficients::CQ2, "03051313|3233"},
    {BWilsonCoefficients::CPQ1, "03051313|3130"},
    {BWilsonCoefficients::CPQ2, "03051313|3133"}, 
};

const std::map<std::string, BWilsonCoefficients> WCoefMapper::inverse_flha_mapping = {
    {"03040405|6161", BWilsonCoefficients::C1},
    {"03040405|4141", BWilsonCoefficients::C2},
    {"03050707|4133", BWilsonCoefficients::C3},
    {"03050707|6153", BWilsonCoefficients::C4},
    {"03050707|4536", BWilsonCoefficients::C5},
    {"03050707|6556", BWilsonCoefficients::C6},
    {"00000305|4422", BWilsonCoefficients::C7},
    {"00000305|6421", BWilsonCoefficients::C8},
    {"03051313|4133", BWilsonCoefficients::C9},
    {"03051313|4137", BWilsonCoefficients::C10},
    {"03040405|6262", BWilsonCoefficients::CP1},
    {"03040405|4242", BWilsonCoefficients::CP2},
    {"03050707|4233", BWilsonCoefficients::CP3},
    {"03050707|6253", BWilsonCoefficients::CP4},
    {"03050707|4636", BWilsonCoefficients::CP5},
    {"03050707|6656", BWilsonCoefficients::CP6},
    {"00000305|4322", BWilsonCoefficients::CP7},
    {"00000305|4321", BWilsonCoefficients::CP8},
    {"03051313|4233", BWilsonCoefficients::CP9},
    {"03051313|4234", BWilsonCoefficients::CP10},
    {"03051313|3230", BWilsonCoefficients::CQ1},
    {"03051313|3233", BWilsonCoefficients::CQ2},
    {"03051313|3130", BWilsonCoefficients::CPQ1},
    {"03051313|3133", BWilsonCoefficients::CPQ2}, 
};

const std::map<std::string, BWilsonCoefficients> WCoefMapper::inverse_mapping = {
    {"C1", BWilsonCoefficients::C1},
    {"C2", BWilsonCoefficients::C2},
    {"C3", BWilsonCoefficients::C3},
    {"C4", BWilsonCoefficients::C4},
    {"C5", BWilsonCoefficients::C5},
    {"C6", BWilsonCoefficients::C6},
    {"C7", BWilsonCoefficients::C7},
    {"C8", BWilsonCoefficients::C8},
    {"C9", BWilsonCoefficients::C9},
    {"C10", BWilsonCoefficients::C10},
    {"CP1", BWilsonCoefficients::CP1},
    {"CP2", BWilsonCoefficients::CP2},
    {"CP3", BWilsonCoefficients::CP3},
    {"CP4", BWilsonCoefficients::CP4},
    {"CP5", BWilsonCoefficients::CP5},
    {"CP6", BWilsonCoefficients::CP6},
    {"CP7", BWilsonCoefficients::CP7},
    {"CP8", BWilsonCoefficients::CP8},
    {"CP9", BWilsonCoefficients::CP9},
    {"CP10", BWilsonCoefficients::CP10},
    {"CQ1", BWilsonCoefficients::CQ1},
    {"CQ2", BWilsonCoefficients::CQ2},
    {"CPQ1", BWilsonCoefficients::CPQ1},
    {"CPQ2", BWilsonCoefficients::CPQ2} 
};

const std::vector<BWilsonCoefficients> WCoefMapper::B_group = {
    BWilsonCoefficients::C1, 
    BWilsonCoefficients::C2,
    BWilsonCoefficients::C3,
    BWilsonCoefficients::C4,
    BWilsonCoefficients::C5,
    BWilsonCoefficients::C6,
    BWilsonCoefficients::C7,
    BWilsonCoefficients::C8,
    BWilsonCoefficients::C9,
    BWilsonCoefficients::C10 
};

const std::vector<BWilsonCoefficients> WCoefMapper::B_prime_group = {
    BWilsonCoefficients::CP1, 
    BWilsonCoefficients::CP2,
    BWilsonCoefficients::CP3,
    BWilsonCoefficients::CP4,
    BWilsonCoefficients::CP5,
    BWilsonCoefficients::CP6,
    BWilsonCoefficients::CP7,
    BWilsonCoefficients::CP8,
    BWilsonCoefficients::CP9,
    BWilsonCoefficients::CP10,
    BWilsonCoefficients::CPQ1,
    BWilsonCoefficients::CPQ2, 
};

const std::vector<BWilsonCoefficients> WCoefMapper::B_scalar_group = {
    BWilsonCoefficients::CQ1, 
    BWilsonCoefficients::CQ2 
};

const std::map<WilsonGroups, std::string> GroupMapper::mapping = {
    {WilsonGroups::BCoefficients, "BCoefficients"},
    {WilsonGroups::BPrimeCoefficients, "BPrimeCoefficients"},
    {WilsonGroups::BScalarCoefficients, "BScalarCoefficients"},
    {WilsonGroups::BCoefficients_THDM, "BCoefficients_THDM"},
    {WilsonGroups::BPrimeCoefficients_THDM, "BPrimeCoefficients_THDM"},
    {WilsonGroups::BScalarCoefficients_THDM, "BScalarCoefficients_THDM"},
    {WilsonGroups::BCoefficients_SUSY, "BCoefficients_SUSY"},
    {WilsonGroups::BPrimeCoefficients_SUSY, "BPrimeCoefficients_SUSY"},
    {WilsonGroups::BScalarCoefficients_SUSY, "BScalarCoefficients_SUSY"},
    {WilsonGroups::BCoefficients_MARTY, "BCoefficients_MARTY"},
    {WilsonGroups::BPrimeCoefficients_MARTY, "BPrimeCoefficients_MARTY"},
    {WilsonGroups::BScalarCoefficients_MARTY, "BScalarCoefficients_MARTY"},
}; 

const std::map<std::string, WilsonGroups> GroupMapper::inverse_mapping = {
    {"BCoefficients", WilsonGroups::BCoefficients},
    {"BPrimeCoefficients", WilsonGroups::BPrimeCoefficients},
    {"BScalarCoefficients", WilsonGroups::BScalarCoefficients},
    {"BCoefficients_THDM", WilsonGroups::BCoefficients_THDM},
    {"BPrimeCoefficients_THDM", WilsonGroups::BPrimeCoefficients_THDM},
    {"BScalarCoefficients_THDM", WilsonGroups::BScalarCoefficients_THDM},
    {"BCoefficients_SUSY", WilsonGroups::BCoefficients_SUSY},
    {"BPrimeCoefficients_SUSY", WilsonGroups::BPrimeCoefficients_SUSY},
    {"BScalarCoefficients_SUSY", WilsonGroups::BScalarCoefficients_SUSY},
    {"BCoefficients_MARTY", WilsonGroups::BCoefficients_MARTY},
    {"BPrimeCoefficients_MARTY", WilsonGroups::BPrimeCoefficients_MARTY},
    {"BScalarCoefficients_MARTY", WilsonGroups::BScalarCoefficients_MARTY},
}; 

const std::map<Observables, std::string> ObservableMapper::mapping = {
    {Observables::BR_BS_MUMU, "BR_Bs__mu_mu"},
    {Observables::BR_BS_MUMU_UNTAG, "BRuntag_Bs__mu_mu"},
    {Observables::BR_BD_MUMU, "BR_Bd__mu_mu"},
    // {Observables::BR_BU_TAUNU, "BR_Bu__tau_nu"},
    {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, "IA_B_K*__gamma"}
};

const std::map<std::string, Observables> ObservableMapper::inverse_mapping = {
    {"BR_Bs__mu_mu", Observables::BR_BS_MUMU},
    {"BRuntag_Bs__mu_mu", Observables::BR_BS_MUMU_UNTAG},
    {"BR_Bd__mu_mu", Observables::BR_BD_MUMU},
    // {"BR_Bu__tau_nu", Observables::BR_BU_TAUNU},
    {"IA_B_K*__gamma", Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA}
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