#include "General.h"

ObservableMapper* ObservableMapper::instance;

const std::map<CoefficientOrder, std::string> OrderMapper::mapping = {
    {CoefficientOrder::NONE, "None"}, 
    {CoefficientOrder::LO, "LO"}, 
    {CoefficientOrder::NLO, "NLO"}, 
    {CoefficientOrder::NNLO, "NNLO"}
};

const std::map<std::string, CoefficientOrder> OrderMapper::inverse_mapping = {
    {"None", CoefficientOrder::NONE}, 
    {"LO", CoefficientOrder::LO}, 
    {"NLO", CoefficientOrder::NLO}, 
    {"NNLO", CoefficientOrder::NNLO}
};

const std::map<WilsonCoefficientList, std::string> WCoefMapper::mapping = {
    {WilsonCoefficientList::C1, "C1"},
    {WilsonCoefficientList::C2, "C2"},
    {WilsonCoefficientList::C3, "C3"},
    {WilsonCoefficientList::C4, "C4"},
    {WilsonCoefficientList::C5, "C5"},
    {WilsonCoefficientList::C6, "C6"},
    {WilsonCoefficientList::C7, "C7"},
    {WilsonCoefficientList::C8, "C8"},
    {WilsonCoefficientList::C9, "C9"},
    {WilsonCoefficientList::C10, "C10"},
    {WilsonCoefficientList::CP1, "CP1"},
    {WilsonCoefficientList::CP2, "CP2"},
    {WilsonCoefficientList::CP3, "CP3"},
    {WilsonCoefficientList::CP4, "CP4"},
    {WilsonCoefficientList::CP5, "CP5"},
    {WilsonCoefficientList::CP6, "CP6"},
    {WilsonCoefficientList::CP7, "CP7"},
    {WilsonCoefficientList::CP8, "CP8"},
    {WilsonCoefficientList::CP9, "CP9"},
    {WilsonCoefficientList::CP10, "CP10"},
    {WilsonCoefficientList::CQ1, "CQ1"},
    {WilsonCoefficientList::CQ2, "CQ2"},
    {WilsonCoefficientList::CPQ1, "CPQ1"},
    {WilsonCoefficientList::CPQ2, "CPQ2"}, 
};

const std::map<WilsonCoefficientList, std::string> WCoefMapper::flha_mapping = {
    {WilsonCoefficientList::C1, "03040405|6161"},
    {WilsonCoefficientList::C2, "03040405|4141"},
    {WilsonCoefficientList::C3, "03050707|4133"},
    {WilsonCoefficientList::C4, "03050707|6153"},
    {WilsonCoefficientList::C5, "03050707|4536"},
    {WilsonCoefficientList::C6, "03050707|6556"},
    {WilsonCoefficientList::C7, "00000305|4422"},
    {WilsonCoefficientList::C8, "00000305|6421"},
    {WilsonCoefficientList::C9, "03051313|4133"},
    {WilsonCoefficientList::C10, "03051313|4137"},
    {WilsonCoefficientList::CP1, "03040405|6262"},
    {WilsonCoefficientList::CP2, "03040405|4242"},
    {WilsonCoefficientList::CP3, "03050707|4233"},
    {WilsonCoefficientList::CP4, "03050707|6253"},
    {WilsonCoefficientList::CP5, "03050707|4636"},
    {WilsonCoefficientList::CP6, "03050707|6656"},
    {WilsonCoefficientList::CP7, "00000305|4322"},
    {WilsonCoefficientList::CP8, "00000305|4321"},
    {WilsonCoefficientList::CP9, "03051313|4233"},
    {WilsonCoefficientList::CP10, "03051313|4234"},
    {WilsonCoefficientList::CQ1, "03051313|3230"},
    {WilsonCoefficientList::CQ2, "03051313|3233"},
    {WilsonCoefficientList::CPQ1, "03051313|3130"},
    {WilsonCoefficientList::CPQ2, "03051313|3133"}, 
};

const std::map<std::string, WilsonCoefficientList> WCoefMapper::inverse_mapping = {
    {"C1", WilsonCoefficientList::C1},
    {"C2", WilsonCoefficientList::C2},
    {"C3", WilsonCoefficientList::C3},
    {"C4", WilsonCoefficientList::C4},
    {"C5", WilsonCoefficientList::C5},
    {"C6", WilsonCoefficientList::C6},
    {"C7", WilsonCoefficientList::C7},
    {"C8", WilsonCoefficientList::C8},
    {"C9", WilsonCoefficientList::C9},
    {"C10", WilsonCoefficientList::C10},
    {"CP1", WilsonCoefficientList::CP1},
    {"CP2", WilsonCoefficientList::CP2},
    {"CP3", WilsonCoefficientList::CP3},
    {"CP4", WilsonCoefficientList::CP4},
    {"CP5", WilsonCoefficientList::CP5},
    {"CP6", WilsonCoefficientList::CP6},
    {"CP7", WilsonCoefficientList::CP7},
    {"CP8", WilsonCoefficientList::CP8},
    {"CP9", WilsonCoefficientList::CP9},
    {"CP10", WilsonCoefficientList::CP10},
    {"CQ1", WilsonCoefficientList::CQ1},
    {"CQ2", WilsonCoefficientList::CQ2},
    {"CPQ1", WilsonCoefficientList::CPQ1},
    {"CPQ2", WilsonCoefficientList::CPQ2} 
};

const std::vector<WilsonCoefficientList> WCoefMapper::B_group = {
    WilsonCoefficientList::C1, 
    WilsonCoefficientList::C2,
    WilsonCoefficientList::C3,
    WilsonCoefficientList::C4,
    WilsonCoefficientList::C5,
    WilsonCoefficientList::C6,
    WilsonCoefficientList::C7,
    WilsonCoefficientList::C8,
    WilsonCoefficientList::C9,
    WilsonCoefficientList::C10 
};

const std::vector<WilsonCoefficientList> WCoefMapper::B_prime_group = {
    WilsonCoefficientList::CP1, 
    WilsonCoefficientList::CP2,
    WilsonCoefficientList::CP3,
    WilsonCoefficientList::CP4,
    WilsonCoefficientList::CP5,
    WilsonCoefficientList::CP6,
    WilsonCoefficientList::CP7,
    WilsonCoefficientList::CP8,
    WilsonCoefficientList::CP9,
    WilsonCoefficientList::CP10,
    WilsonCoefficientList::CPQ1,
    WilsonCoefficientList::CPQ2, 
};

const std::vector<WilsonCoefficientList> WCoefMapper::B_scalar_group = {
    WilsonCoefficientList::CQ1, 
    WilsonCoefficientList::CQ2 
};

const std::map<WilsonGroups, std::string> GroupMapper::mapping = {
    {WilsonGroups::BCoefficients, "BCoefficients"},
    {WilsonGroups::BPrimeCoefficients, "BPrimeCoefficients"},
    {WilsonGroups::BScalarCoefficients, "BScalarCoefficients"},
}; 

const std::map<std::string, WilsonGroups> GroupMapper::inverse_mapping = {
    {"BCoefficients", WilsonGroups::BCoefficients},
    {"BPrimeCoefficients", WilsonGroups::BPrimeCoefficients},
    {"BScalarCoefficients", WilsonGroups::BScalarCoefficients},
}; 