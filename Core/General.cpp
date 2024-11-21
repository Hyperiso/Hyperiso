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

const std::map<WilsonCoefficient, std::string> WCoefMapper::mapping = {
    {WilsonCoefficient::C1, "C1"},
    {WilsonCoefficient::C2, "C2"},
    {WilsonCoefficient::C3, "C3"},
    {WilsonCoefficient::C4, "C4"},
    {WilsonCoefficient::C5, "C5"},
    {WilsonCoefficient::C6, "C6"},
    {WilsonCoefficient::C7, "C7"},
    {WilsonCoefficient::C8, "C8"},
    {WilsonCoefficient::C9, "C9"},
    {WilsonCoefficient::C10, "C10"},
    {WilsonCoefficient::CP1, "CP1"},
    {WilsonCoefficient::CP2, "CP2"},
    {WilsonCoefficient::CP3, "CP3"},
    {WilsonCoefficient::CP4, "CP4"},
    {WilsonCoefficient::CP5, "CP5"},
    {WilsonCoefficient::CP6, "CP6"},
    {WilsonCoefficient::CP7, "CP7"},
    {WilsonCoefficient::CP8, "CP8"},
    {WilsonCoefficient::CP9, "CP9"},
    {WilsonCoefficient::CP10, "CP10"},
    {WilsonCoefficient::CQ1, "CQ1"},
    {WilsonCoefficient::CQ2, "CQ2"},
    {WilsonCoefficient::CPQ1, "CPQ1"},
    {WilsonCoefficient::CPQ2, "CPQ2"}, 
};

const std::map<WilsonCoefficient, std::string> WCoefMapper::mapping = {
    {WilsonCoefficient::C1, "03040405|6161"},
    {WilsonCoefficient::C2, "03040405|4141"},
    {WilsonCoefficient::C3, "03050707|4133"},
    {WilsonCoefficient::C4, "03050707|6153"},
    {WilsonCoefficient::C5, "03050707|4536"},
    {WilsonCoefficient::C6, "03050707|6556"},
    {WilsonCoefficient::C7, "00000305|4422"},
    {WilsonCoefficient::C8, "00000305|6421"},
    {WilsonCoefficient::C9, "03051313|4133"},
    {WilsonCoefficient::C10, "03051313|4137"},
    {WilsonCoefficient::CP1, "03040405|6262"},
    {WilsonCoefficient::CP2, "03040405|4242"},
    {WilsonCoefficient::CP3, "03050707|4233"},
    {WilsonCoefficient::CP4, "03050707|6253"},
    {WilsonCoefficient::CP5, "03050707|4636"},
    {WilsonCoefficient::CP6, "03050707|6656"},
    {WilsonCoefficient::CP7, "00000305|4322"},
    {WilsonCoefficient::CP8, "00000305|4321"},
    {WilsonCoefficient::CP9, "03051313|4233"},
    {WilsonCoefficient::CP10, "03051313|4234"},
    {WilsonCoefficient::CQ1, "03051313|3230"},
    {WilsonCoefficient::CQ2, "03051313|3233"},
    {WilsonCoefficient::CPQ1, "03051313|3130"},
    {WilsonCoefficient::CPQ2, "03051313|3133"}, 
};

const std::map<std::string, WilsonCoefficient> WCoefMapper::inverse_mapping = {
    {"C1", WilsonCoefficient::C1},
    {"C2", WilsonCoefficient::C2},
    {"C3", WilsonCoefficient::C3},
    {"C4", WilsonCoefficient::C4},
    {"C5", WilsonCoefficient::C5},
    {"C6", WilsonCoefficient::C6},
    {"C7", WilsonCoefficient::C7},
    {"C8", WilsonCoefficient::C8},
    {"C9", WilsonCoefficient::C9},
    {"C10", WilsonCoefficient::C10},
    {"CP1", WilsonCoefficient::CP1},
    {"CP2", WilsonCoefficient::CP2},
    {"CP3", WilsonCoefficient::CP3},
    {"CP4", WilsonCoefficient::CP4},
    {"CP5", WilsonCoefficient::CP5},
    {"CP6", WilsonCoefficient::CP6},
    {"CP7", WilsonCoefficient::CP7},
    {"CP8", WilsonCoefficient::CP8},
    {"CP9", WilsonCoefficient::CP9},
    {"CP10", WilsonCoefficient::CP10},
    {"CQ1", WilsonCoefficient::CQ1},
    {"CQ2", WilsonCoefficient::CQ2},
    {"CPQ1", WilsonCoefficient::CPQ1},
    {"CPQ2", WilsonCoefficient::CPQ2} 
};

const std::vector<WilsonCoefficient> WCoefMapper::B_group = {
    WilsonCoefficient::C1, 
    WilsonCoefficient::C2,
    WilsonCoefficient::C3,
    WilsonCoefficient::C4,
    WilsonCoefficient::C5,
    WilsonCoefficient::C6,
    WilsonCoefficient::C7,
    WilsonCoefficient::C8,
    WilsonCoefficient::C9,
    WilsonCoefficient::C10 
};

const std::vector<WilsonCoefficient> WCoefMapper::B_prime_group = {
    WilsonCoefficient::CP1, 
    WilsonCoefficient::CP2,
    WilsonCoefficient::CP3,
    WilsonCoefficient::CP4,
    WilsonCoefficient::CP5,
    WilsonCoefficient::CP6,
    WilsonCoefficient::CP7,
    WilsonCoefficient::CP8,
    WilsonCoefficient::CP9,
    WilsonCoefficient::CP10,
    WilsonCoefficient::CPQ1,
    WilsonCoefficient::CPQ2, 
};

const std::vector<WilsonCoefficient> WCoefMapper::B_scalar_group = {
    WilsonCoefficient::CQ1, 
    WilsonCoefficient::CQ2 
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