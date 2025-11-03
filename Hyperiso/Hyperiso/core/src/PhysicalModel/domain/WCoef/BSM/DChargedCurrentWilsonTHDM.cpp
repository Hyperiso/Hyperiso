#include "DChargedCurrentWilsonTHDM.h"

C_S1_cs_THDM::C_S1_cs_THDM()
    : WilsonCoefficient("C_S1_cs_THDM", GroupMapper::str(WGroup::CC_cs, ScaleType::MATCHING))
{
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_BSM", 8},              // ld
            {ParameterType::SM, "MASS", 15},   // m_tau
            {ParameterType::BSM, "MASS", 37},  // m_H
            {ParameterType::BSM, "YE", {3, 3}},// l_tau
            {ParameterType::SM, "QCD", {5, 1}} // m_b
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_S1_cs_THDM::compute_LO(const ParamSrc& src) {
    double ld     = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double mH     = src.get_val(ParameterType::BSM, "MASS", 37);
    double m_tau  = src.get_val(ParameterType::SM, "MASS", 15);
    double l_tau  = src.get_val(ParameterType::BSM, "YE", {3, 3});
    double m_b    = src.get_val(ParameterType::SM, "QCD", {5, 1});

    return -m_b * m_tau * ld * l_tau / std::pow(mH, 2);
}

C_S2_cs_THDM::C_S2_cs_THDM()
    : WilsonCoefficient("C_S2_cs_THDM", GroupMapper::str(WGroup::CC_cs, ScaleType::MATCHING))
{
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_BSM", 7},              // lu
            {ParameterType::SM, "MASS", 15},   // m_tau
            {ParameterType::SM, "MASS", 4},    // m_c (charm)
            {ParameterType::BSM, "MASS", 37},  // m_H
            {ParameterType::BSM, "YE", {3, 3}} // l_tau
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_S2_cs_THDM::compute_LO(const ParamSrc& src) {
    double lu     = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double mH     = src.get_val(ParameterType::BSM, "MASS", 37);
    double m_c    = src.get_val(ParameterType::SM, "MASS", 4);    // m_c
    double m_tau  = src.get_val(ParameterType::SM, "MASS", 15);   // m_tau
    double l_tau  = src.get_val(ParameterType::BSM, "YE", {3, 3}); // l_tau

    return -m_c * m_tau * lu * l_tau / std::pow(mH, 2);
}



C_S1_cd_THDM::C_S1_cd_THDM()
    : WilsonCoefficient("C_S1_cd_THDM", GroupMapper::str(WGroup::CC_cd, ScaleType::MATCHING))
{
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_BSM", 8},              // ld
            {ParameterType::SM, "MASS", 15},   // m_tau
            {ParameterType::BSM, "MASS", 37},  // m_H
            {ParameterType::BSM, "YE", {3, 3}},// l_tau
            {ParameterType::SM, "QCD", {5, 1}} // m_b
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_S1_cd_THDM::compute_LO(const ParamSrc& src) {
    double ld     = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 8);
    double mH     = src.get_val(ParameterType::BSM, "MASS", 37);
    double m_tau  = src.get_val(ParameterType::SM, "MASS", 15);
    double l_tau  = src.get_val(ParameterType::BSM, "YE", {3, 3});
    double m_b    = src.get_val(ParameterType::SM, "QCD", {5, 1});

    return -m_b * m_tau * ld * l_tau / std::pow(mH, 2);
}

C_S2_cd_THDM::C_S2_cd_THDM()
    : WilsonCoefficient("C_S2_cd_THDM", GroupMapper::str(WGroup::CC_cd, ScaleType::MATCHING))
{
    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_SI_BSM", 7},              // lu
            {ParameterType::SM, "MASS", 15},   // m_tau
            {ParameterType::SM, "MASS", 4},    // m_c (charm)
            {ParameterType::BSM, "MASS", 37},  // m_H
            {ParameterType::BSM, "YE", {3, 3}} // l_tau
        },
        compute_LO,
        get_lhaid_from_name(QCDOrder::LO)
    };
}

double C_S2_cd_THDM::compute_LO(const ParamSrc& src) {
    double lu     = src.get_val(ParameterType::WILSON, "WPARAM_SI_BSM", 7);
    double mH     = src.get_val(ParameterType::BSM, "MASS", 37);
    double m_c    = src.get_val(ParameterType::SM, "MASS", 4);    // m_c
    double m_tau  = src.get_val(ParameterType::SM, "MASS", 15);   // m_tau
    double l_tau  = src.get_val(ParameterType::BSM, "YE", {3, 3}); // l_tau

    return -m_c * m_tau * lu * l_tau / std::pow(mH, 2);
}