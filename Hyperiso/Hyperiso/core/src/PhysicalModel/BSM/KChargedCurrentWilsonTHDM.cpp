#include "KChargedCurrentWilsonTHDM.h"

C_S1_su_THDM::C_S1_su_THDM()
    : WilsonCoefficient("C_S1_su_THDM", GroupMapper::str(WGroup::BCC_su, ScaleType::MATCHING))
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
        LhaID(4051516, 3231, 0, 1)
    };
}

double C_S1_su_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double ld     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 8})->get_val();
    double mH     = src.at({ParameterType::BSM, "MASS", 37})->get_val();
    double m_tau  = src.at({ParameterType::SM, "MASS", 15})->get_val();
    double l_tau  = src.at({ParameterType::BSM, "YE", {3, 3}})->get_val();
    double m_b    = src.at({ParameterType::SM, "QCD", {5, 1}})->get_val();

    return -m_b * m_tau * ld * l_tau / std::pow(mH, 2);
}

C_S2_su_THDM::C_S2_su_THDM()
    : WilsonCoefficient("C_S2_su_THDM", GroupMapper::str(WGroup::BCC_su, ScaleType::MATCHING))
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
        LhaID(4051516, 3131, 0, 1)
    };
}

double C_S2_su_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {
    double lu     = src.at({ParameterType::WILSON, "WPARAM_SI_BSM", 7})->get_val();
    double mH     = src.at({ParameterType::BSM, "MASS", 37})->get_val();
    double m_c    = src.at({ParameterType::SM, "MASS", 4})->get_val();    // m_c
    double m_tau  = src.at({ParameterType::SM, "MASS", 15})->get_val();   // m_tau
    double l_tau  = src.at({ParameterType::BSM, "YE", {3, 3}})->get_val(); // l_tau

    return -m_c * m_tau * lu * l_tau / std::pow(mH, 2);
}
