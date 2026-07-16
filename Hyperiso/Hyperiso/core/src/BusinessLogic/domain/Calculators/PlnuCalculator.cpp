#include "PlnuCalculator.h"

PlnuCalculator::PlnuCalculator(int P_id, int l_id, complex_t C_A, complex_t C_P, 
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p) : iobspp_sm(p) {
    if (!allowed_P.contains(P_id)) {
        LOG_ERROR("Logic Error", P_id, "is not a pseudoscalar meson");
    }

    this->m_P = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", P_id}, DataType::VALUE);
    this->tau_P = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", P_id}, DataType::VALUE) / HBAR;
    this->f_P = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {P_id, 1}}, DataType::VALUE);
    this->m_l = (*p)(ParamId{ParameterType::SM, "MASS", l_id}, DataType::VALUE);
    this->m_qu = (*p)(ParamId{ParameterType::SM, "MASS", allowed_P.at(P_id)[0]}, DataType::VALUE);
    this->m_qd = allowed_P.at(P_id)[1] == 5 ? (*p)(ParamId{ParameterType::SM, "SMINPUTS", 5}, DataType::VALUE) : (*p)(ParamId{ParameterType::SM, "MASS", allowed_P.at(P_id)[1]}, DataType::VALUE);
    this->V_sq = std::pow(std::abs((*p)(ParamId{ParameterType::SM, "VCKM", {allowed_P.at(P_id)[2], allowed_P.at(P_id)[3]}}, DataType::VALUE)), 2);
    this->C_A = C_A;
    this->C_P = C_P;
}

double PlnuCalculator::BR_0_SM() {
    double G_F = (*iobspp_sm)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE);
    double r_l = std::pow(m_l / m_P, 2);

    return std::pow(G_F * f_P * m_l * (1 - r_l), 2) * V_sq * m_P * tau_P / (8 * PI);
}

double PlnuCalculator::R_SM_BSM() {
    return std::pow(std::abs(
        C_A - std::pow(m_P, 2) / (m_l * (m_qu + m_qd)) * C_P
    ), 2);
}
