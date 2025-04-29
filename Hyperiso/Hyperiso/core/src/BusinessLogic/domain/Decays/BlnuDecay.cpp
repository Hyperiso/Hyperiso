#include "BlnuDecay.h"


complex_t BlnuDecay::R(double m_B, double m_b, double m_tau) {
    complex_t C_A = w_proxy->getFM(WGroup::Blnu, WCoef::CBlnu_A, QCDOrder::LO);
    complex_t C_P = w_proxy->getFM(WGroup::Blnu, WCoef::CBlnu_P, QCDOrder::LO);

    return std::pow(std::abs(C_A + std::pow(m_B, 2) * C_P / (m_b * m_tau)), 2);
}

double BlnuDecay::ckm(complex_t V_ub) {
    return std::pow(std::abs(V_ub), 2);
}

double BlnuDecay::pref(double G_F,
                       double f_B,
                       double tau_B,
                       double m_B,
                       double m_tau)
{
    double beta = 1 - std::pow(m_tau / m_B, 2);
    return std::pow(G_F * f_B * m_tau * beta, 2) * tau_B * m_B / (8 * PI * HBAR);
}

double BlnuDecay::BR_B_taunu(double pref, double ckm, double R) {
    return pref * ckm * R;
}

void BlnuDecay::build_op_tree() {

    // SM Parameters
    auto G_F = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 2));
    auto m_tau = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 15));
    auto V_ub = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(0, 2)));
    auto m_b = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "QCD", LhaID(5, 1)));

    // Flavor Parameters
    auto m_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 521));
    auto life_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FLIFE", 521));
    auto f_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FCONST", LhaID(521, 1)));

    auto dummy = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 511));

    // Operator nodes
    auto R_tau_nu = std::make_shared<OperatorNode>("R_tau_nu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->R(values[0], values[1], values[2]); });
    R_tau_nu->addChildren({m_B, m_b, m_tau, dummy});
    roots.emplace(Observables::R_TAU_NU, R_tau_nu);
    auto ckm = std::make_shared<OperatorNode>("ckm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->ckm(values[0]); });
    ckm->addChildren({V_ub});
    auto prefactor = std::make_shared<OperatorNode>("prefactor", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->pref(values[0], values[1], values[2], values[3], values[4]); });
    prefactor->addChildren({G_F, f_B, life_B, m_B, m_tau});
    auto BR_Bu_tau_nu = std::make_shared<OperatorNode>("BR_Bu__tau_nu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_B_taunu(values[0], values[1], values[2]); });
    BR_Bu_tau_nu->addChildren({prefactor, ckm, R_tau_nu});
    roots.emplace(Observables::BR_BU_TAU_NU, BR_Bu_tau_nu);
}
