#include "BDlnuDecay.h"

double BDlnuDecay::ckm(double V_cb_r, double V_cb_i) {
    return std::pow(std::abs(complex_t(V_cb_r, V_cb_i)), 2);
}

double BDlnuDecay::pref(double G_F, double tau_B, double m_B, double m_D, double G1) {
    return std::pow(G_F * m_B * m_D * (1 + m_D / m_B) * G1 / (4 * PI), 2) * m_D * tau_B / (3 * PI * HBAR);
}

double BDlnuDecay::t(double rD, double w) {
    return 1 + rD * (rD - 2 * w);
}

double BDlnuDecay::G(double rho2, double w) {
    double z = (std::sqrt(1 + w) - RT2) / (std::sqrt(w + 1) + RT2);
    return 1 + z * (-8 * rho2 + z * (51 * rho2 - 10 - z * (252 * rho2 - 84)));
}

double BDlnuDecay::f_1(double rD, double rl, double rho2, double w) {
    double a = rl / t(rD, w); 
    return std::pow(w * w - 1, 1.5) * std::pow(1 - a, 2) * (1 + a / 2) * G(rho2, w);
}

double BDlnuDecay::f_2(double rD, double rl, double rho2, double w) {
    double tw = t(rD, w); 
    return std::pow(w * w - 1, 1.5) * (1 + w) * std::pow(1 - rl / tw, 2) * G(rho2, w) / (tw * (1 - w));
}

double BDlnuDecay::f_3(double rD, double rl, double rho2, double w) {
    double tw = t(rD, w); 
    return std::pow(w * w - 1, 1.5) * (1 + w) * std::pow(1 - rl / tw, 2) * G(rho2, w) / (1 - w);
}

double BDlnuDecay::f_4(double rD, double rl, double rho2, double w) {
    double tw = t(rD, w); 
    return std::pow(w * w - 1, 1.5) * (1 + w) * std::pow(1 - rl / tw, 2) * tw * G(rho2, w) / (1 - w);
}

double BDlnuDecay::w_max(double rD, double rl) {
    return (1 + rD * rD - rl) / (2 * rD);
}

double BDlnuDecay::I1(double rD, double rl, double rho2, double w_u) {
    auto f = [this, rD, rl, rho2] (double w) {
        return f_1(rD, rl, rho2, w);
    };

    return integrate(f, 1, w_u, 1e-3);
}

double BDlnuDecay::I2(double rD, double rl, double rho2, double w_u) {
    auto f = [this, rD, rl, rho2] (double w) {
        return f_2(rD, rl, rho2, w);
    };

    return integrate(f, 1, w_u, 1e-3);
}

double BDlnuDecay::I3(double rD, double rl, double rho2, double w_u) {
    auto f = [this, rD, rl, rho2] (double w) {
        return f_3(rD, rl, rho2, w);
    };

    return integrate(f, 1, w_u, 1e-3);
}

double BDlnuDecay::I4(double rD, double rl, double rho2, double w_u) {
    auto f = [this, rD, rl, rho2] (double w) {
        return f_4(rD, rl, rho2, w);
    };

    return integrate(f, 1, w_u, 1e-3);
}

double BDlnuDecay::W1(double I1) {
    complex_t C_A = get_wilsons()->getFM(WGroup::Blnu, WCoef::CBlnu_A, QCDOrder::LO);
    return std::pow(std::abs(C_A), 2) * I1;
}

double BDlnuDecay::W2(double I2, double Delta2, double rl) {
    complex_t C_A = get_wilsons()->getFM(WGroup::Blnu, WCoef::CBlnu_A, QCDOrder::LO);
    return -3 * Delta2 * rl * std::pow(std::abs(C_A), 2) * I2 / 2;
}

double BDlnuDecay::W3(double I3, double Delta2, double m_l, double m_bc) {
    complex_t C_A = get_wilsons()->getFM(WGroup::Blnu, WCoef::CBlnu_A, QCDOrder::LO);
    complex_t C_P = get_wilsons()->getFM(WGroup::Blnu, WCoef::CBlnu_P, QCDOrder::LO);
    return -3 * m_l * Delta2 * std::real(C_A * std::conj(C_P)) * I3 / m_bc;
}

double BDlnuDecay::W4(double I4, double Delta2, double m_B, double m_bc) {
    complex_t C_P = get_wilsons()->getFM(WGroup::Blnu, WCoef::CBlnu_P, QCDOrder::LO);
    return -3 * Delta2 * std::pow(std::abs(C_P) * m_B / m_bc, 2) * I4 / 2;
}

double BDlnuDecay::W(double W1, double W2, double W3, double W4) {
    return W1 + W2 + W3 + W4;
}

double BDlnuDecay::BR_B_Dlnu(double pref, double ckm, double W) {
    return pref * ckm * W;
}

void BDlnuDecay::build_op_tree() {
    // SM Parameters
    auto G_F = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 2));
    auto alpha_s_MZ = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 3));
    auto M_Z = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 4));
    auto mt_pole = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 6));
    auto mb_mb = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 5));

    auto m_tau = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 15));
    auto m_e = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 11));
    auto m_d = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 1));
    auto m_u = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 2));
    auto m_s = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 3));
    auto m_c = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 4));
    auto V_cb_r = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "RECKM", 12));
    auto V_cb_i = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "IMCKM", 12));

    // Flavor Parameters
    auto m_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 521));
    auto life_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FLIFE", 521));
    auto m_D = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 421));

    // Formfactor parameters
    auto g1 = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Dlnu", 1));
    auto rho2 = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Dlnu", 2));
    auto delta = std::make_shared<ParameterNode>(ParamId(ParameterType::FF, "B_Dlnu", 3));

    // Operator nodes
    auto qcd = std::make_shared<OperatorNode>("qcd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return 0; });
    qcd->addChildren({alpha_s_MZ, M_Z, mt_pole, mb_mb, m_u, m_d, m_s, m_c});
    auto m_b_muh = std::make_shared<OperatorNode>("m_b", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return QCDHelper::msbar_mass(5, winfo.hadronic_scale); });
    m_b_muh->addChildren({qcd});
    auto m_c_muh = std::make_shared<OperatorNode>("m_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return QCDHelper::msbar_mass(4, winfo.hadronic_scale); });
    m_c_muh->addChildren({qcd});
    auto m_bc = std::make_shared<OperatorNode>("m_b - m_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] - values[1]; });
    m_bc->addChildren({m_b_muh, m_c_muh});
    auto r_D = std::make_shared<OperatorNode>("r_D", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    r_D->addChildren({m_D, m_B});
    auto r_tau = std::make_shared<OperatorNode>("r_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return std::pow(values[0] / values[1], 2); });
    r_tau->addChildren({m_tau, m_B});
    auto r_e = std::make_shared<OperatorNode>("r_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return std::pow(values[0] / values[1], 2); });
    r_e->addChildren({m_e, m_B});
    auto delta2 = std::make_shared<OperatorNode>("Delta(w)^2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return std::pow(values[0], 2); });
    delta2->addChildren({delta});
    auto w_tau = std::make_shared<OperatorNode>("w_max_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->w_max(values[0], values[1]); });
    w_tau->addChildren({r_D, r_tau});
    auto w_e = std::make_shared<OperatorNode>("w_max_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->w_max(values[0], values[1]); });
    w_e->addChildren({r_D, r_e});
    auto i1_tau = std::make_shared<OperatorNode>("I1_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->I1(values[0], values[1], values[2], values[3]); });
    i1_tau->addChildren({r_D, r_tau, rho2, w_tau});
    auto i1_e = std::make_shared<OperatorNode>("I1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->I1(values[0], values[1], values[2], values[3]); });
    i1_e->addChildren({r_D, r_e, rho2, w_e});
    auto i2_tau = std::make_shared<OperatorNode>("I2_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->I2(values[0], values[1], values[2], values[3]); });
    i2_tau->addChildren({r_D, r_tau, rho2, w_tau});
    auto i2_e = std::make_shared<OperatorNode>("I2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->I2(values[0], values[1], values[2], values[3]); });
    i2_e->addChildren({r_D, r_e, rho2, w_e});
    auto i3_tau = std::make_shared<OperatorNode>("I3_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->I3(values[0], values[1], values[2], values[3]); });
    i3_tau->addChildren({r_D, r_tau, rho2, w_tau});
    auto i3_e = std::make_shared<OperatorNode>("I3_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->I3(values[0], values[1], values[2], values[3]); });
    i3_e->addChildren({r_D, r_e, rho2, w_e});
    auto i4_tau = std::make_shared<OperatorNode>("I4_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->I4(values[0], values[1], values[2], values[3]); });
    i4_tau->addChildren({r_D, r_tau, rho2, w_tau});
    auto i4_e = std::make_shared<OperatorNode>("I4_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->I4(values[0], values[1], values[2], values[3]); });
    i4_e->addChildren({r_D, r_e, rho2, w_e});
    auto w1tau = std::make_shared<OperatorNode>("w1tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W1(values[0]); });
    w1tau->addChildren({i1_tau});
    auto w1e = std::make_shared<OperatorNode>("w1e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W1(values[0]); });
    w1e->addChildren({i1_e});
    auto w2tau = std::make_shared<OperatorNode>("w2tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W2(values[0], values[1], values[2]); });
    w2tau->addChildren({i2_tau, delta2, r_tau});
    auto w2e = std::make_shared<OperatorNode>("w2e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W2(values[0], values[1], values[2]); });
    w2e->addChildren({i2_e, delta2, r_e});
    auto w3tau = std::make_shared<OperatorNode>("w3tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W3(values[0], values[1], values[2], values[3]); });
    w3tau->addChildren({i3_tau, delta2, m_tau, m_bc});
    auto w3e = std::make_shared<OperatorNode>("w3e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W3(values[0], values[1], values[2], values[3]); });
    w3e->addChildren({i3_e, delta2, m_e, m_bc});
    auto w4tau = std::make_shared<OperatorNode>("w4tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W4(values[0], values[1], values[2], values[3]); });
    w4tau->addChildren({i4_tau, delta2, m_B, m_bc});
    auto w4e = std::make_shared<OperatorNode>("w4e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W4(values[0], values[1], values[2], values[3]); });
    w4e->addChildren({i4_e, delta2, m_B, m_bc});
    auto wtau = std::make_shared<OperatorNode>("w_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W(values[0], values[1], values[2], values[3]); });
    wtau->addChildren({w1tau, w2tau, w3tau, w4tau});
    auto we = std::make_shared<OperatorNode>("w_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->W(values[0], values[1], values[2], values[3]); });
    we->addChildren({w1e, w2e, w3e, w4e});
    auto ckm = std::make_shared<OperatorNode>("ckm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->ckm(values[0], values[1]); });
    ckm->addChildren({V_cb_r, V_cb_i});
    auto prefactor = std::make_shared<OperatorNode>("prefactor", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->pref(values[0], values[1], values[2], values[3], values[4]); });
    prefactor->addChildren({G_F, life_B, m_B, m_D, g1});
    auto br_tau = std::make_shared<OperatorNode>("BR(B > D0 tau nu_tau)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_B_Dlnu(values[0], values[1], values[2]); });
    br_tau->addChildren({prefactor, ckm, wtau});
    roots.emplace(Observables::BR_B__D_TAU_NU, br_tau);

    auto br_e = std::make_shared<OperatorNode>("BR(B > D0 e nu_e)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return this->BR_B_Dlnu(values[0], values[1], values[2]); });
    br_e->addChildren({prefactor, ckm, we});
    auto xi = std::make_shared<OperatorNode>("xi__D_l_nu", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    xi->addChildren({br_tau, br_e});
    roots.emplace(Observables::XI__D_L_NU, xi);

}
