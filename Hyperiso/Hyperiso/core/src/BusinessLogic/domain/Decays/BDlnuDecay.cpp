#include "BDlnuDecay.h"

double BDlnuDecay::ckm(complex_t V_cb) {
    return std::pow(std::abs(V_cb), 2);
}

double BDlnuDecay::pref(double G_F, double tau_B, double m_B, double m_D, double V11) {
    return std::pow(G_F * m_B * m_B * V11, 2) * m_D * tau_B / (96 * PI3 * HBAR);
}

double BDlnuDecay::t(double rD, double w) {
    return 1 + rD * (rD - 2 * w);
}

double BDlnuDecay::lambda_D(double rD, double w) {
    return 4 * rD * rD * (w * w - 1);
}

double BDlnuDecay::x_l(double rl, double rD, double w) {
    return rl * rl / t(rD, w);
}

double BDlnuDecay::phi(double rl, double rD, double w) {
    return t(rD, w) * std::sqrt(lambda_D(rD, w)) * std::pow(1 - x_l(rl, rD, w), 2);
}

double BDlnuDecay::w_max(double rD, double rl) {
    return (1 + rD * rD - rl * rl) / (2 * rD);
}

double BDlnuDecay::V_1(double w, double rho_D2) {
    double z = (std::sqrt(1 + w) - RT2) / (std::sqrt(1 + w) + RT2);
    return 1 + z * (-8 * rho_D2 + z * ((51 * rho_D2 - 10) - z * (252 * rho_D2 - 84)));
}

double BDlnuDecay::S_1(double w, double rho_D2, double Delta) {
    double u = w - 1;
    return V_1(w, rho_D2) * (1 + Delta * (-0.019 + u * (0.041 - 0.015 * u)));
}

double BDlnuDecay::H_V0(double w, double rD, double rho_D2) {
    return std::sqrt(rD * (w * w - 1) / t(rD, w)) * (1 + rD) * V_1(w, rho_D2);
}

double BDlnuDecay::H_Vt(double w, double rD, double rho_D2, double Delta) {
    return std::sqrt(rD / t(rD, w)) * (1 - rD) * (1 + w) * S_1(w, rho_D2, Delta);
}

double BDlnuDecay::H_S(double w, double rD, double rqm, double rho_D2, double Delta) {
    return std::sqrt(rD) * (1 - rD) * (1 + w) * S_1(w, rho_D2, Delta) / rqm;
}

double BDlnuDecay::H_T(double w, double rD, double rqp, double rho_D2, double Delta) {
    double a = std::sqrt(rD * (w * w - 1)) * rqp / (t(rD, w) * (1 + rD));
    return -a * (std::pow(1 + rD, 2) * V_1(w, rho_D2) - 2 * rD * (1 + w) * S_1(w, rho_D2, Delta));
}

double BDlnuDecay::F_V0_1(double rD, double rl, double rho_D2, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rho_D2] (double w) {
        return phi(rl, rD, w) * std::pow(H_V0(w, rD, rho_D2), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_V0_2(double rD, double rl, double rho_D2, double w_m, bool flag) {
    if (!flag) return 0;

    auto f = [this, rD, rl, rho_D2] (double w) {
        return phi(rl, rD, w) * x_l(rl, rD, w) * std::pow(H_V0(w, rD, rho_D2), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_Vt(double rD, double rl, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;

    auto f = [this, rD, rl, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * x_l(rl, rD, w) * std::pow(H_Vt(w, rD, rho_D2, Delta), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_S(double rD, double rl, double rqm, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rqm, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * std::pow(H_S(w, rD, rqm, rho_D2, Delta), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_T_1(double rD, double rl, double rqp, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rqp, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * std::pow(H_T(w, rD, rqp, rho_D2, Delta), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_T_2(double rD, double rl, double rqp, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rqp, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * x_l(rl, rD, w) * std::pow(H_T(w, rD, rqp, rho_D2, Delta), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_V0_Vt(double rD, double rl, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * x_l(rl, rD, w) * H_V0(w, rD, rho_D2) * H_Vt(w, rD, rho_D2, Delta);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_V0_S(double rD, double rl, double rqm, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rqm, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * std::sqrt(x_l(rl, rD, w)) * H_V0(w, rD, rho_D2) * H_S(w, rD, rqm, rho_D2, Delta);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_V0_T(double rD, double rl, double rqp, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rqp, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * std::sqrt(x_l(rl, rD, w)) * H_V0(w, rD, rho_D2) * H_T(w, rD, rqp, rho_D2, Delta);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_Vt_S(double rD, double rl, double rqm, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rqm, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * std::sqrt(x_l(rl, rD, w)) * H_Vt(w, rD, rho_D2, Delta) * H_S(w, rD, rqm, rho_D2, Delta);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_Vt_T(double rD, double rl, double rqp, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rqp, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * std::sqrt(x_l(rl, rD, w)) * H_Vt(w, rD, rho_D2, Delta) * H_T(w, rD, rqp, rho_D2, Delta);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_S_T(double rD, double rl, double rqp, double rqm, double rho_D2, double Delta, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rl, rqp, rqm, rho_D2, Delta] (double w) {
        return phi(rl, rD, w) * H_S(w, rD, rqm, rho_D2, Delta) * H_T(w, rD, rqp, rho_D2, Delta);
    };

    return integrate(f, 1, w_m, 1e-3);
}

complex_t BDlnuDecay::C_V() {
    return w_proxy->getFM(WGroup::BCLNU, WCoef::C_V1, QCDOrder::LO) + w_proxy->getFM(WGroup::BCLNU, WCoef::C_V2, QCDOrder::LO);
}

complex_t BDlnuDecay::C_S() {
    return w_proxy->getFM(WGroup::BCLNU, WCoef::C_S1, QCDOrder::LO) + w_proxy->getFM(WGroup::BCLNU, WCoef::C_S2, QCDOrder::LO);
}

complex_t BDlnuDecay::C_T() {
    return w_proxy->getFM(WGroup::BCLNU, WCoef::C_T, QCDOrder::LO);
}

double BDlnuDecay::c_flag(complex_t C) {
    if (fpeq(std::abs(C), 0.)) {
        return 0;
    }
    return 1;
}

double BDlnuDecay::Gamma_m(double F_V0_1, double F_T_2, double G_V0_T, complex_t C_V, complex_t C_T) {
    double c_vv = std::pow(std::abs(C_V), 2);
    double c_tt = 16 * std::pow(std::abs(C_T), 2);
    double c_vt = -8 * std::real(C_V * std::conj(C_T));

    return c_vv * F_V0_1 + c_tt * F_T_2 + c_vt * G_V0_T;
}

double BDlnuDecay::Gamma_p(double F_V0_2, double F_Vt, double F_S, double F_T_1, double G_Vt_S, double G_V0_T,
                           complex_t C_V, complex_t C_S, complex_t C_T)
{
    double c_vv = 0.5 * std::pow(std::abs(C_V), 2);
    double c_ss = 1.5 * std::pow(std::abs(C_S), 2);
    double c_tt = 8 * std::pow(std::abs(C_T), 2);
    double c_vs = 3 * std::real(C_V * std::conj(C_S));
    double c_vt = -4 * std::real(C_V * std::conj(C_T));

    return c_vv * (F_V0_2 + 3 * F_Vt) + c_ss * F_S + c_tt * F_T_1 + c_vs * G_Vt_S + c_vt * G_V0_T;
}

double BDlnuDecay::Gamma(double gamma_p, double gamma_m) {
    return gamma_p + gamma_m;
}

double BDlnuDecay::B_theta(double G_V0_Vt, double G_V0_S, double G_Vt_T, double G_S_T, complex_t C_V, complex_t C_S, complex_t C_T) {
    double c_vv = 1.5 * std::pow(std::abs(C_V), 2);
    double c_vs = 1.5 * std::real(C_V * std::conj(C_S));
    double c_vt = -6 * std::real(C_V * std::conj(C_T));
    double c_st = -6 * std::real(C_S * std::conj(C_T));

    return c_vv * G_V0_Vt + c_vs * G_V0_S + c_vt * G_Vt_T + c_st * G_S_T;
}

double BDlnuDecay::BR_B_Dtaunu(double pref, double ckm, double width) {
    return pref * ckm * width;
}

void BDlnuDecay::build_op_tree() {
    // SM Parameters
    auto G_F = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "SMINPUTS", 2));
    auto m_tau = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 15));
    auto m_e = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 11));
    auto m_d = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 1));
    auto m_u = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 2));
    auto m_s = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 3));
    auto m_c = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "MASS", 4));
    auto V_cb = std::make_shared<ParameterNode>(ParamId(ParameterType::SM, "VCKM", LhaID(1, 2)));

    // Flavor Parameters
    auto m_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 521));
    auto life_B = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FLIFE", 521));
    auto m_D = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 421));

    // Formfactor parameters
    auto V11 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Dlnu", 1));
    auto rho_D2 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Dlnu", 2));
    auto Delta = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Dlnu", 3));

    // Scales and Wilsons
    auto hadronic_scale = std::make_shared<ParameterNode>(ParamId(ParameterType::WILSON, "B_SCALE", 1));

    // Operator nodes
    auto m_b_muh = std::make_shared<OperatorNode>("m_b", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return QCDHelper::msbar_mass(5, values[0]); });
    m_b_muh->addChild(hadronic_scale);
    auto m_c_muh = std::make_shared<OperatorNode>("m_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return QCDHelper::msbar_mass(4, values[0]); });
    m_c_muh->addChild(hadronic_scale);

    auto r_D = std::make_shared<OperatorNode>("r_D", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    r_D->addChildren({m_D, m_B});
    auto r_tau = std::make_shared<OperatorNode>("r_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    r_tau->addChildren({m_tau, m_B});
    auto r_e = std::make_shared<OperatorNode>("r_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    r_e->addChildren({m_e, m_B});
    auto r_qp = std::make_shared<OperatorNode>("r_qp", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return (values[0] + values[1]) / values[2]; });
    r_qp->addChildren({m_b_muh, m_c_muh, m_B});
    auto r_qm = std::make_shared<OperatorNode>("r_qm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return (values[0] - values[1]) / values[2]; });
    r_qm->addChildren({m_b_muh, m_c_muh, m_B});
    auto w_t = std::make_shared<OperatorNode>("w_t", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return w_max(values[0], values[1]); });
    w_t->addChildren({r_D, r_tau});
    auto w_e = std::make_shared<OperatorNode>("w_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return w_max(values[0], values[1]); });
    w_e->addChildren({r_D, r_e});

    auto nC_V = std::make_shared<OperatorNode>("C_V", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_V(); });
    auto nC_S = std::make_shared<OperatorNode>("C_S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_S(); });
    auto nC_T = std::make_shared<OperatorNode>("C_T", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_T(); });

    auto flag_V = std::make_shared<OperatorNode>("flag_V", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_flag(values[0]); });
    flag_V->addChild(nC_V);
    auto flag_S = std::make_shared<OperatorNode>("flag_S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_flag(values[0]); });
    flag_S->addChild(nC_S);
    auto flag_T = std::make_shared<OperatorNode>("flag_T", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_flag(values[0]); });
    flag_T->addChild(nC_T);

    auto n_F_V0_1 = std::make_shared<OperatorNode>("F_V0_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_V0_1(values[0], values[1], values[2], values[3], (bool)values[4]); });
    n_F_V0_1->addChildren({r_D, r_tau, rho_D2, w_t, flag_V});

    auto n_F_V0_1_e = std::make_shared<OperatorNode>("F_V0_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_V0_1(values[0], values[1], values[2], values[3], (bool)values[4]); });
    n_F_V0_1_e->addChildren({r_D, r_e, rho_D2, w_e, flag_V});

    auto n_F_V0_2 = std::make_shared<OperatorNode>("F_V0_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_V0_2(values[0], values[1], values[2], values[3], (bool)values[4]); });
    n_F_V0_2->addChildren({r_D, r_tau, rho_D2, w_t, flag_V});

    auto n_F_V0_2_e = std::make_shared<OperatorNode>("F_V0_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_V0_2(values[0], values[1], values[2], values[3], (bool)values[4]); });
    n_F_V0_2_e->addChildren({r_D, r_e, rho_D2, w_e, flag_V});

    auto n_F_Vt = std::make_shared<OperatorNode>("F_Vt", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vt(values[0], values[1], values[2], values[3], values[4], (bool)values[5]); });
    n_F_Vt->addChildren({r_D, r_tau, rho_D2, Delta, w_t, flag_V});

    auto n_F_Vt_e = std::make_shared<OperatorNode>("F_Vt_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vt(values[0], values[1], values[2], values[3], values[4], (bool)values[5]); });
    n_F_Vt_e->addChildren({r_D, r_e, rho_D2, Delta, w_e, flag_V});

    auto n_F_S = std::make_shared<OperatorNode>("F_S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_S(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6]); });
    n_F_S->addChildren({r_D, r_tau, r_qm, rho_D2, Delta, w_t, flag_S});

    auto n_F_S_e = std::make_shared<OperatorNode>("F_S_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_S(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6]); });
    n_F_S_e->addChildren({r_D, r_e, r_qm, rho_D2, Delta, w_e, flag_S});

    auto n_F_T_1 = std::make_shared<OperatorNode>("F_T_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_T_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6]); });
    n_F_T_1->addChildren({r_D, r_tau, r_qp, rho_D2, Delta, w_t, flag_T});

    auto n_F_T_1_e = std::make_shared<OperatorNode>("F_T_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_T_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6]); });
    n_F_T_1_e->addChildren({r_D, r_e, r_qp, rho_D2, Delta, w_e, flag_T});

    auto n_F_T_2 = std::make_shared<OperatorNode>("F_T_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_T_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6]); });
    n_F_T_2->addChildren({r_D, r_tau, r_qp, rho_D2, Delta, w_t, flag_T});

    auto n_F_T_2_e = std::make_shared<OperatorNode>("F_T_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_T_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6]); });
    n_F_T_2_e->addChildren({r_D, r_e, r_qp, rho_D2, Delta, w_e, flag_T});

    auto n_G_V0_Vt = std::make_shared<OperatorNode>("G_V0_Vt", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_V0_Vt(values[0], values[1], values[2], values[3], values[4], (bool)values[5]); });
    n_G_V0_Vt->addChildren({r_D, r_tau, rho_D2, Delta, w_t, flag_V});

    auto n_G_V0_S = std::make_shared<OperatorNode>("G_V0_S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_V0_S(values[0], values[1], values[2], values[3], values[4], values[5], (bool)(values[6]) && (bool)values[7]); });
    n_G_V0_S->addChildren({r_D, r_tau, r_qm, rho_D2, Delta, w_t, flag_V, flag_S});

    auto n_G_V0_T = std::make_shared<OperatorNode>("G_V0_T", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_V0_T(values[0], values[1], values[2], values[3], values[4], values[5], (bool)(values[6]) && (bool)values[7]); });
    n_G_V0_T->addChildren({r_D, r_tau, r_qp, rho_D2, Delta, w_t, flag_V, flag_T});

    auto n_G_V0_T_e = std::make_shared<OperatorNode>("G_V0_T_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_V0_T(values[0], values[1], values[2], values[3], values[4], values[5], (bool)(values[6]) && (bool)values[7]); });
    n_G_V0_T_e->addChildren({r_D, r_e, r_qp, rho_D2, Delta, w_e, flag_V, flag_T});

    auto n_G_Vt_S = std::make_shared<OperatorNode>("G_Vt_S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vt_S(values[0], values[1], values[2], values[3], values[4], values[5], (bool)(values[6]) && (bool)values[7]); });
    n_G_Vt_S->addChildren({r_D, r_tau, r_qm, rho_D2, Delta, w_t, flag_V, flag_S});

    auto n_G_Vt_S_e = std::make_shared<OperatorNode>("G_Vt_S_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vt_S(values[0], values[1], values[2], values[3], values[4], values[5], (bool)(values[6]) && (bool)values[7]); });
    n_G_Vt_S_e->addChildren({r_D, r_e, r_qm, rho_D2, Delta, w_e, flag_V, flag_S});

    auto n_G_Vt_T = std::make_shared<OperatorNode>("G_Vt_T", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vt_T(values[0], values[1], values[2], values[3], values[4], values[5], (bool)(values[6]) && (bool)values[7]); });
    n_G_Vt_T->addChildren({r_D, r_tau, r_qp, rho_D2, Delta, w_t, flag_V, flag_T});

    auto n_G_S_T = std::make_shared<OperatorNode>("G_S_T", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_S_T(values[0], values[1], values[2], values[3], values[4], values[5], values[6], (bool)(values[7]) && (bool)values[8]); });
    n_G_S_T->addChildren({r_D, r_tau, r_qp, r_qm, rho_D2, Delta, w_t, flag_S, flag_T});
    
    auto gamma_m = std::make_shared<OperatorNode>("Gamma_m", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_m(values[0], values[1], values[2], values[3], values[4]); });
    gamma_m->addChildren({n_F_V0_1, n_F_T_2, n_G_V0_T, nC_V, nC_T});

    auto gamma_m_e = std::make_shared<OperatorNode>("Gamma_m_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_m(values[0], values[1], values[2], values[3], values[4]); });
    gamma_m_e->addChildren({n_F_V0_1_e, n_F_T_2_e, n_G_V0_T_e, nC_V, nC_T});

    auto gamma_p = std::make_shared<OperatorNode>("Gamma_p", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_p(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8]); });
    gamma_p->addChildren({n_F_V0_2, n_F_Vt, n_F_S, n_F_T_1, n_G_Vt_S, n_G_V0_T, nC_V, nC_S, nC_T});

    auto gamma_p_e = std::make_shared<OperatorNode>("Gamma_p_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_p(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8]); });
    gamma_p_e->addChildren({n_F_V0_2_e, n_F_Vt_e, n_F_S_e, n_F_T_1_e, n_G_Vt_S_e, n_G_V0_T_e, nC_V, nC_S, nC_T});

    auto gamma = std::make_shared<OperatorNode>("Gamma", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    gamma->addChildren({gamma_m, gamma_p});

    auto gamma_e = std::make_shared<OperatorNode>("Gamma_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    gamma_e->addChildren({gamma_m_e, gamma_p_e});

    auto b_theta = std::make_shared<OperatorNode>("b_theta", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return B_theta(values[0], values[1], values[2], values[3], values[4], values[5], values[6]); });
    b_theta->addChildren({n_G_V0_Vt, n_G_V0_S, n_G_Vt_T, n_G_S_T, nC_V, nC_S, nC_T});

    auto prefactor = std::make_shared<OperatorNode>("G0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref(values[0], values[1], values[2], values[3], values[4]); });
    prefactor->addChildren({G_F, life_B, m_B, m_D, V11});

    auto nckm = std::make_shared<OperatorNode>("ckm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return ckm(values[0]); });
    nckm->addChildren({V_cb});

    auto BR = std::make_shared<OperatorNode>("BR(B > D tau nu)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return BR_B_Dtaunu(values[0], values[1], values[2]); });
    BR->addChildren({prefactor, nckm, gamma});

    roots.emplace(Observables::BR_B__D_TAU_NU, BR);

    auto R_D = std::make_shared<OperatorNode>("R(D)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    R_D->addChildren({gamma, gamma_e});

    roots.emplace(Observables::R_D, R_D);

    auto A_FB = std::make_shared<OperatorNode>("A_FB", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    A_FB->addChildren({b_theta, gamma});

    roots.emplace(Observables::A_FB_B__D_TAU_NU, A_FB);

    auto P_tau = std::make_shared<OperatorNode>("P_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return (values[0] - values[1]) / values[2]; });
    P_tau->addChildren({gamma_p, gamma_m, gamma});

    roots.emplace(Observables::P_TAU_B__D_TAU_NU, P_tau);

}