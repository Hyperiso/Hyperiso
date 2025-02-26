#include "BDstarlnuDecay.h"

double BDstarlnuDecay::ckm(double V_cb_r, double V_cb_i) {
    return std::pow(std::abs(complex_t(V_cb_r, V_cb_i)), 2);
}

double BDstarlnuDecay::pref(double G_F, double tau_B, double m_B, double m_D, double h_A1_1) {
    return std::pow(G_F * m_B * m_B * h_A1_1, 2) * m_D * tau_B / (96 * PI3 * HBAR);
}

double BDstarlnuDecay::t(double rD, double w) {
    return 1 + rD * (rD - 2 * w);
}

double BDstarlnuDecay::lambda_D(double rD, double w) {
    return 4 * rD * rD * (w * w - 1);
}

double BDstarlnuDecay::x_l(double rl, double rD, double w) {
    return rl * rl / t(rD, w);
}

double BDstarlnuDecay::phi(double rl, double rD, double w, double rho_D2) {
    return t(rD, w) * std::sqrt(lambda_D(rD, w)) * std::pow((1 - x_l(rl, rD, w)) * h_A1(w, rho_D2), 2);
}

double BDstarlnuDecay::w_max(double rD, double rl) {
    return (1 + rD * rD - rl * rl) / (2 * rD);
}

double BDstarlnuDecay::h_A1(double w, double rho_D2) {
    double z = (std::sqrt(1 + w) - RT2) / (std::sqrt(1 + w) + RT2);
    return 1 + z * (-8 * rho_D2 + z * ((53 * rho_D2 - 15) - z * (231 * rho_D2 - 91)));
}

double BDstarlnuDecay::R_1(double w, double R_11) {
    double u = w - 1;
    return R_11 + u * (-0.12 + 0.05 * u);
}

double BDstarlnuDecay::R_2(double w, double R_21) {
    double u = w - 1;
    return R_21 + u * (0.11 - 0.06 * u);
}

double BDstarlnuDecay::R_3(double w) {
    double u = w - 1;
    return 1.22 + u * (-0.052 + 0.026 * u);
}

double BDstarlnuDecay::H_Vp(double w, double rt_rD, double R_11) {
    return rt_rD * (w + 1 - R_1(w, R_11) * std::sqrt(w * w - 1));
}

double BDstarlnuDecay::H_Vm(double w, double rt_rD, double R_11) {
    return rt_rD * (w + 1 + R_1(w, R_11) * std::sqrt(w * w - 1));
}

double BDstarlnuDecay::H_V0(double w, double rD, double rt_rD, double R_21) {
    return rt_rD / std::sqrt(t(rD, w)) * ((rD - w) * (1 + w) + R_2(w, R_21) * (w * w - 1));
}

double BDstarlnuDecay::H_Vt(double w, double rD, double rt_rD, double one_m_rD_sq, double R_21) {
    return std::sqrt((w * w - 1) / t(rD, w)) / (2 * rt_rD) * (-2 * rD * (1 + w) + one_m_rD_sq * R_2(w, R_21) - t(rD, w) * R_3(w));
}

double BDstarlnuDecay::H_S(double w, double rD, double rt_rD, double one_m_rD_sq, double rqp, double R_21) {
    return std::sqrt(w * w - 1) / (2 * rt_rD * rqp) * (-2 * rD * (1 + w) + one_m_rD_sq * R_2(w, R_21) - t(rD, w) * R_3(w));
}

double BDstarlnuDecay::H_Tp(double w, double rD, double rt_rD, double rqp, double rqm, double R_11) {
    return rt_rD / std::sqrt(t(rD, w)) * (rqp * std::sqrt(w * w - 1) * R_1(w, R_11) + rqm * (1 + w));
}

double BDstarlnuDecay::H_Tm(double w, double rD, double rt_rD, double rqp, double rqm, double R_11) {
    return rt_rD / std::sqrt(t(rD, w)) * (rqp * std::sqrt(w * w - 1) * R_1(w, R_11) - rqm * (1 + w));
}

double BDstarlnuDecay::H_T0(double w, double rD, double rt_rD, double one_m_rD_sq, double rqp, double rqm, double R_11) {
    return rt_rD * (1 + w) / (one_m_rD_sq * t(rD, w)) * (2 * rqp * one_m_rD_sq * (w - 1) * R_1(w, R_11)
                                                            - rqm * (w - 1) * t(rD, w) * (R_3(w) - 1)
                                                            - rqm * (1 + rD) * (std::pow(1 - rD, 2) + 2 * (w - 1)));
}

double BDstarlnuDecay::F_V0_1(double rD, double rt_rD, double rl, double rho_D2, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rho_D2, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * std::pow(H_V0(w, rD, rt_rD, R_21), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_V0_2(double rD, double rt_rD, double rl, double rho_D2, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rho_D2, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * std::pow(H_V0(w, rD, rt_rD, R_21), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vp_1(double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::pow(H_Vp(w, rt_rD, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vp_2(double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * std::pow(H_Vp(w, rt_rD, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vm_1(double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::pow(H_Vm(w, rt_rD, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vm_2(double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * std::pow(H_Vm(w, rt_rD, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vt(double rD, double rt_rD, double one_m_rD_sq, double rl,double rho_D2, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rho_D2, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * std::pow(H_Vt(w, rD, rt_rD, one_m_rD_sq, R_21), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_S(double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rho_D2, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rqp, rho_D2, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * std::pow(H_S(w, rD, rt_rD, one_m_rD_sq, rqp, R_21), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_T0_1(double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::pow(H_T0(w, rD, rt_rD, one_m_rD_sq, rqp, rqm, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_T0_2(double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * std::pow(H_T0(w, rD, rt_rD, one_m_rD_sq, rqp, rqm, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Tp_1(double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::pow(H_Tp(w, rD, rt_rD, rqp, rqm, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Tp_2(double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * std::pow(H_Tp(w, rD, rt_rD, rqp, rqm, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Tm_1(double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::pow(H_Tm(w, rD, rt_rD, rqp, rqm, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Tm_2(double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * std::pow(H_Tm(w, rD, rt_rD, rqp, rqm, R_11), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vp_Vm_1(double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * H_Vp(w, rt_rD, R_11) * H_Vm(w, rt_rD, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vp_Vm_2(double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * H_Vp(w, rt_rD, R_11) * H_Vm(w, rt_rD, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vp_Tp(double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_Vp(w, rt_rD, R_11) * H_Tp(w, rD, rt_rD, rqp, rqm, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vp_Tm(double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_Vp(w, rt_rD, R_11) * H_Tm(w, rD, rt_rD, rqp, rqm, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vm_Tp(double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_Vm(w, rt_rD, R_11) * H_Tp(w, rD, rt_rD, rqp, rqm, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vm_Tm(double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, rl, rqp, rqm, rho_D2, R_11] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_Vm(w, rt_rD, R_11) * H_Tm(w, rD, rt_rD, rqp, rqm, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_V0_Vt(double rD, double rt_rD, double one_m_rD_sq, double rl, double rho_D2, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rho_D2, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * x_l(rl, rD, w) * H_V0(w, rD, rt_rD, R_21) * H_Vt(w, rD, rt_rD, one_m_rD_sq, R_21);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_V0_S(double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rho_D2, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rqp, rho_D2, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_V0(w, rD, rt_rD, R_21) * H_S(w, rD, rt_rD, one_m_rD_sq, rqp, R_21);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_V0_T0(double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rqp, rqm, rho_D2, R_11, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_V0(w, rD, rt_rD, R_21) * H_T0(w, rD, rt_rD, one_m_rD_sq, rqp, rqm, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vt_S(double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rho_D2, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rqp, rho_D2, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_Vt(w, rD, rt_rD, one_m_rD_sq, R_21) * H_S(w, rD, rt_rD, one_m_rD_sq, rqp, R_21);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vt_T0(double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rqp, rqm, rho_D2, R_11, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_Vt(w, rD, rt_rD, one_m_rD_sq, R_21) * H_T0(w, rD, rt_rD, one_m_rD_sq, rqp, rqm, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_S_T0(double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double R_21, double w_m, bool flag) {
    if (!flag) return 0;
    
    auto f = [this, rD, rt_rD, one_m_rD_sq, rl, rqp, rqm, rho_D2, R_11, R_21] (double w) {
        return phi(rl, rD, w, rho_D2) * std::sqrt(x_l(rl, rD, w)) * H_S(w, rD, rt_rD, one_m_rD_sq, rqp, R_21) * H_T0(w, rD, rt_rD, one_m_rD_sq, rqp, rqm, R_11);
    };

    return integrate(f, 1, w_m, 1e-3);
}

complex_t BDstarlnuDecay::C_V1() {
    auto wil = get_wilsons();
    return wil->getFM(WGroup::BCLNU, WCoef::C_V1, QCDOrder::LO);
}

complex_t BDstarlnuDecay::C_V2() {
    auto wil = get_wilsons();
    return wil->getFM(WGroup::BCLNU, WCoef::C_V2, QCDOrder::LO);
}

complex_t BDstarlnuDecay::C_A() {
    auto wil = get_wilsons();
    return wil->getFM(WGroup::BCLNU, WCoef::C_V1, QCDOrder::LO) - wil->getFM(WGroup::BCLNU, WCoef::C_V2, QCDOrder::LO);
}

complex_t BDstarlnuDecay::C_P() {
    auto wil = get_wilsons();
    return wil->getFM(WGroup::BCLNU, WCoef::C_S1, QCDOrder::LO) - wil->getFM(WGroup::BCLNU, WCoef::C_S2, QCDOrder::LO);
}

complex_t BDstarlnuDecay::C_T() {
    auto wil = get_wilsons();
    return wil->getFM(WGroup::BCLNU, WCoef::C_T, QCDOrder::LO);
}

double BDstarlnuDecay::c_flag(complex_t C) {
    if (fpeq(std::abs(C), 0.)) {
        return 0;
    }
    return 1;
}

double BDstarlnuDecay::Gamma_tau_m(double F_Vp_1, double F_Vm_1, double F_V0_1, double F_Tp_2, double F_Tm_2, double F_T0_2, 
                                   double G_Vp_Vm_1, double G_T0_V0, double G_Tp_Vp, double G_Tm_Vm, double G_Tp_Vm, double G_Tm_Vp,
                                   complex_t C_V1, complex_t C_V2, complex_t C_T)
{
    double c_vsq = std::pow(std::abs(C_V1), 2) + std::pow(std::abs(C_V2), 2);
    double c_vv = -2 * std::real(C_V1 * std::conj(C_V2));
    double c_tt = 16 * std::pow(std::abs(C_T), 2);
    double c_v1t = -8 * std::real(C_V1 * std::conj(C_T));
    double c_v2t = 8 * std::real(C_V2 * std::conj(C_T));

    return c_vsq * (F_Vp_1 + F_Vm_1 + F_V0_1) 
            + c_vv * (F_V0_1 + 2 * G_Vp_Vm_1) 
            + c_tt * (F_Tp_2 + F_Tm_2 + F_T0_2) 
            + c_v1t * (G_T0_V0 + G_Tp_Vp - G_Tm_Vm)
            + c_v2t * (G_T0_V0 + G_Tp_Vm - G_Tm_Vp);
}

double BDstarlnuDecay::Gamma_tau_p(double F_Vp_2, double F_Vm_2, double F_V0_2, double F_Vt, double F_S, double F_Tp_1, double F_Tm_1, double F_T0_1,
                                   double G_Vp_Vm_2,double G_Vt_S, double G_T0_V0, double G_Tp_Vp, double G_Tm_Vm, double G_Tp_Vm, double G_Tm_Vp,
                                   complex_t C_V1, complex_t C_V2, complex_t C_A, complex_t C_P, complex_t C_T)
{
    double c_vsq = 0.5 * std::pow(std::abs(C_V1), 2) + std::pow(std::abs(C_V2), 2);
    double c_vv = -std::real(C_V1 * std::conj(C_V2));
    double c_ss = 1.5 * std::pow(std::abs(C_P), 2);
    double c_tt = 8 * std::pow(std::abs(C_T), 2);
    double c_vs = 3 * std::real(C_A * std::conj(C_P));
    double c_v1t = -4 * std::real(C_V1 * std::conj(C_T));
    double c_v2t = 4 * std::real(C_V2 * std::conj(C_T));

    return c_vsq * (F_Vp_2 + F_Vm_2 + F_V0_2 + 3 * F_Vt) 
            + c_vv * (F_V0_2 + 2 * G_Vp_Vm_2 + 3 * F_Vt) 
            + c_ss * F_S 
            + c_tt * (F_Tp_1 + F_Tm_1 + F_T0_1) 
            + c_vs * G_Vt_S
            + c_v1t * (G_T0_V0 + G_Tp_Vp - G_Tm_Vm)
            + c_v2t * (G_T0_V0 + G_Tp_Vm - G_Tm_Vp);
}

double BDstarlnuDecay::Gamma_D_0(double F_V0_1, double F_V0_2, double F_Vt, double F_S, double F_T0_1, double F_T0_2,
                                 double G_Vt_S, double G_V0_T0,
                                 complex_t C_A, complex_t C_P, complex_t C_T)
{
    double c_aa = 0.5 * std::pow(std::abs(C_A), 2);
    double c_pp = 1.5 * std::pow(std::abs(C_P), 2);
    double c_tt = 8 * std::pow(std::abs(C_T), 2);
    double c_ap = 3 * std::real(C_A * std::conj(C_P));
    double c_at = -12 * std::real(C_A * std::conj(C_T));

    return c_aa * (2 * F_V0_1 + F_V0_2 + 3 * F_Vt) 
            + c_pp * F_S 
            + c_tt * (F_T0_1 + 2 * F_T0_2) 
            + c_ap * G_Vt_S 
            + c_at * G_V0_T0;
}

double BDstarlnuDecay::B_theta(double F_Vp_1, double F_Vm_1, double F_Tp_2, double F_Tm_2,
                               double G_V0_Vt,double G_V0_S, double G_Vt_T0, double G_Tp_Vp, double G_Tm_Vm, double G_Tp_Vm, double G_Tm_Vp, double G_T0_S,
                               complex_t C_V1, complex_t C_V2, complex_t C_A, complex_t C_P, complex_t C_T)
{
    double c_vv = 0.75 * std::pow(std::abs(C_V1), 2) - std::pow(std::abs(C_V2), 2);
    double c_aa = 1.5 * std::pow(std::abs(C_A), 2);
    double c_tt = 12 * std::pow(std::abs(C_T), 2);
    double c_ap = 1.5 * std::real(C_A * std::conj(C_P));
    double c_v1t = -6 * std::real(C_V1 * std::conj(C_T));
    double c_v2t = 6 * std::real(C_V2 * std::conj(C_T));
    double c_pt = -6 * std::real(C_P * std::conj(C_T));

    return c_vv * (F_Vp_1 - F_Vm_1) 
            + c_aa * G_V0_Vt 
            + c_tt * (F_Tp_2 - F_Tm_2) 
            + c_ap * G_V0_S 
            + c_v1t * (G_Vt_T0 + G_Tp_Vp + G_Tm_Vm)
            + c_v2t * (G_Vt_T0 + G_Tp_Vm + G_Tm_Vp)
            + c_pt * G_T0_S;
}

void BDstarlnuDecay::build_op_tree() {
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
    auto m_D_star = std::make_shared<ParameterNode>(ParamId(ParameterType::FLAVOR, "FMASS", 423));

    // Formfactor parameters
    auto h_A1_1 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Dslnu", 1));
    auto rho_D2 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Dslnu", 2));
    auto R_11 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Dslnu", 3));
    auto R_21 = std::make_shared<ParameterNode>(ParamId(ParameterType::DECAY, "B_Dslnu", 4));

    // Operator nodes
    auto qcd = std::make_shared<OperatorNode>("qcd", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return 0; });
    qcd->addChildren({alpha_s_MZ, M_Z, mt_pole, mb_mb, m_u, m_d, m_s, m_c});
    auto m_b_muh = std::make_shared<OperatorNode>("m_b", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return QCDHelper::msbar_mass(5, winfo.hadronic_scale); });
    m_b_muh->addChildren({qcd});
    auto m_c_muh = std::make_shared<OperatorNode>("m_c", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return QCDHelper::msbar_mass(4, winfo.hadronic_scale); });
    m_c_muh->addChildren({qcd});
    auto r_D = std::make_shared<OperatorNode>("r_D", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    r_D->addChildren({m_D_star, m_B});
    auto rt_rD = std::make_shared<OperatorNode>("sqrt(r_D)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return std::sqrt(values[0]); });
    rt_rD->addChildren({r_D});
    auto one_m_rD2 = std::make_shared<OperatorNode>("1 - r_D^2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return 1. - values[0] * values[0]; });
    one_m_rD2->addChildren({r_D});
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

    auto nC_V1 = std::make_shared<OperatorNode>("C_V1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_V1(); });
    auto nC_V2 = std::make_shared<OperatorNode>("C_V2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_V2(); });
    auto nC_A = std::make_shared<OperatorNode>("C_A", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_A(); });
    auto nC_P = std::make_shared<OperatorNode>("C_P", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_P(); });
    auto nC_T = std::make_shared<OperatorNode>("C_T", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return C_T(); });

    auto flag_V1 = std::make_shared<OperatorNode>("flag_V1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_flag(values[0]); });
    flag_V1->addChild(nC_V1);
    auto flag_V2 = std::make_shared<OperatorNode>("flag_V2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_flag(values[0]); });
    flag_V2->addChild(nC_V2);
    auto flag_A = std::make_shared<OperatorNode>("flag_A", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_flag(values[0]); });
    flag_A->addChild(nC_A);
    auto flag_P = std::make_shared<OperatorNode>("flag_P", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_flag(values[0]); });
    flag_P->addChild(nC_P);
    auto flag_T = std::make_shared<OperatorNode>("flag_T", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return c_flag(values[0]); });
    flag_T->addChild(nC_T);

    auto n_F_V0_1 = std::make_shared<OperatorNode>("F_V0_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_V0_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_V0_1->addChildren({r_D, rt_rD, r_tau, rho_D2, R_21, w_t, flag_V1, flag_V2});

    auto n_F_V0_1_e = std::make_shared<OperatorNode>("F_V0_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_V0_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_V0_1_e->addChildren({r_D, rt_rD, r_e, rho_D2, R_21, w_e, flag_V1, flag_V2});

    auto n_F_V0_2 = std::make_shared<OperatorNode>("F_V0_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_V0_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_V0_2->addChildren({r_D, rt_rD, r_tau, rho_D2, R_21, w_t, flag_V1, flag_V2});

    auto n_F_V0_2_e = std::make_shared<OperatorNode>("F_V0_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_V0_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_V0_2_e->addChildren({r_D, rt_rD, r_e, rho_D2, R_21, w_e, flag_V1, flag_V2});

    auto n_F_Vp_1 = std::make_shared<OperatorNode>("F_Vp_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vp_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_Vp_1->addChildren({r_D, rt_rD, r_tau, rho_D2, R_11, w_t, flag_V1, flag_V2});

    auto n_F_Vp_1_e = std::make_shared<OperatorNode>("F_Vp_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vp_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_Vp_1_e->addChildren({r_D, rt_rD, r_e, rho_D2, R_11, w_e, flag_V1, flag_V2});

    auto n_F_Vp_2 = std::make_shared<OperatorNode>("F_Vp_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vp_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_Vp_2->addChildren({r_D, rt_rD, r_tau, rho_D2, R_11, w_t, flag_V1, flag_V2});

    auto n_F_Vp_2_e = std::make_shared<OperatorNode>("F_Vp_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vp_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_Vp_2_e->addChildren({r_D, rt_rD, r_e, rho_D2, R_11, w_e, flag_V1, flag_V2});

    auto n_F_Vm_1 = std::make_shared<OperatorNode>("F_Vm_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vm_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_Vm_1->addChildren({r_D, rt_rD, r_tau, rho_D2, R_11, w_t, flag_V1, flag_V2});

    auto n_F_Vm_1_e = std::make_shared<OperatorNode>("F_Vm_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vm_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_Vm_1_e->addChildren({r_D, rt_rD, r_e, rho_D2, R_11, w_e, flag_V1, flag_V2});

    auto n_F_Vm_2 = std::make_shared<OperatorNode>("F_Vm_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vm_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_Vm_2->addChildren({r_D, rt_rD, r_tau, rho_D2, R_11, w_t, flag_V1, flag_V2});

    auto n_F_Vm_2_e = std::make_shared<OperatorNode>("F_Vm_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vm_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] || (bool)values[7]); });
    n_F_Vm_2_e->addChildren({r_D, rt_rD, r_e, rho_D2, R_11, w_e, flag_V1, flag_V2});

    auto n_F_Vt = std::make_shared<OperatorNode>("F_Vt", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vt(values[0], values[1], values[2], values[3], values[4], values[5], values[6], (bool)values[7] || (bool)values[8]); });
    n_F_Vt->addChildren({r_D, rt_rD, one_m_rD2, r_tau, rho_D2, R_21, w_t, flag_V1, flag_V2});

    auto n_F_Vt_e = std::make_shared<OperatorNode>("F_Vt_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Vt(values[0], values[1], values[2], values[3], values[4], values[5], values[6], (bool)values[7] || (bool)values[8]); });
    n_F_Vt_e->addChildren({r_D, rt_rD, one_m_rD2, r_e, rho_D2, R_21, w_e, flag_V1, flag_V2});

    auto n_F_S = std::make_shared<OperatorNode>("F_S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_S(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_S->addChildren({r_D, rt_rD, one_m_rD2, r_tau, r_qp, rho_D2, R_21, w_e, flag_P});

    auto n_F_S_e = std::make_shared<OperatorNode>("F_S_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_S(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_S_e->addChildren({r_D, rt_rD, one_m_rD2, r_e, r_qp, rho_D2, R_21, w_e, flag_P});

    auto n_F_T0_1 = std::make_shared<OperatorNode>("F_T0_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_T0_1(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], (bool)values[9]); });
    n_F_T0_1->addChildren({r_D, rt_rD, one_m_rD2, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_T});

    auto n_F_T0_1_e = std::make_shared<OperatorNode>("F_T0_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_T0_1(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], (bool)values[9]); });
    n_F_T0_1_e->addChildren({r_D, rt_rD, one_m_rD2, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_T});

    auto n_F_T0_2 = std::make_shared<OperatorNode>("F_T0_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_T0_2(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], (bool)values[9]); });
    n_F_T0_2->addChildren({r_D, rt_rD, one_m_rD2, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_T});

    auto n_F_T0_2_e = std::make_shared<OperatorNode>("F_T0_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_T0_2(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], (bool)values[9]); });
    n_F_T0_2_e->addChildren({r_D, rt_rD, one_m_rD2, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_T});

    auto n_F_Tp_1 = std::make_shared<OperatorNode>("F_Tp_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Tp_1(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_Tp_1->addChildren({r_D, rt_rD, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_T});

    auto n_F_Tp_1_e = std::make_shared<OperatorNode>("F_Tp_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Tp_1(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_Tp_1_e->addChildren({r_D, rt_rD, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_T});

    auto n_F_Tp_2 = std::make_shared<OperatorNode>("F_Tp_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Tp_2(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_Tp_2->addChildren({r_D, rt_rD, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_T});

    auto n_F_Tp_2_e = std::make_shared<OperatorNode>("F_Tp_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Tp_2(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_Tp_2_e->addChildren({r_D, rt_rD, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_T});

    auto n_F_Tm_1 = std::make_shared<OperatorNode>("F_Tm_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Tm_1(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_Tm_1->addChildren({r_D, rt_rD, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_T});

    auto n_F_Tm_1_e = std::make_shared<OperatorNode>("F_Tm_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Tm_1(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_Tm_1_e->addChildren({r_D, rt_rD, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_T});

    auto n_F_Tm_2 = std::make_shared<OperatorNode>("F_Tm_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Tm_2(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_Tm_2->addChildren({r_D, rt_rD, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_T});

    auto n_F_Tm_2_e = std::make_shared<OperatorNode>("F_Tm_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return F_Tm_2(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8]); });
    n_F_Tm_2_e->addChildren({r_D, rt_rD, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_T});

    auto n_G_Vp_Vm_1 = std::make_shared<OperatorNode>("G_Vp_Vm_1", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vp_Vm_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] && (bool)values[7]); });
    n_G_Vp_Vm_1->addChildren({r_D, rt_rD, r_tau, rho_D2, R_11, w_t, flag_V1, flag_V2});

    auto n_G_Vp_Vm_1_e = std::make_shared<OperatorNode>("G_Vp_Vm_1_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vp_Vm_1(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] && (bool)values[7]); });
    n_G_Vp_Vm_1_e->addChildren({r_D, rt_rD, r_e, rho_D2, R_11, w_e, flag_V1, flag_V2});

    auto n_G_Vp_Vm_2 = std::make_shared<OperatorNode>("G_Vp_Vm_2", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vp_Vm_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] && (bool)values[7]); });
    n_G_Vp_Vm_2->addChildren({r_D, rt_rD, r_tau, rho_D2, R_11, w_t, flag_V1, flag_V2});

    auto n_G_Vp_Vm_2_e = std::make_shared<OperatorNode>("G_Vp_Vm_2_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vp_Vm_2(values[0], values[1], values[2], values[3], values[4], values[5], (bool)values[6] && (bool)values[7]); });
    n_G_Vp_Vm_2_e->addChildren({r_D, rt_rD, r_e, rho_D2, R_11, w_e, flag_V1, flag_V2});

    auto n_G_Vp_Tp = std::make_shared<OperatorNode>("G_Vp_Tp", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vp_Tp(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vp_Tp->addChildren({r_D, rt_rD, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_V1, flag_T});

    auto n_G_Vp_Tp_e = std::make_shared<OperatorNode>("G_Vp_Tp_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vp_Tp(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vp_Tp_e->addChildren({r_D, rt_rD, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_V1, flag_T});

    auto n_G_Vp_Tm = std::make_shared<OperatorNode>("G_Vp_Tm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vp_Tm(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vp_Tm->addChildren({r_D, rt_rD, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_V2, flag_T});

    auto n_G_Vp_Tm_e = std::make_shared<OperatorNode>("G_Vp_Tm_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vp_Tm(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vp_Tm_e->addChildren({r_D, rt_rD, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_V2, flag_T});

    auto n_G_Vm_Tp = std::make_shared<OperatorNode>("G_Vm_Tp", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vm_Tp(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vm_Tp->addChildren({r_D, rt_rD, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_V2, flag_T});

    auto n_G_Vm_Tp_e = std::make_shared<OperatorNode>("G_Vm_Tp_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vm_Tp(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vm_Tp_e->addChildren({r_D, rt_rD, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_V2, flag_T});

    auto n_G_Vm_Tm = std::make_shared<OperatorNode>("G_Vm_Tm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vm_Tm(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vm_Tm->addChildren({r_D, rt_rD, r_tau, r_qp, r_qm, rho_D2, R_11, w_t, flag_V1, flag_T});

    auto n_G_Vm_Tm_e = std::make_shared<OperatorNode>("G_Vm_Tm_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vm_Tm(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vm_Tm_e->addChildren({r_D, rt_rD, r_e, r_qp, r_qm, rho_D2, R_11, w_e, flag_V1, flag_T});

    auto n_G_V0_Vt = std::make_shared<OperatorNode>("G_V0_Vt", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_V0_Vt(values[0], values[1], values[2], values[3], values[4], values[5], values[6], (bool)values[7]); });
    n_G_V0_Vt->addChildren({r_D, rt_rD, one_m_rD2, r_tau, rho_D2, R_21, w_t, flag_A});

    auto n_G_V0_S = std::make_shared<OperatorNode>("G_V0_S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_V0_S(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_V0_S->addChildren({r_D, rt_rD, one_m_rD2, r_tau, r_qp, rho_D2, R_21, w_t, flag_A, flag_P});

    auto n_G_V0_T0 = std::make_shared<OperatorNode>("G_V0_T0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_V0_T0(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], (bool)values[10] && (bool)values[12] || (bool)values[11] && (bool)values[12]); });
    n_G_V0_T0->addChildren({r_D, rt_rD, one_m_rD2, r_tau, r_qp, r_qm, rho_D2, R_11, R_21, w_t, flag_V1, flag_V2, flag_T});

    auto n_G_V0_T0_e = std::make_shared<OperatorNode>("G_V0_T0_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_V0_T0(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], (bool)values[10] && (bool)values[12] || (bool)values[11] && (bool)values[12]); });
    n_G_V0_T0_e->addChildren({r_D, rt_rD, one_m_rD2, r_e, r_qp, r_qm, rho_D2, R_11, R_21, w_e, flag_V1, flag_V2, flag_T});

    auto n_G_Vt_S = std::make_shared<OperatorNode>("G_Vt_S", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vt_S(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vt_S->addChildren({r_D, rt_rD, one_m_rD2, r_tau, r_qp, rho_D2, R_21, w_t, flag_A, flag_P});

    auto n_G_Vt_S_e = std::make_shared<OperatorNode>("G_Vt_S_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vt_S(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], (bool)values[8] && (bool)values[9]); });
    n_G_Vt_S_e->addChildren({r_D, rt_rD, one_m_rD2, r_e, r_qp, rho_D2, R_21, w_e, flag_A, flag_P});

    auto n_G_Vt_T0 = std::make_shared<OperatorNode>("G_Vt_T0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_Vt_T0(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], (bool)values[10] && (bool)values[12] || (bool)values[11] && (bool)values[12]); });
    n_G_Vt_T0->addChildren({r_D, rt_rD, one_m_rD2, r_tau, r_qp, r_qm, rho_D2, R_11, R_21, w_t, flag_V1, flag_V2, flag_T});

    auto n_G_S_T0 = std::make_shared<OperatorNode>("G_S_T0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return G_S_T0(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], (bool)values[10] && (bool)values[11]); });
    n_G_S_T0->addChildren({r_D, rt_rD, one_m_rD2, r_tau, r_qp, r_qm, rho_D2, R_11, R_21, w_t, flag_P, flag_T});

    auto gamma_tau_m = std::make_shared<OperatorNode>("Gamma_tau_m", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_tau_m(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14]); });
    gamma_tau_m->addChildren({n_F_Vp_1, n_F_Vm_1, n_F_V0_1, n_F_Tp_2, n_F_Tm_2, n_F_T0_2, n_G_Vp_Vm_1, n_G_V0_T0, n_G_Vp_Tp, n_G_Vm_Tm, n_G_Vm_Tp, n_G_Vp_Tm, nC_V1, nC_V2, nC_T});

    auto gamma_e_m = std::make_shared<OperatorNode>("Gamma_tau_m", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_tau_m(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14]); });
    gamma_e_m->addChildren({n_F_Vp_1_e, n_F_Vm_1_e, n_F_V0_1_e, n_F_Tp_2_e, n_F_Tm_2_e, n_F_T0_2_e, n_G_Vp_Vm_1_e, n_G_V0_T0_e, n_G_Vp_Tp_e, n_G_Vm_Tm_e, n_G_Vm_Tp_e, n_G_Vp_Tm_e, nC_V1, nC_V2, nC_T});

    auto gamma_tau_p = std::make_shared<OperatorNode>("Gamma_tau_p", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_tau_p(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15], values[16], values[17], values[18], values[19]); });
    gamma_tau_p->addChildren({n_F_Vp_2, n_F_Vm_2, n_F_V0_2, n_F_Vt, n_F_S, n_F_Tp_1, n_F_Tm_1, n_F_T0_1, n_G_Vp_Vm_2, n_G_Vt_S, n_G_V0_T0, n_G_Vp_Tp, n_G_Vm_Tm, n_G_Vm_Tp, n_G_Vp_Tm, nC_V1, nC_V2, nC_A, nC_P, nC_T});

    auto gamma_e_p = std::make_shared<OperatorNode>("Gamma_e_p", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_tau_p(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15], values[16], values[17], values[18], values[19]); });
    gamma_e_p->addChildren({n_F_Vp_2_e, n_F_Vm_2_e, n_F_V0_2_e, n_F_Vt_e, n_F_S_e, n_F_Tp_1_e, n_F_Tm_1_e, n_F_T0_1_e, n_G_Vp_Vm_2_e, n_G_Vt_S_e, n_G_V0_T0_e, n_G_Vp_Tp_e, n_G_Vm_Tm_e, n_G_Vm_Tp_e, n_G_Vp_Tm_e, nC_V1, nC_V2, nC_A, nC_P, nC_T});

    auto gamma = std::make_shared<OperatorNode>("Gamma", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    gamma->addChildren({gamma_tau_m, gamma_tau_p});

    auto gamma_e = std::make_shared<OperatorNode>("Gamma_e", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] + values[1]; });
    gamma_e->addChildren({gamma_e_m, gamma_e_p});

    auto gamma_D_0 = std::make_shared<OperatorNode>("Gamma_D_0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return Gamma_D_0(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10]); });
    gamma_D_0->addChildren({n_F_V0_1, n_F_V0_2, n_F_Vt, n_F_S, n_F_T0_1, n_F_T0_2, n_G_Vt_S, n_G_V0_T0, nC_A, nC_P, nC_T});

    auto b_theta = std::make_shared<OperatorNode>("b_theta", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return B_theta(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15], values[16]); });
    b_theta->addChildren({n_F_Vp_1, n_F_Vm_1, n_F_Tp_2, n_F_Tm_2, n_G_V0_Vt, n_G_V0_S, n_G_Vt_T0, n_G_Vp_Tp, n_G_Vm_Tm, n_G_Vm_Tp, n_G_Vp_Tm, n_G_S_T0, nC_V1, nC_V2, nC_A, nC_P, nC_T});

    auto prefactor = std::make_shared<OperatorNode>("G0", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return pref(values[0], values[1], values[2], values[3], values[4]); });
    prefactor->addChildren({G_F, life_B, m_B, m_D_star, h_A1_1});

    auto nckm = std::make_shared<OperatorNode>("ckm", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return ckm(values[0], values[1]); });
    nckm->addChildren({V_cb_r, V_cb_i});

    auto BR = std::make_shared<OperatorNode>("BR(B > D* tau nu)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] * values[1] * values[2]; });
    BR->addChildren({prefactor, nckm, gamma});

    roots.emplace(Observables::BR_B__DSTAR_TAU_NU, BR);

    auto R_D = std::make_shared<OperatorNode>("R(D*)", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    R_D->addChildren({gamma, gamma_e});

    roots.emplace(Observables::R_DSTAR, R_D);

    auto A_FB = std::make_shared<OperatorNode>("A_FB", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    A_FB->addChildren({b_theta, gamma});

    roots.emplace(Observables::A_FB_B__DSTAR_TAU_NU, A_FB);

    auto P_tau = std::make_shared<OperatorNode>("P_tau", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return (values[0] - values[1]) / values[2]; });
    P_tau->addChildren({gamma_tau_p, gamma_tau_m, gamma});

    roots.emplace(Observables::P_TAU_B__DSTAR_TAU_NU, P_tau);

    auto P_D = std::make_shared<OperatorNode>("P_D", [this] ([[maybe_unused]] const std::vector<scalar_t>& values) { return values[0] / values[1]; });
    P_D->addChildren({gamma_D_0, gamma});

    roots.emplace(Observables::P_D_B__DSTAR_TAU_NU, P_D);
}
