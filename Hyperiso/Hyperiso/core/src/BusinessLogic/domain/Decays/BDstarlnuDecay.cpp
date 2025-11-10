#include "BDstarlnuDecay.h"

void BDstarlnuDecay::load_params() {
    ObsParameterProxy p;
    double V_cb2 = std::pow(std::abs(p(ParamId{ParameterType::SM, "VCKM", {1, 2}})), 2);

    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.m_e = p(ParamId{ParameterType::SM, "MASS", 11});
    cache.m_tau = p(ParamId{ParameterType::SM, "MASS", 15});
    cache.m_B = p(ParamId{ParameterType::FLAVOR, "FMASS", 511});
    cache.m_D_star = p(ParamId{ParameterType::FLAVOR, "FMASS", 413});
    cache.tau_B = p(ParamId{ParameterType::FLAVOR, "FLIFE", 511});
    cache.h_A1_1 = p(ParamId{ParameterType::DECAY, "B_Dslnu", 1});
    cache.rho_D2 = p(ParamId{ParameterType::DECAY, "B_Dslnu", 2});
    cache.R_11 = p(ParamId{ParameterType::DECAY, "B_Dslnu", 3});
    cache.R_21 = p(ParamId{ParameterType::DECAY, "B_Dslnu", 4});
    cache.r_D = cache.m_D_star / cache.m_B;
    cache.sqrt_rD = std::sqrt(cache.r_D);
    cache.one_m_rD2 = 1. - std::pow(cache.r_D, 2);
    cache.r_e = cache.m_e / cache.m_B;
    cache.r_tau = cache.m_tau / cache.m_B;
    double m_b = ObsQCDProxy()(MassConfig(5, cache.m_B, MassType::POLE, MassType::POLE));
    double m_c = ObsQCDProxy()(MassConfig(4, cache.m_B, MassType::POLE, MassType::POLE));
    cache.r_qp = (m_b + m_c) / cache.m_B;
    cache.r_qm = (m_b - m_c) / cache.m_B;
    cache.w_e = w_max(cache.r_e);
    cache.w_tau = w_max(cache.r_tau);
    cache.BR_pref = std::pow(cache.G_F * cache.m_B * cache.m_B * cache.h_A1_1, 2) * cache.m_D_star * cache.tau_B * V_cb2 / (96 * PI3 * HBAR);
    cache.C_V1 = w_proxy->getFM(WGroup::CC_bc, WCoef::C_V1_bc, QCDOrder::LO);
    cache.C_V2 = w_proxy->getFM(WGroup::CC_bc, WCoef::C_V2_bc, QCDOrder::LO);
    cache.C_A = cache.C_V1 - cache.C_V2;
    cache.C_P = w_proxy->getFM(WGroup::CC_bc, WCoef::C_S1_bc, QCDOrder::LO) - w_proxy->getFM(WGroup::CC_bc, WCoef::C_S2_bc, QCDOrder::LO);
    cache.C_T = w_proxy->getFM(WGroup::CC_bc, WCoef::C_T_bc, QCDOrder::LO);
    cache.C_V1_flag = !fpeq(std::abs(cache.C_V1), 0.0);
    cache.C_V2_flag = !fpeq(std::abs(cache.C_V2), 0.0);
    cache.C_A_flag = !fpeq(std::abs(cache.C_A), 0.0);
    cache.C_P_flag = !fpeq(std::abs(cache.C_P), 0.0);
    cache.C_T_flag = !fpeq(std::abs(cache.C_T), 0.0);
    cache.Gamma_p = 0.0;
    cache.Gamma_m = 0.0;

    printf("m_B = %.4e\n", cache.m_B);
    printf("m_D = %.4e\n", cache.m_D_star);
    printf("hA_1(1) = %.4e\n", cache.h_A1_1);
    printf("R1(1) = %.4e\n", cache.R_11);
    printf("R2(1) = %.4e\n", cache.R_21);
    // printf("R3(1) = %.4e\n", R3_1);
}

double BDstarlnuDecay::t(double w) {
    return 1 + cache.r_D * (cache.r_D - 2 * w);
}

double BDstarlnuDecay::lambda_D(double w) {
    return 4 * cache.r_D * cache.r_D * (w * w - 1);
}

double BDstarlnuDecay::x_l(double rl, double w) {
    return rl * rl / t(w);
}

double BDstarlnuDecay::phi(double rl, double w) {
    return t(w) * std::sqrt(lambda_D(w)) * std::pow((1 - x_l(rl, w)) * h_A1(w), 2);
}

double BDstarlnuDecay::w_max(double rl) {
    return (1 + cache.r_D * cache.r_D - rl * rl) / (2 * cache.r_D);
}

double BDstarlnuDecay::h_A1(double w) {
    double z = (std::sqrt(1 + w) - RT2) / (std::sqrt(1 + w) + RT2);
    return 1 + z * (-8 * cache.rho_D2 + z * ((53 * cache.rho_D2 - 15) - z * (231 * cache.rho_D2 - 91)));
}

double BDstarlnuDecay::R_1(double w) {
    double u = w - 1;
    return cache.R_11 + u * (-0.12 + 0.05 * u);
}

double BDstarlnuDecay::R_2(double w) {
    double u = w - 1;
    return cache.R_21 + u * (-0.11 - 0.06 * u);
}

double BDstarlnuDecay::R_3(double w) {
    double u = w - 1;
    return 0.97 + u * (-0.052 + 0.026 * u);
}

double BDstarlnuDecay::H_Vp(double w) {
    return cache.sqrt_rD * (w + 1 - R_1(w) * std::sqrt(w * w - 1));
}

double BDstarlnuDecay::H_Vm(double w) {
    return cache.sqrt_rD * (w + 1 + R_1(w) * std::sqrt(w * w - 1));
}

double BDstarlnuDecay::H_V0(double w) {
    return cache.sqrt_rD / std::sqrt(t(w)) * ((cache.r_D - w) * (1 + w) + R_2(w) * (w * w - 1));
}

double BDstarlnuDecay::H_Vt(double w) {
    return std::sqrt((w * w - 1) / t(w)) / (2 * cache.sqrt_rD) * (-2 * cache.r_D * (1 + w) + cache.one_m_rD2 * R_2(w) - t(w) * R_3(w));
}

double BDstarlnuDecay::H_S(double w) {
    return std::sqrt(w * w - 1) / (2 * cache.sqrt_rD * cache.r_qp) * (-2 * cache.r_D * (1 + w) + cache.one_m_rD2 * R_2(w) - t(w) * R_3(w));
}

double BDstarlnuDecay::H_Tp(double w) {
    return cache.sqrt_rD / std::sqrt(t(w)) * (cache.r_qp * std::sqrt(w * w - 1) * R_1(w) + cache.r_qm * (1 + w));
}

double BDstarlnuDecay::H_Tm(double w) {
    return cache.sqrt_rD / std::sqrt(t(w)) * (cache.r_qp * std::sqrt(w * w - 1) * R_1(w) - cache.r_qm * (1 + w));
}

double BDstarlnuDecay::H_T0(double w) {
    return cache.sqrt_rD * (1 + w) / (cache.one_m_rD2 * t(w)) * (2 * cache.r_qp * cache.one_m_rD2 * (w - 1) * R_1(w)
                                                            - cache.r_qm * (w - 1) * t(w) * (R_3(w) - 1)
                                                            - cache.r_qm * (1 + cache.r_D) * (std::pow(1 - cache.r_D, 2) + 2 * (w - 1)));
}

double BDstarlnuDecay::F_V0_1(double rl, double w_m) {
    if (!cache.C_V1_flag && !cache.C_V2_flag) return 0;
    auto f = [this, rl] (double w) { return phi(rl, w) * std::pow(H_V0(w), 2); };
    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_V0_2(double rl, double w_m) {
    if (!cache.C_V1_flag && !cache.C_V2_flag) return 0;

    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_V0(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vp_1(double rl, double w_m) {
    if (!cache.C_V1_flag && !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_Vp(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vp_2(double rl, double w_m) {
    if (!cache.C_V1_flag && !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_Vp(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vm_1(double rl, double w_m) {
    if (!cache.C_V1_flag && !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_Vm(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vm_2(double rl, double w_m) {
    if (!cache.C_V1_flag && !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_Vm(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Vt(double rl, double w_m) {
    if (!cache.C_V1_flag && !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_Vt(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_S(double rl, double w_m) {
    if (!cache.C_P_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_S(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_T0_1(double rl, double w_m) {
    if (!cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_T0(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_T0_2(double rl, double w_m) {
    if (!cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_T0(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Tp_1(double rl, double w_m) {
    if (!cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_Tp(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Tp_2(double rl, double w_m) {
    if (!cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_Tp(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Tm_1(double rl, double w_m) {
    if (!cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_Tm(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::F_Tm_2(double rl, double w_m) {
    if (!cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_Tm(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vp_Vm_1(double rl, double w_m) {
    if (!cache.C_V1_flag || !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * H_Vp(w) * H_Vm(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vp_Vm_2(double rl, double w_m) {
    if (!cache.C_V1_flag || !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * H_Vp(w) * H_Vm(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vp_Tp(double rl, double w_m) {
    if (!cache.C_V1_flag || !cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_Vp(w) * H_Tp(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vp_Tm(double rl, double w_m) {
    if (!cache.C_T_flag || !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_Vp(w) * H_Tm(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vm_Tp(double rl, double w_m) {
    if (!cache.C_T_flag || !cache.C_V2_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_Vm(w) * H_Tp(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vm_Tm(double rl, double w_m) {
    if (!cache.C_V1_flag || !cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_Vm(w) * H_Tm(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_V0_Vt(double rl, double w_m) {
    if (!cache.C_A_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * H_V0(w) * H_Vt(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_V0_S(double rl, double w_m) {
    if (!cache.C_A_flag || !cache.C_P_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_V0(w) * H_S(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_V0_T0(double rl, double w_m) {
    if (!(cache.C_V1_flag && cache.C_T_flag || cache.C_V2_flag && cache.C_T_flag)) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_V0(w) * H_T0(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vt_S(double rl, double w_m) {
    if (!cache.C_A_flag || !cache.C_P_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_Vt(w) * H_S(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_Vt_T0(double rl, double w_m) {
    if (!(cache.C_V1_flag && cache.C_T_flag || cache.C_V2_flag && cache.C_T_flag)) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_Vt(w) * H_T0(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::G_S_T0(double rl, double w_m) {
    if (!cache.C_T_flag || !cache.C_P_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_S(w) * H_T0(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDstarlnuDecay::Gamma_tau_m(double rl, double w_m) {
    double c_vsq = std::pow(std::abs(cache.C_V1), 2) + std::pow(std::abs(cache.C_V2), 2);
    double c_vv = -2 * std::real(cache.C_V1 * std::conj(cache.C_V2));
    double c_tt = 16 * std::pow(std::abs(cache.C_T), 2);
    double c_v1t = -8 * std::real(cache.C_V1 * std::conj(cache.C_T));
    double c_v2t = 8 * std::real(cache.C_V2 * std::conj(cache.C_T));

    double f_v0_1 = F_V0_1(rl, w_m);
    double g_v0_t0 = G_V0_T0(rl, w_m);
    return c_vsq * (F_Vp_1(rl, w_m) + F_Vm_1(rl, w_m) + f_v0_1) 
            + c_vv * (f_v0_1 + 2 * G_Vp_Vm_1(rl, w_m)) 
            + c_tt * (F_Tp_2(rl, w_m) + F_Tm_2(rl, w_m) + F_T0_2(rl, w_m)) 
            + c_v1t * (g_v0_t0 + G_Vp_Tp(rl, w_m) - G_Vm_Tm(rl, w_m))
            + c_v2t * (g_v0_t0 + G_Vm_Tp(rl, w_m) - G_Vp_Tm(rl, w_m));
}

double BDstarlnuDecay::Gamma_tau_p(double rl, double w_m) {
    double c_vsq = 0.5 * std::pow(std::abs(cache.C_V1), 2) + std::pow(std::abs(cache.C_V2), 2);
    double c_vv = -std::real(cache.C_V1 * std::conj(cache.C_V2));
    double c_ss = 1.5 * std::pow(std::abs(cache.C_P), 2);
    double c_tt = 8 * std::pow(std::abs(cache.C_T), 2);
    double c_vs = 3 * std::real(cache.C_A * std::conj(cache.C_P));
    double c_v1t = -4 * std::real(cache.C_V1 * std::conj(cache.C_T));
    double c_v2t = 4 * std::real(cache.C_V2 * std::conj(cache.C_T));

    double f_v0_2 = F_V0_2(rl, w_m);
    double g_v0_t0 = G_V0_T0(rl, w_m);
    return c_vsq * (F_Vp_2(rl, w_m) + F_Vm_2(rl, w_m) + f_v0_2 + 3 * F_Vt(rl, w_m)) 
            + c_vv * (f_v0_2 + 2 * G_Vp_Vm_2(rl, w_m) + 3 * F_Vt(rl, w_m)) 
            + c_ss * F_S (rl, w_m)
            + c_tt * (F_Tp_1(rl, w_m) + F_Tm_1(rl, w_m) + F_T0_1(rl, w_m)) 
            + c_vs * G_Vt_S(rl, w_m)
            + c_v1t * (g_v0_t0 + G_Vp_Tp(rl, w_m) - G_Vm_Tm(rl, w_m))
            + c_v2t * (g_v0_t0 + G_Vm_Tp(rl, w_m) - G_Vp_Tm(rl, w_m));
}

double BDstarlnuDecay::BR() {
    if (cache.Gamma_p == 0)
        cache.Gamma_p = Gamma_tau_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = Gamma_tau_m(cache.r_tau, cache.w_tau);

    return cache.BR_pref * (cache.Gamma_p + cache.Gamma_m);
}

double BDstarlnuDecay::A_FB() {
    double c_vv = 0.75 * std::pow(std::abs(cache.C_V1), 2) - std::pow(std::abs(cache.C_V2), 2);
    double c_aa = 1.5 * std::pow(std::abs(cache.C_A), 2);
    double c_tt = 12 * std::pow(std::abs(cache.C_T), 2);
    double c_ap = 1.5 * std::real(cache.C_A * std::conj(cache.C_P));
    double c_v1t = -6 * std::real(cache.C_V1 * std::conj(cache.C_T));
    double c_v2t = 6 * std::real(cache.C_V2 * std::conj(cache.C_T));
    double c_pt = -6 * std::real(cache.C_P * std::conj(cache.C_T));

    double g_vt_t0 = G_Vt_T0(cache.r_tau, cache.w_tau);
    double b_theta = c_vv * (F_Vp_1(cache.r_tau, cache.w_tau) - F_Vm_1(cache.r_tau, cache.w_tau)) 
            + c_aa * G_V0_Vt(cache.r_tau, cache.w_tau)
            + c_tt * (F_Tp_2(cache.r_tau, cache.w_tau) - F_Tm_2(cache.r_tau, cache.w_tau)) 
            + c_ap * G_V0_S(cache.r_tau, cache.w_tau)
            + c_v1t * (g_vt_t0 + G_Vp_Tp(cache.r_tau, cache.w_tau) + G_Vm_Tm(cache.r_tau, cache.w_tau))
            + c_v2t * (g_vt_t0 + G_Vm_Tp(cache.r_tau, cache.w_tau) + G_Vp_Tm(cache.r_tau, cache.w_tau))
            + c_pt * G_S_T0(cache.r_tau, cache.w_tau);

    if (cache.Gamma_p == 0)
        cache.Gamma_p = Gamma_tau_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = Gamma_tau_m(cache.r_tau, cache.w_tau);

    return b_theta / (cache.Gamma_p + cache.Gamma_m);
}

double BDstarlnuDecay::R_Dstar() {
    if (cache.Gamma_p == 0)
        cache.Gamma_p = Gamma_tau_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = Gamma_tau_m(cache.r_tau, cache.w_tau);

    double gamma_e = Gamma_tau_p(cache.r_e, cache.w_e) + Gamma_tau_m(cache.r_e, cache.w_e);

    return (cache.Gamma_p + cache.Gamma_m) / gamma_e;
}

double BDstarlnuDecay::P_tau() {
    if (cache.Gamma_p == 0)
        cache.Gamma_p = Gamma_tau_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = Gamma_tau_m(cache.r_tau, cache.w_tau);

    return (cache.Gamma_p - cache.Gamma_m) / (cache.Gamma_p + cache.Gamma_m);
}

double BDstarlnuDecay::P_D() {
    double c_aa = 0.5 * std::pow(std::abs(cache.C_A), 2);
    double c_pp = 1.5 * std::pow(std::abs(cache.C_P), 2);
    double c_tt = 8 * std::pow(std::abs(cache.C_T), 2);
    double c_ap = 3 * std::real(cache.C_A * std::conj(cache.C_P));
    double c_at = -12 * std::real(cache.C_A * std::conj(cache.C_T));

    double gamma_D_0 = c_aa * (2 * F_V0_1(cache.r_tau, cache.w_tau) + F_V0_2(cache.r_tau, cache.w_tau) + 3 * F_Vt(cache.r_tau, cache.w_tau)) 
            + c_pp * F_S(cache.r_tau, cache.w_tau)
            + c_tt * (F_T0_1(cache.r_tau, cache.w_tau) + 2 * F_T0_2(cache.r_tau, cache.w_tau)) 
            + c_ap * G_Vt_S(cache.r_tau, cache.w_tau)
            + c_at * G_V0_T0(cache.r_tau, cache.w_tau);

    if (cache.Gamma_p == 0)
        cache.Gamma_p = Gamma_tau_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = Gamma_tau_m(cache.r_tau, cache.w_tau);

    return gamma_D_0 / (cache.Gamma_p + cache.Gamma_m);
}


std::vector<ObservableValue> BDstarlnuDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_B__DSTAR_TAU_NU:   
        value = BR();
        break;
    case Observables::A_FB_B__DSTAR_TAU_NU:   
        value = A_FB();
        break;
    case Observables::R_DSTAR:   
        value = R_Dstar();
        break;
    case Observables::P_TAU_B__DSTAR_TAU_NU:   
        value = P_tau();
        break;
    case Observables::P_D_B__DSTAR_TAU_NU:   
        value = P_D();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BDstarlnuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
