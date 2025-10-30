#include "BDlnuDecay.h"

double BDlnuDecay::t(double w) {
    return 1 + cache.r_D * (cache.r_D - 2 * w);
}

double BDlnuDecay::lambda_D(double w) {
    return 4 * cache.r_D * cache.r_D * (w * w - 1);
}

double BDlnuDecay::x_l(double rl, double w) {
    return rl * rl / t(w);
}

double BDlnuDecay::phi(double rl, double w) {
    return t(w) * std::sqrt(lambda_D(w)) * std::pow(1 - x_l(rl, w), 2);
}

double BDlnuDecay::w_max(double rD, double rl) {
    return (1 + rD * rD - rl * rl) / (2 * rD);
}

double BDlnuDecay::V_1(double w) {
    double z = (std::sqrt(1 + w) - RT2) / (std::sqrt(1 + w) + RT2);
    return 1 + z * (-8 * cache.rho_D2 + z * ((51 * cache.rho_D2 - 10) - z * (252 * cache.rho_D2 - 84)));
}

double BDlnuDecay::S_1(double w) {
    double u = w - 1;
    return V_1(w) * (1 + cache.Delta * (-0.019 + u * (0.041 - 0.015 * u)));
}

double BDlnuDecay::H_V0(double w) {
    return std::sqrt(cache.r_D * (w * w - 1) / t(w)) * (1 + cache.r_D) * V_1(w);
}

double BDlnuDecay::H_Vt(double w) {
    return std::sqrt(cache.r_D / t(w)) * (1 - cache.r_D) * (1 + w) * S_1(w);
}

double BDlnuDecay::H_S(double w) {
    return std::sqrt(cache.r_D) * (1 - cache.r_D) * (1 + w) * S_1(w) / cache.r_qm;
}

double BDlnuDecay::H_T(double w) {
    double a = std::sqrt(cache.r_D * (w * w - 1)) * cache.r_qp / (t(w) * (1 + cache.r_D));
    return -a * (std::pow(1 + cache.r_D, 2) * V_1(w) - 2 * cache.r_D * (1 + w) * S_1(w));
}

double BDlnuDecay::F_V0_1(double rl, double w_m) {
    if (!cache.C_V_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_V0(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_V0_2(double rl, double w_m) {
    if (!cache.C_V_flag) return 0;

    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_V0(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_Vt(double rl, double w_m) {
    if (!cache.C_V_flag) return 0;

    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_Vt(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_S(double rl, double w_m) {
    if (!cache.C_S_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_S(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_T_1(double rl, double w_m) {
    if (!cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::pow(H_T(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::F_T_2(double rl, double w_m) {
    if (!cache.C_T_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * std::pow(H_T(w), 2);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_V0_Vt(double rl, double w_m) {
    if (!cache.C_V_flag) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * x_l(rl, w) * H_V0(w) * H_Vt(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_V0_S(double rl, double w_m) {
    if (!(cache.C_V_flag && cache.C_S_flag)) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_V0(w) * H_S(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_V0_T(double rl, double w_m) {
    if (!(cache.C_V_flag && cache.C_T_flag)) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_V0(w) * H_T(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_Vt_S(double rl, double w_m) {
    if (!(cache.C_V_flag && cache.C_S_flag)) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_Vt(w) * H_S(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_Vt_T(double rl, double w_m) {
    if (!(cache.C_V_flag && cache.C_T_flag)) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * std::sqrt(x_l(rl, w)) * H_Vt(w) * H_T(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::G_S_T(double rl, double w_m) {
    if (!(cache.C_S_flag && cache.C_T_flag)) return 0;
    
    auto f = [this, rl] (double w) {
        return phi(rl, w) * H_S(w) * H_T(w);
    };

    return integrate(f, 1, w_m, 1e-3);
}

double BDlnuDecay::gamma_m(double r_l, double w_m) {
    double c_vv = std::pow(std::abs(cache.C_V), 2);
    double c_tt = 16 * std::pow(std::abs(cache.C_T), 2);
    double c_vt = -8 * std::real(cache.C_V * std::conj(cache.C_T));

    return c_vv * F_V0_1(r_l, w_m) + c_tt * F_T_2(r_l, w_m) + c_vt * G_V0_T(r_l, w_m);
}

double BDlnuDecay::gamma_p(double r_l, double w_m) {
    double c_vv = 0.5 * std::pow(std::abs(cache.C_V), 2);
    double c_ss = 1.5 * std::pow(std::abs(cache.C_S), 2);
    double c_tt = 8 * std::pow(std::abs(cache.C_T), 2);
    double c_vs = 3 * std::real(cache.C_V * std::conj(cache.C_S));
    double c_vt = -4 * std::real(cache.C_V * std::conj(cache.C_T));

    return c_vv * (F_V0_2(r_l, w_m) + 3 * F_Vt(r_l, w_m)) + c_ss * F_S(r_l, w_m) + c_tt * F_T_1(r_l, w_m) + c_vs * G_Vt_S(r_l, w_m) + c_vt * G_V0_T(r_l, w_m);
}

double BDlnuDecay::BR() {
    if (cache.Gamma_p == 0)
        cache.Gamma_p = gamma_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = gamma_m(cache.r_tau, cache.w_tau);

    return cache.BR_pref * (cache.Gamma_p + cache.Gamma_m);
}

double BDlnuDecay::A_FB() {
    double c_vv = 1.5 * std::pow(std::abs(cache.C_V), 2);
    double c_vs = 1.5 * std::real(cache.C_V * std::conj(cache.C_S));
    double c_vt = -6 * std::real(cache.C_V * std::conj(cache.C_T));
    double c_st = -6 * std::real(cache.C_S * std::conj(cache.C_T));
    double b_theta = c_vv * G_V0_Vt(cache.r_tau, cache.w_tau) + c_vs * G_V0_S(cache.r_tau, cache.w_tau) + c_vt * G_Vt_T(cache.r_tau, cache.w_tau) + c_st * G_S_T(cache.r_tau, cache.w_tau);

    if (cache.Gamma_p == 0)
        cache.Gamma_p = gamma_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = gamma_m(cache.r_tau, cache.w_tau);

    return b_theta / (cache.Gamma_p + cache.Gamma_m);
}

double BDlnuDecay::R_D() {
    double gamma_e = gamma_p(cache.r_e, cache.w_e) + gamma_m(cache.r_e, cache.w_e);

    if (cache.Gamma_p == 0)
        cache.Gamma_p = gamma_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = gamma_m(cache.r_tau, cache.w_tau);

    return (cache.Gamma_p + cache.Gamma_m) / gamma_e;
}

double BDlnuDecay::P_tau() {
    if (cache.Gamma_p == 0)
        cache.Gamma_p = gamma_p(cache.r_tau, cache.w_tau);

    if (cache.Gamma_m == 0)
        cache.Gamma_m = gamma_m(cache.r_tau, cache.w_tau);

    return (cache.Gamma_p - cache.Gamma_m) / (cache.Gamma_p + cache.Gamma_m);
}

void BDlnuDecay::load_params() {
    ObsParameterProxy p;
    double m_b = ObsQCDProxy()(MassConfig(5, w_config.hadronic_scale, MassType::POLE, MassType::POLE));
    double m_c = ObsQCDProxy()(MassConfig(4, w_config.hadronic_scale, MassType::POLE, MassType::POLE));
    double V_cb2 = std::pow(std::abs(p(ParamId{ParameterType::SM, "VCKM", {1, 2}})), 2);

    cache.G_F = p(ParamId{ParameterType::SM, "SMINPUTS", 2});
    cache.m_e = p(ParamId{ParameterType::SM, "MASS", 11});
    cache.m_tau = p(ParamId{ParameterType::SM, "MASS", 15});
    cache.m_B = p(ParamId{ParameterType::FLAVOR, "FMASS", 521});
    cache.m_D = p(ParamId{ParameterType::FLAVOR, "FMASS", 421});
    cache.tau_B = p(ParamId{ParameterType::FLAVOR, "FLIFE", 521});
    cache.V11 = p(ParamId{ParameterType::DECAY, "B_Dlnu", 1});
    cache.rho_D2 = p(ParamId{ParameterType::DECAY, "B_Dlnu", 2});
    cache.Delta= p(ParamId{ParameterType::DECAY, "B_Dlnu", 3});
    cache.r_D = cache.m_D / cache.m_B;
    cache.r_e = cache.m_e / cache.m_B;
    cache.r_tau = cache.m_tau / cache.m_B;
    cache.r_qp = (m_b + m_c) / cache.m_B;
    cache.r_qp = (m_b - m_c) / cache.m_B;
    cache.w_e = w_max(cache.r_D, cache.r_e);
    cache.w_tau = w_max(cache.r_D, cache.r_tau);
    cache.BR_pref = std::pow(cache.G_F * cache.m_B * cache.m_B * cache.V11, 2) * cache.m_D * cache.tau_B * V_cb2 / (96 * PI3 * HBAR);
    cache.C_V = w_proxy->getFM(WGroup::CC_bc, WCoef::C_V1_bc, QCDOrder::LO) + w_proxy->getFM(WGroup::CC_bc, WCoef::C_V2_bc, QCDOrder::LO);
    cache.C_S = w_proxy->getFM(WGroup::CC_bc, WCoef::C_S1_bc, QCDOrder::LO) + w_proxy->getFM(WGroup::CC_bc, WCoef::C_S2_bc, QCDOrder::LO);
    cache.C_T = w_proxy->getFM(WGroup::CC_bc, WCoef::C_T_bc, QCDOrder::LO);
    cache.C_V_flag = !fpeq(std::abs(cache.C_V), 0.0);
    cache.C_S_flag = !fpeq(std::abs(cache.C_S), 0.0);
    cache.C_T_flag = !fpeq(std::abs(cache.C_T), 0.0);
    cache.Gamma_p = 0.0;
    cache.Gamma_m = 0.0;
}

std::vector<ObservableValue> BDlnuDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_B__D_TAU_NU:   
        value = BR();
        break;
    case Observables::A_FB_B__D_TAU_NU:   
        value = A_FB();
        break;
    case Observables::R_D:   
        value = R_D();
        break;
    case Observables::P_TAU_B__D_TAU_NU:   
        value = P_tau();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BDlnuDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}
