#include "BKstarGammaDecay.h"

using Charge = BKstarGammaConfig::B_Charge;

void BKstarGammaDecay::load_params() {
    cache.mu_b = w_config.hadronic_scale;
    cache.mu_h = sqrt(cache.mu_b * (*p)(ParamId{ParameterType::DECAY, "B_Ks", 14}, DataType::VALUE));
    fill_wilson_cache();
    
    cache.alpha_em = (*p)(ParamId{ParameterType::SM, "EW", {1, 2}}, DataType::VALUE);
    cache.m_b_m_b = (*p)(ParamId{ParameterType::SM, "SMINPUTS", 5}, DataType::VALUE);
    
    cache.lambda_hat_u = std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {0, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {0, 2}}, DataType::VALUE) 
                            / (std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {1, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {1, 2}}, DataType::VALUE));
    cache.m_b_mu_b = (*iobs_qcdp)(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    // cache.z = std::pow((*ports.iobs_qcdp)(MassConfig(4, cache.mu_b, MassType::MSBAR, MassType::POLE)) / cache.m_b_mu_b, 2);
    cache.z = 9.9551e-02; // ASK : Discrepancy between Isospin Asymmetry and BR. Need to homogeneize with QCDfCalculator
    cache.f_Ks_par = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {323, 1}}, DataType::VALUE);
    cache.f_B = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {521, 1}}, DataType::VALUE);
    cache.alpha_s_mu_b = (*iobs_qcdp)(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    cache.alpha_s_mu_h = (*iobs_qcdp)(AlphasConfig(cache.mu_h, MassType::POLE, MassType::POLE));
    cache.lambda_B = (*p)(ParamId{ParameterType::DECAY, "B_Ks", 13}, DataType::VALUE) / (1. - cache.alpha_s_mu_h * log(pow(cache.mu_h, 2)) * 1.8 / (3. * PI));
    // cache.T1_B_Ks = (*p)(ParamId{ParameterType::DECAY, "B_Ks", 16}); // ASK : Make it coherent with the FF choice ?
    cache.mu_0 = (*p)(ParamId{ParameterType::DECAY, "B_Ks", 17}, DataType::VALUE);
    if (fpeq(cache.mu_0, -1.)) cache.mu_0 = cache.mu_b;
    cache.L_b = log(cache.mu_b / (*p)(ParamId{ParameterType::SM, "QCD", {5, 3}}, DataType::VALUE));
    cache.C_F = (*iobs_qcdp).get_constants()->C_F;
    cache.Nc = (*iobs_qcdp).get_constants()->Nc;
    cache.n_f = 5.0;
    double eta = cache.alpha_s_mu_h / (*iobs_qcdp)(AlphasConfig(1.0, MassType::POLE, MassType::POLE));
    cache.f_Ks_perp = (*p)(ParamId{ParameterType::FLAVOR, "FCONST", {323, 2}}, DataType::VALUE) * pow(eta, cache.C_F / (*iobs_qcdp).get_constants()->beta[4][0]);
    
    load_cfg_dependent_params();
    
    // printf("mu_b = %.4e\n", cache.mu_b);
    // printf("mu_h = %.4e\n", cache.mu_h);
    // printf("alpha_em = %.4e\n", cache.alpha_em);
    // printf("m_b(m_b) = %.4e\n", cache.m_b_m_b);
    // printf("tau_B = %.4e\n", cache.tau_B);
    // printf("m_B = %.4e\n", cache.m_B);
    // printf("m_K* = %.4e\n", cache.m_Ks);
    // printf("N' = %.4e + %.4e i\n", cache.N_prime.real(), cache.N_prime.imag());
    // printf("m_b(mu_b) = %.4e\n", cache.m_b_mu_b);
    // printf("z = %.4e\n", cache.z);
    // printf("f_K*_par = %.4e\n", cache.f_Ks_par);
    // printf("f_K*_perp = %.4e\n", cache.f_Ks_perp);
    // printf("f_B = %.4e\n", cache.f_B);
    // printf("lambda_B = %.4e\n", cache.lambda_B);
    // printf("L_b = %.4e\n", cache.L_b);
    // printf("alpha_s(mu_b) = %.4e\n", cache.alpha_s_mu_b);
    // printf("alpha_s(mu_h) = %.4e\n", cache.alpha_s_mu_h);
    // printf("mu_0 = %.4e\n", cache.mu_0);

    // complex_t HVp = H_V(1, false);
    // complex_t HVm = H_V(-1, false);
    // complex_t HVpbar = H_V(1, true);
    // complex_t HVmbar = H_V(-1, true);
    // complex_t T_perp_p = cache.qcdf_calculator.T_perp_p(0.0, false);
    // complex_t T_perp_m = cache.qcdf_calculator.T_perp_m(0.0, false);
    // complex_t pref = 2. * I * cache.N_prime * cache.m_b_m_b * (std::pow(cache.m_B, 2) - std::pow(cache.m_Ks, 2)) / (cache.m_B * std::sqrt(4 * PI * cache.alpha_em));
    // printf("T_perp_p = %.4e + %.4e i\n", T_perp_p.real(), T_perp_p.imag());
    // printf("T_perp_m = %.4e + %.4e i\n", T_perp_m.real(), T_perp_m.imag());
    // printf("T_1(0) = %.4e\n", cache.ff_calculator.get(BV_FF::T1, 0.0));
    
    // printf("pref = %.4e + %.4e i\n", pref.real(), pref.imag());

    // printf("HV+ = %.4e + %.4e i\n", HVp.real(), HVp.imag());
    // printf("HV- = %.4e + %.4e i\n", HVm.real(), HVm.imag());
    // printf("HV+bar = %.4e + %.4e i\n", HVpbar.real(), HVpbar.imag());
    // printf("HV-bar = %.4e + %.4e i\n", HVmbar.real(), HVmbar.imag());
}

std::vector<ObservableValue> BKstarGammaDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_B0__KSTAR0_GAMMA:   
        set_cfg_flags(BKstarGammaConfig::B_Charge::B_0);
        value = BR();
        break;
    case Observables::IA_B0__KSTAR0_GAMMA:   
        set_cfg_flags(BKstarGammaConfig::B_Charge::B_0);
        value = delta_0();
        break;
    case Observables::BR_B__KSTAR_GAMMA:   
        set_cfg_flags(BKstarGammaConfig::B_Charge::B_PLUS);
        value = BR();
        break;
    case Observables::IA_B__KSTAR_GAMMA:  
        set_cfg_flags(BKstarGammaConfig::B_Charge::B_PLUS); 
        value = delta_0();
        break;
    default:
        LOG_ERROR("IndexError", "Observable", ObservableMapper::str(obs), "doesn't belong to the decay", DecayMapper::str(this->id));
    }

    return {ObservableValue(ObservableMapper::to_id(obs), value)};
}

std::vector<ObservableValue> BKstarGammaDecay::compute_observable(ObservableId obs) {
    return compute_observable(ObservableMapper::enum_of(obs).value());
}

void BKstarGammaDecay::fill_wilson_cache() {
    cache.C.clear();
    cache.C_trad.clear();
    this->w_proxy->set_basis(WilsonBasis::B_STANDARD);

    auto b_wilsons  = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);

    WCoef bp_cached[1] {WCoef::CP7};

    for (const auto& [id, val] : b_wilsons) {
        cache.C[id] = val;
    }
    for (auto id : bp_cached) {
        cache.C[id] = bp_wilsons.at(id);
    }

    this->w_proxy->set_basis(WilsonBasis::B_TRADITIONAL);
    auto b_wilsons_trad  = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons_trad = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);

    // TODO : Check values for trad basis ?
    cache.C_trad[WCoef::C1] = -8.3741e-2;
    cache.C_trad[WCoef::C2] =  1.0281;
    cache.C_trad[WCoef::C3] =  9.5354e-3;
    cache.C_trad[WCoef::C4] = -2.8921e-2;
    cache.C_trad[WCoef::C5] =  6.9522e-3;
    cache.C_trad[WCoef::C6] = -3.0791e-2;
    cache.C_trad[WCoef::C7] = -2.8010e-1;
    cache.C_trad[WCoef::C8] = -1.6378e-1;

    for (const auto& [id, val] : b_wilsons_trad) {
        cache.C_trad[id] = val;
    }

    cache.C_trad[WCoef::C7] += bp_wilsons_trad[WCoef::CP7];
    cache.C_trad[WCoef::C8] += bp_wilsons_trad[WCoef::CP8];

    ObsParameterMutator().set(ParamId{ParameterType::WILSON, "B_SCALE", 1}, cache.mu_h);
    cache.C2_h = w_proxy->getFR(WGroup::B, WCoef::C2, w_config.order);
    cache.C8_h = w_proxy->getFR(WGroup::B, WCoef::C8, w_config.order) + w_proxy->getFR(WGroup::BPrime, WCoef::CP8, w_config.order);

    // TODO : Check values for trad basis ?
    cache.C2_h =  9.8536e-1;
    cache.C8_h = -1.7263e-1;

    ObsParameterMutator().set(ParamId{ParameterType::WILSON, "B_SCALE", 1}, cache.mu_b);
    this->w_proxy->set_basis(WilsonBasis::B_STANDARD);
}

void BKstarGammaDecay::load_cfg_dependent_params() {
    cache.ff_calculator = BVFFCalculator(
        cfg.charge == Charge::B_0 ? 511 : 521,
        cfg.charge == Charge::B_0 ? 313 : 323,
        p,
        cfg.ff_src
    );

    cache.qcdf_calculator = BVQCDfCalculator(
        cfg.charge == Charge::B_0 ? 511 : 521,
        cfg.charge == Charge::B_0 ? 313 : 323,
        cache.mu_b,
        cache.C,
        std::make_shared<BVFFCalculator>(cache.ff_calculator),
        B_FF_Type::FULL,
        p,
        iobs_qcdp
    );

    cache.tau_B = (*p)(ParamId{ParameterType::FLAVOR, "FLIFE", cfg.charge == Charge::B_0 ? 511 : 521}, DataType::VALUE) / HBAR;
    cache.m_B = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 511 : 521}, DataType::VALUE);
    cache.m_Ks = (*p)(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 313 : 323}, DataType::VALUE);
    cache.N_prime = -(*p)(ParamId{ParameterType::SM, "SMINPUTS", 2}, DataType::VALUE) * cache.m_B * std::conj((*p)(ParamId{ParameterType::SM, "VCKM", {2, 1}}, DataType::VALUE)) * (*p)(ParamId{ParameterType::SM, "VCKM", {2, 2}}, DataType::VALUE) * cache.alpha_em / (PI * RT2);
}

void BKstarGammaDecay::set_cfg_flags(BKstarGammaConfig::B_Charge charge) {
    if (cfg.charge != charge) {
        cfg.charge = charge;
        load_cfg_dependent_params();
    }
}

complex_t BKstarGammaDecay::H_V(double sign, bool bar) {
    complex_t pref = 2. * I * cache.N_prime * cache.m_b_m_b * (std::pow(cache.m_B, 2) - std::pow(cache.m_Ks, 2)) / (cache.m_B * std::sqrt(4 * PI * cache.alpha_em));
    if (bar) pref = -std::conj(pref);

    double T1_0 = cache.ff_calculator.get(BV_FF::T1, 0.0);
    complex_t T_perp_p = cache.qcdf_calculator.T_perp_p(0.0, bar);
    complex_t T_perp_m = cache.qcdf_calculator.T_perp_m(0.0, bar);

    complex_t C7 = -sign * (sign == 1 ? cache.C[WCoef::CP7] : cache.C[WCoef::C7]); 

    return pref * (C7 * T1_0 + (T_perp_m - sign * T_perp_p) / 2.0);
}

complex_t BKstarGammaDecay::K1() {
    double F_p = cache.qcdf_calculator.F_perp(0.0);

    // printf("F_perp = %.4e\n", F_p);
    // printf("G_perp = %.4e + %.4e i\n", std::real(cache.qcdf_calculator.G_perp()), std::imag(cache.qcdf_calculator.G_perp()));
    // printf("X_perp = %.4e + %.4e i\n", std::real(cache.qcdf_calculator.X_perp(0.0)), std::imag(cache.qcdf_calculator.X_perp(0.0)));

    return -(cache.C_trad[WCoef::C6] + cache.C_trad[WCoef::C5] / cache.Nc) * F_p
           + cache.C_F * cache.alpha_s_mu_b / (4 * cache.Nc * PI) * (
                pow(cache.m_b_mu_b / cache.m_B, 2) * cache.C_trad[WCoef::C8] * cache.qcdf_calculator.X_perp(0.0)
              - cache.C_trad[WCoef::C2] * ((-4 * cache.L_b + 2) * F_p / 3 - cache.qcdf_calculator.G_perp()) 
              + F_p * log(cache.mu_b / cache.mu_0) * (
                    8. * cache.C_trad[WCoef::C3] / 3. 
                  + 4 * cache.n_f * (cache.C_trad[WCoef::C4] + cache.C_trad[WCoef::C6]) / 3. 
                  - 8. * (cache.Nc * cache.C_trad[WCoef::C6] + cache.C_trad[WCoef::C5]))
            );
}

complex_t BKstarGammaDecay::K2(int q) {
    // printf("H_perp = %.4e + %.4e i\n", std::real(cache.qcdf_calculator.H_perp()), std::imag(cache.qcdf_calculator.H_perp()));
    // printf("log(mu_b / mu_0) = %.4e\n", log(cache.mu_b / cache.mu_0));
    // printf("log(mu_b / mb_1S) = %.4e\n", cache.L_b);

    complex_t k2 = cache.C_trad[WCoef::C4] + cache.C_trad[WCoef::C3] / cache.Nc 
                    + cache.C_F * cache.alpha_s_mu_b / (4 * cache.Nc * PI) * (
                        cache.C_trad[WCoef::C2] * ((2 - 4 * cache.L_b) / 3. - cache.qcdf_calculator.H_perp()) 
                      + log(cache.mu_b / cache.mu_0) * (
                        -44. * cache.C_trad[WCoef::C3] / 3. 
                        -4. * cache.n_f * (cache.C_trad[WCoef::C4] + cache.C_trad[WCoef::C6]) / 3.)
                    );
    if (q == 2) 
        k2 += cache.lambda_hat_u * (cache.C_trad[WCoef::C2] + cache.C_trad[WCoef::C1] / cache.Nc);

    return k2;
}

double BKstarGammaDecay::BR() {
    double pref = (std::pow(cache.m_B, 2) - std::pow(cache.m_Ks, 2)) / 16. / PI / std::pow(cache.m_B, 3);
    double gamma = pref * (std::pow(std::abs(H_V(1, false)), 2) + std::pow(std::abs(H_V(-1, false)), 2));
    double gamma_bar = pref * (std::pow(std::abs(H_V(1, true)), 2) + std::pow(std::abs(H_V(-1, true)), 2));

    return cache.tau_B * (gamma + gamma_bar) / 2;
}

double BKstarGammaDecay::delta_0() {
    complex_t a7c = cache.C_trad[WCoef::C7] + cache.alpha_s_mu_b * cache.C_F * (cache.C_trad[WCoef::C2] * BV::G2(cache.z, cache.L_b) + cache.C_trad[WCoef::C8] * BV::G8(cache.L_b)) / (4. * PI) 
                    + cache.alpha_s_mu_h * cache.C_F * (cache.C8_h * cache.qcdf_calculator.H_8() + cache.C2_h * cache.qcdf_calculator.H_2()) / (4 * PI);

    // printf("H8a7 = %.4e\n", cache.qcdf_calculator.H_8());
	// printf("H2a7 = %.4e + %.4e i\n", std::real(cache.qcdf_calculator.H_2()), std::imag(cache.qcdf_calculator.H_2()));
	// printf("G8a7 = %.4e + %.4e i\n", std::real(BV::G8(cache.L_b)), std::imag(BV::G8(cache.L_b)));
	// printf("G2a7 = %.4e + %.4e i\n", std::real(BV::G2(cache.z, cache.L_b)), std::imag(BV::G2(cache.z, cache.L_b)));

    complex_t pref = 4. * PI2 * cache.f_B / (cache.m_b_mu_b * cache.ff_calculator.get(BV_FF::T1, 0.0) * a7c);
    complex_t t1 = cache.f_Ks_perp * K1() / cache.m_b_mu_b;
    complex_t f2 = cache.f_Ks_par * cache.m_Ks / (6. * cache.lambda_B * cache.m_B);
    complex_t bd = -pref * (t1 + f2 * K2(1));
    complex_t bu = 2. * pref * (t1 + f2 * K2(2));

    // printf("f_B = %.4e\n", cache.f_B);
    // printf("m_b(mu_b) = %.4e\n", cache.m_b_mu_b);
    // printf("T_1(0) = %.4e\n", cache.ff_calculator.get(BV_FF::T1, 0.0));
    // printf("f_K*_perp = %.4e\n", cache.f_Ks_perp);
    // printf("f_K*_par = %.4e\n", cache.f_Ks_par);
    // printf("m_K* = %.4e\n", cache.m_Ks);
    // printf("lambda_B = %.4e\n", cache.lambda_B);
    // printf("m_B = %.4e\n", cache.m_B);

    // printf("K1 = %.4e + %.4e i\n", std::real(K1()), std::imag(K1()));
    // printf("K2u = %.4e + %.4e i\n", std::real(K2(2)), std::imag(K2(2)));
    // printf("K2d = %.4e + %.4e i\n", std::real(K2(1)), std::imag(K2(1)));
    // printf("a7c = %.4e + %.4e i\n", std::real(a7c), std::imag(a7c));
    // printf("bu = %.4e + %.4e i\n", std::real(bu), std::imag(bu));
    // printf("bd = %.4e + %.4e i\n", std::real(bd), std::imag(bd));

    return std::real(bd - bu);
}
