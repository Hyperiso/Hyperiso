#include "BKstarGammaDecay.h"

using Charge = BKstarGammaConfig::B_Charge;

void BKstarGammaDecay::load_params() {
    ObsParameterProxy p;
    cache.mu_b = w_config.hadronic_scale;
    cache.mu_h = sqrt(cache.mu_b * p(ParamId{ParameterType::DECAY, "B_Ks", 14}));

    fill_wilson_cache();

    cache.ff_calculator = BVFFCalculator(
        cfg.charge == Charge::B_0 ? 511 : 521,
        cfg.charge == Charge::B_0 ? 313 : 323,
        cfg.ff_src
    );

    cache.qcdf_calculator = BVQCDfCalculator(
        cfg.charge == Charge::B_0 ? 511 : 521,
        cfg.charge == Charge::B_0 ? 313 : 323,
        w_config.hadronic_scale,
        cache.C,
        std::make_shared<BVFFCalculator>(cache.ff_calculator),
        B_FF_Type::FULL
    );
    
    cache.alpha_em = 1.0 / p(ParamId{ParameterType::SM, "SMINPUTS", 1});
    cache.m_b_m_b = p(ParamId{ParameterType::SM, "SMINPUTS", 5});
    cache.tau_B = p(ParamId{ParameterType::FLAVOR, "FLIFE", cfg.charge == Charge::B_0 ? 511 : 521}) / HBAR;
    cache.m_B = p(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 511 : 521});
    cache.m_Ks = p(ParamId{ParameterType::FLAVOR, "FMASS", cfg.charge == Charge::B_0 ? 313 : 323});
    cache.N_prime = -p(ParamId{ParameterType::SM, "SMINPUTS", 2}) * cache.m_B * std::conj(p(ParamId{ParameterType::SM, "VCKM", {2, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {2, 2}}) * cache.alpha_em / (PI * RT2);
    cache.lambda_hat_u = std::conj(p(ParamId{ParameterType::SM, "VCKM", {0, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {0, 2}}) 
                            / (std::conj(p(ParamId{ParameterType::SM, "VCKM", {1, 1}})) * p(ParamId{ParameterType::SM, "VCKM", {1, 2}}));
    cache.m_b_mu_b = ObsQCDProxy()(MassConfig(5, cache.mu_b, MassType::MSBAR, MassType::POLE));
    cache.z = std::pow(ObsQCDProxy()(MassConfig(4, cache.mu_b, MassType::MSBAR, MassType::POLE)) / cache.m_b_mu_b, 2);
    cache.f_Ks_par = p(ParamId{ParameterType::FLAVOR, "FCONST", {323, 1}});
    cache.f_B = p(ParamId{ParameterType::FLAVOR, "FCONST", {521, 1}});
    cache.lambda_B = p(ParamId{ParameterType::DECAY, "B_Ks", 13});
    // cache.T1_B_Ks = p(ParamId{ParameterType::DECAY, "B_Ks", 16});
    cache.mu_0 = p(ParamId{ParameterType::DECAY, "B_Ks", 17});
    if (fpeq(cache.mu_0, -1.)) cache.mu_0 = cache.mu_b;
    cache.L_b = log(cache.mu_b / p(ParamId{ParameterType::SM, "QCD", {5, 3}}));
    cache.C_F = ObsQCDProxy().get_constants()->C_F;
    cache.Nc = ObsQCDProxy().get_constants()->Nc;
    cache.n_f = 5.0; // TODO : Link with get_nf vs. hard-coded ?
    cache.alpha_s_mu_b = ObsQCDProxy()(AlphasConfig(cache.mu_b, MassType::POLE, MassType::POLE));
    cache.alpha_s_mu_h = ObsQCDProxy()(AlphasConfig(cache.mu_h, MassType::POLE, MassType::POLE));
    double eta = cache.alpha_s_mu_b / ObsQCDProxy()(AlphasConfig(1.0, MassType::POLE, MassType::POLE));
    cache.f_Ks_perp = p(ParamId{ParameterType::FLAVOR, "FCONST", {323, 2}}) * pow(eta, cache.C_F / ObsQCDProxy().get_constants()->beta[5][0]);
}

std::vector<ObservableValue> BKstarGammaDecay::compute_observable(Observables obs) {
    double value;
    switch (obs) {
    case Observables::BR_B__KSTAR_GAMMA:   
        value = BR();
        break;
    case Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA:   
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
    auto b_wilsons = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);
    WCoef bp_cached[1] {WCoef::CP7};

    for (auto p : b_wilsons) cache.C.emplace(p); 
    for (auto id : bp_cached) cache.C.emplace(std::pair{id, bp_wilsons.at(id)});

    this->w_proxy->set_basis(WilsonBasis::B_TRADITIONAL);
    auto b_wilsons_trad = w_proxy->getAFR(WGroup::B, this->w_config.order);
    auto bp_wilsons_trad = w_proxy->getAFR(WGroup::BPrime, this->w_config.order);

    cache.C_trad.emplace(WCoef::C1, b_wilsons_trad.at(WCoef::C1) /* + bp_wilsons_trad.at(WCoef::CP1) */ );
    cache.C_trad.emplace(WCoef::C2, b_wilsons_trad.at(WCoef::C2) /* + bp_wilsons_trad.at(WCoef::CP2) */ );
    cache.C_trad.emplace(WCoef::C3, b_wilsons_trad.at(WCoef::C3) /* + bp_wilsons_trad.at(WCoef::CP3) */ );
    cache.C_trad.emplace(WCoef::C4, b_wilsons_trad.at(WCoef::C4) /* + bp_wilsons_trad.at(WCoef::CP4) */ );
    cache.C_trad.emplace(WCoef::C5, b_wilsons_trad.at(WCoef::C5) /* + bp_wilsons_trad.at(WCoef::CP5) */ );
    cache.C_trad.emplace(WCoef::C6, b_wilsons_trad.at(WCoef::C6) /* + bp_wilsons_trad.at(WCoef::CP6) */ );
    cache.C_trad.emplace(WCoef::C7, b_wilsons_trad.at(WCoef::C7) + bp_wilsons_trad.at(WCoef::CP7));
    cache.C_trad.emplace(WCoef::C8, b_wilsons_trad.at(WCoef::C8) + bp_wilsons_trad.at(WCoef::CP8));

    ObsParameterMutator().set(ParamId{ParameterType::WILSON, "B_SCALE", 1}, cache.mu_h);
    cache.C2_h = w_proxy->getFR(WGroup::B, WCoef::C2, w_config.order) /* + w_proxy->getFR(WGroup::BPrime, WCoef::CP2, w_config.order) */;
    cache.C8_h = w_proxy->getFR(WGroup::B, WCoef::C8, w_config.order) + w_proxy->getFR(WGroup::BPrime, WCoef::CP8, w_config.order);
}

complex_t BKstarGammaDecay::H_V(double sign, bool bar) {
    complex_t pref = 2. * I * cache.N_prime * cache.m_b_m_b * (std::pow(cache.m_B, 2) - std::pow(cache.m_Ks, 2)) / (cache.m_B * std::sqrt(4 * PI * cache.alpha_em));
    if (bar) pref = std::conj(pref);

    double T1_0 = cache.ff_calculator.get(BV_FF::T1, 0.0);
    complex_t T_perp_p = cache.qcdf_calculator.T_perp_p(0.0, bar);
    complex_t T_perp_m = cache.qcdf_calculator.T_perp_m(0.0, bar);

    complex_t C7 = sign * (sign == 1 ? cache.C[WCoef::CP7] : cache.C[WCoef::C7]); 

    return pref * (C7 * T1_0 + (T_perp_m - sign * T_perp_p) / 2.0);
}

complex_t BKstarGammaDecay::K1() {
    double F_p = cache.qcdf_calculator.F_perp(0.0);
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
    complex_t k2 = cache.C_trad[WCoef::C4] + cache.C_trad[WCoef::C3] / cache.Nc 
                    + cache.C_F * cache.alpha_s_mu_b / (4 * cache.Nc * PI) * (
                        cache.C_trad[WCoef::C2] * ((2 + 4 * cache.L_b) / 3. - cache.qcdf_calculator.H_perp()) 
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

    return pref * cache.tau_B * (gamma + gamma_bar) / 2;
}

double BKstarGammaDecay::delta_0() {
    complex_t a7c = cache.C_trad[WCoef::C7] + cache.alpha_s_mu_b * cache.C_F * (cache.C_trad[WCoef::C2] * BV::G2(cache.z, cache.L_b) + cache.C_trad[WCoef::C8] * BV::G8(cache.L_b)) / (4. * PI) 
                    + cache.alpha_s_mu_h * cache.C_F * (cache.C8_h * cache.qcdf_calculator.H_8() - cache.C2_h * cache.qcdf_calculator.H_2()) / (4 * PI);

    complex_t pref = 4 * PI2 * cache.f_B / (cache.m_b_mu_b * cache.ff_calculator.get(BV_FF::T1, 0.0) * a7c);
    complex_t t1 = cache.f_Ks_perp * K1() / cache.m_b_mu_b;
    complex_t f2 = cache.f_Ks_par * cache.m_Ks / (6 * cache.lambda_B * cache.m_B);
    complex_t bd = -pref * (t1 + f2 * K2(1));
    complex_t bu = 2. * pref * (t1 + f2 * K2(2));
    return std::real(bd - bu);
}
