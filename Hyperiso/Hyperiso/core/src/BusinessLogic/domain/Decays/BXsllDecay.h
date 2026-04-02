#ifndef __BXSLLDECAY_H__
#define __BXSLLDECAY_H__

#include "DecayParent.h"
#include "ObsQCDProxy.h"
#include "Include.h"
#include "Math.h"
#include "DefaultConfig.h"

struct BXsllConfig : public DecayConfig {
    enum class Lepton {E, MU, TAU};
    Lepton gen {Lepton::MU};
};

struct BXsllDecayCache {
    double alpha_em;
    double m_b_1S;
    double m_l_hat, m_c_hat, m_D_hat;
    double z;
    double L_b, L_b_5GeV, L_l;
    double alpha_s_mu_b;
    double pref_dB0_ds, pref_dB_ds;
    double pref_A0_0, pref_A0_1;
    double pref_delta_mb2, pref_delta_mb3, pref_delta_mc2, pref_delta_brems, pref_delta_em;

    std::array<double, 3> rand_err;
    
    std::map<WCoef, complex_t> C;
    std::map<WCoef, complex_t> C_LO;

    std::array<complex_t, 100> F_17_lookup;
    std::array<complex_t, 100> F_27_lookup;
    std::array<complex_t, 100> F_19_lookup;
    std::array<complex_t, 100> F_29_lookup;
    std::array<double, 100> delta_brems_lookup;

    std::array<double, 6> cc_res_mass;
    std::array<double, 6> cc_res_br;
    std::array<double, 6> cc_res_width_tot;
    std::array<double, 6> cc_res_width_had;
};

/**
 * @brief Decay parent for the inclusive B > X_s l+ l- decays. Implements the integrated branching ratio in both q² \in [1, 6] GeV² and q² > 14.4 GeV², as well as the forward-backward asymmetry of the decay. 
 */
class BXsllDecay : public DecayParentConfigurable<BXsllConfig> {

private:
    BXsllDecayCache cache;
    BXsllConfig cfg;

protected:
    // utility
    void set_cfg_flags(BXsllConfig::Lepton gen);
    void fill_wilson_cache();
    void load_cfg_dep_params();

    double f(double z);
    double h(double z);
    double g_rho(double z);
    double g_lambda(double z);
    double kappa(double z);

    double f_7(double s);
    double f_9(double s);
    
    complex_t Gm1(double t);
    complex_t G0(double t);
    complex_t Delta_i_23(double s, double z, double w);
    complex_t Delta_i_27(double s, double z, double w);
    double tau_22(double s, double w, complex_t Delta_23, complex_t Delta_27);
    complex_t tau_27(double s, double w, complex_t Delta_23, complex_t Delta_27);
    complex_t tau_28(double s, double w, complex_t Delta_23, complex_t Delta_27);
    complex_t tau_29(double s, double w, complex_t Delta_23, complex_t Delta_27);
    double tau_77(double s);
    double tau_78(double s);
    double tau_88(double s);
    double tau_89(double s);
    double tau_99(double s);
    double tau_79(double s);
    complex_t tau_210(double s, double z);
    double tau_710(double s);
    double tau_810(double s);
    double tau_910(double s);
    double sigma(double s);
    double sigma_9(double s);
    double sigma_7(double s, double L_mu);
    complex_t F(double r);

    complex_t Sigma_1(double s);
    double Sigma_2(double s);
    complex_t Sigma_3(double s);
    complex_t Sigma_7(double s, double z);
    double omega_22(double s, double L_l);
    complex_t omega_27(double s, double L_l);
    complex_t omega_29(double s, double L_l);
    complex_t omega_210(double s, double L_l, double z);
    double omega_77(double s, double L_l);
    double omega_79(double s, double L_l);
    double omega_710(double s, double L_l);
    double omega_99(double s, double L_l);
    double omega_910(double s, double L_l);
    double omega_1010(double s, double L_l);

    complex_t g(double z, double s);
    double breit_wigner(double s, double m_V, double br, double gamma_tot, double gamma_had);
    double PV_breit_wigner(double s, double m_V, double br, double gamma_tot, double gamma_had);
    double PV_R_cc_cont(double s);
    double R_cc(double s);
    double R_cc_cont(double s);
    double PV_R_cc(double s);
    complex_t g_ld(double z, double s);
    complex_t C9_eff(double s, QCDOrder order, bool prime);
    
    complex_t F_17(double s);
    complex_t F_27(double s);
    complex_t F_19(double s);
    complex_t F_29(double s);

    complex_t C7_new(double s, bool prime);
    complex_t C9_new(double s, bool prime);
    complex_t C10_new(double s, bool prime);

    double W_7(double s);
    double W_9(double s);
    double W_10(double s);
    complex_t W_27(double s);
    complex_t W_29(double s);
    complex_t W_210(double s);
    double W_79(double s);
    double W_710(double s);
    double W_910(double s);

    double dB0_ds(double s, double ml_hat);
    double delta_mb2(double s);
    double delta_mb3(double s);
    double delta_mc2(double s);
    double delta_bremA(double s);
    double delta_bremB_base(double s);
    double delta_bremB(double s);
    double delta_em(double s, double L_l);
    double dB_ds(double s, double ml_hat, double L_l);

    double A_FB_0(double s, double ml_hat);
    double delta_A_mb2(double s);
    double delta_A_mc2(double s);
    double delta_A_brem(double s);
    double delta_A_em(double s, double L_l);
    double A_FB(double s, double ml_hat, double L_l);

    std::vector<ObservableValue> BR_B_Xsll(Observables oid);
    std::vector<ObservableValue> A_FB_B_Xsll(Observables oid);

public:
    BXsllDecay(QCDOrder order, double matching_scale, double hadronic_scale, ObservablePortsConfig& ports) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__Xs_l_l), matching_scale, hadronic_scale, order, ports) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime), GroupMapper::to_id(WGroup::BScalar)};
        this->max_order = QCDOrder::NNLO;
        this->cache = BXsllDecayCache {};
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
    void set_config_spe(BXsllConfig config) override {this->cfg = config;}
};

#endif // __BXSLLDECAY_H__
