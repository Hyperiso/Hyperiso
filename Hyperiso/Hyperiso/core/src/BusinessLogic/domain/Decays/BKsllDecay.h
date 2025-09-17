#ifndef __BKSTARLLDECAY_H__
#define __BKSTARLLDECAY_H__

#include "DecayParent.h"
#include "Include.h"
#include "ObsQCDProxy.h"

enum class FF {A0, A1, A12, V, T1, T2, T23, A2, T3};

struct BKstarllConfig {
    enum class FF_Src {BSZ_SR_LAT, BSZ_SR, GRvDV, GKvD_SR_LAT, GKvD_SR, HLMW};
    enum class Power_Corrections_Impl {SuperIso, Van_Dyk, Khodjamirian};
    enum class FF_Type {FULL, SOFT};
    enum class B_Charge {B_0, B_PLUS};
    enum class Lepton {E, MU, TAU};

    FF_Src ff_src {FF_Src::BSZ_SR_LAT};
    FF_Type ff_type {FF_Type::FULL};
    Power_Corrections_Impl power_corr_impl {Power_Corrections_Impl::SuperIso};
    B_Charge charge {B_Charge::B_PLUS};
    Lepton gen {Lepton::MU};
};

struct BKstarllCache {
    std::map<FF, std::array<double, 3>> alpha_ai;
    std::map<FF, double> m_R;
    std::map<WCoef, complex_t> C;
    std::map<WCoef, complex_t> C_bar;
    double C_F, Nc, beta_0;
    double mu_f;
    double tp, tm, t0;
    double z0;
    double f_K_perp;
    double lambda_B_p;
    double pref_T_perp, pref_T_par;
    double omega_0;
    double e_q, e_d {-1./3}, e_u {2./3};
    double m_b_PS, m_b_pole, m_b_mu_b;
    double m_c_pole, m_c_mu_b;
    complex_t T_par_m_0;
    double a_1_perp, a_2_perp;
    double a_1_par, a_2_par;
    double alpha_s_mu_b, alpha_s_mu_f, alpha_s_mb_pole, alpha_s_1_GeV;
    double eta_f;
    complex_t lambda_hat_u;
    double Delta_M;
    double z_c; // = (mc_pole / mb_PS)^2
    double L_b; // = log(mu_b / mb_PS)
    complex_t N_0;
    double m_l;
    double q2_min, q2_max;
    double q2_low, q2_high;
    double kappa;
};

/**
 * @brief Decay parent for the exclusive B > K* l+ l- decays. Implements the integrated branching ratio and angular observables in several q² bins, as well as the forward-backward asymmetry of the decay. 
 */
class BKstarllDecay : public ConfigurableDecayParent<BKstarllDecay, BKstarllConfig::FF_Src, BKstarllConfig::FF_Type, BKstarllConfig::Power_Corrections_Impl, BKstarllConfig::B_Charge> {
public:
    BKstarllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : ConfigurableDecayParent(DecayMapper::to_id(Decays::B__Kstar_l_l), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar};
        this->max_order = QCDOrder::NNLO;
    }

    void build_op_tree() override;

    void set_config_flag_c(BKstarllConfig::FF_Src src) { cfg.ff_src = src; };
    void set_config_flag_c(BKstarllConfig::FF_Type tp) { cfg.ff_type = tp; };
    void set_config_flag_c(BKstarllConfig::Power_Corrections_Impl impl) { cfg.power_corr_impl = impl; };
    void set_config_flag_c(BKstarllConfig::B_Charge charge) { cfg.charge = charge; };
    void set_config_flag_c(BKstarllConfig::Lepton lepton_gen) { cfg.gen = lepton_gen; };

private:
    BKstarllConfig cfg {};
    BKstarllCache cache;
    bool FF_loaded {false};

    static inline constexpr std::array<std::array<double, 6>, 6> P_bar {{
        { 0.5000, 0.0000, 0.0000,  0.0000, 0.0000,  0.0000},
        {-0.1667, 1.0000, 0.0000,  0.0000, 0.0000,  0.0000},
        { 0.0000, 0.0000, 1.0000, -0.1667, 16.000, -2.6667},
        { 0.0000, 0.0000, 0.0000,  0.5000, 0.0000,  8.0000},
        { 0.0000, 0.0000, 1.0000, -0.1667, 4.0000, -0.6667},
        { 0.0000, 0.0000, 0.0000,  0.5000, 0.0000,  2.0000}
    }};

    static inline constexpr size_t LOOKUP_SIZE = 50;
    static inline std::array<scalar_t, LOOKUP_SIZE> T_perp_p_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> T_perp_m_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> T_par_p_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> T_par_m_lookup;

protected:
    // Auxiliary
    double z(double t);
    double pole(double q2, double m_R);
    complex_t h(double q2, double m_q, double mu_b);
    double phi_Kstar(double u, double a1, double a2);
    complex_t B_0(double s, double m_q);
    complex_t L_1(complex_t x);
    complex_t I_1(double u, double m_q, double q2, double m_Bd);
    complex_t F_27_u(double s_hat, double l);
    complex_t F_19_u(double s_hat, double l);
    complex_t F_29_u(double s_hat, double l);
    complex_t A_Seidel(double s, double m_b, double mu_b);
    complex_t B_Seidel(double s, double m_b, double mu_b);
    complex_t C_Seidel(double s, double mu_b);

    void fill_wilson_cache();
    void fill_wilson_bar_cache();
    std::shared_ptr<ParameterNode> get_lepton_mass();

    // Form Factors
    std::shared_ptr<OperatorNode> load_FF_params();
    double F_a(FF a, double q2);
    double E_K(double q2, double m_B, double m_K);
    double A_2(double q2, double m_B, double m_K);
    double T_3(double q2, double m_B, double m_K);
    double xi_perp(double q2, double m_B, double m_K);
    double xi_par(double q2, double m_B, double m_K);
    double f_perp(double q2, double m_B, double m_K);
    double f_par(double q2, double m_B, double m_K);
    double f_0(double q2, double m_B, double m_K);

    // QCDf 
    double F_perp(double s_hat);
    double X_perp(double s_hat);
    double gv_dga_4(double u, double z3a, double z3v, double w10a, double dtp, double dtm);
    complex_t F_V(double v);

    complex_t Y(double q2);
    complex_t Y_u(double q2);
    double L(double q2);
    complex_t C_perp_0(double q2, double m_B, double sign);
    complex_t C_par_0(double q2, double m_B, double sign);
    complex_t C_perp_f(double q2, double sign);
    complex_t C_par_f(double q2, double sign);
    complex_t C_perp_nf(double q2, double m_B);
    complex_t C_par_nf(double q2, double m_B);

    complex_t t_perp(double u, double m_q, double q2, double E_Kstar, double m_B);
    complex_t t_par(double u, double m_q, double q2, double E_Kstar, double m_B);
    complex_t T_par_m_0(double m_B);
    complex_t T_par_p_p_f(double u, double q2, double m_B, double m_K);
    complex_t T_par_p_m_f(double u, double q2, double m_B, double m_K);
    complex_t T_perp_p_p_f(double u, double q2, double m_B, double m_K);
    complex_t T_perp_p_m_f(double u, double q2, double m_B, double m_K);
    complex_t T_perp_p_nf(double u, double q2, double m_B, double m_K);
    complex_t T_par_p_nf(double u, double q2, double m_B, double m_K);
    complex_t T_par_m_nf(double u, double q2, double m_B, double m_K);

    complex_t inv_lambda_B_m(double q2, double m_B);

    complex_t I_perp_p(double q2, double m_B, double m_K); 
    complex_t I_perp_m(double q2, double m_B, double m_K);
    complex_t I_par_p(double q2, double m_B, double m_K);
    complex_t I_par_m(double q2, double m_B, double m_K);

    complex_t I_HSA_1(double q2, double m_B);
    complex_t I_HSA_2(double q2, double m_B, double z3a, double z3v, double w10a, double dtp, double dtm);
    complex_t delta_T_perp_WA(double q2, double m_B, double m_K, double f_B, double f_K_par);
    complex_t delta_T_perp_HSA(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm);

    complex_t T_perp_p(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm);
    complex_t T_perp_m(double q2, double m_B, double m_K, double f_B, double f_K_par, double z3a, double z3v, double w10a, double dtp, double dtm);
    complex_t T_par_p(double q2, double m_B, double m_K);
    complex_t T_par_m(double q2, double m_B, double m_K);
    complex_t Delta_par(double q2, double m_B, double m_K, double f_B, double f_K_par);

    complex_t T_perp_p_cached(double q2);
    complex_t T_perp_m_cached(double q2);
    complex_t T_par_p_cached(double q2);
    complex_t T_par_m_cached(double q2);

    // Transversity amplitudes
    double beta_l(double q2, double m_l);
    double lambda(double q2, double m_B, double m_K);
    complex_t N(double q2, double m_B, double m_K, double m_l);

    complex_t A_perp_low(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_par_low(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_0_low(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_t_low(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s);
    complex_t A_S_low(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s);

    complex_t C7_eff(double q2, double m_B);
    complex_t C9_eff(double q2, double m_B);
    complex_t A_perp_high(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_par_high(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_0_high(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_t_high(double q2, double m_B, double m_K, double m_l, double m_s);
    complex_t A_S_high(double q2, double m_B, double m_K, double m_l, double m_s);

    complex_t interpolate(double q2, complex_t val_low, complex_t val_high);
    complex_t A_perp(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_par(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_0(double q2, double m_B, double m_K, double m_l, double sign);
    complex_t A_t(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s);
    complex_t A_S(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s);
    
    // Angular coefficients
    double J1s(double q2, double m_B, double m_K, double m_l);
    double J1c(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s);
    double J2s(double q2, double m_B, double m_K, double m_l);
    double J2c(double q2, double m_B, double m_K, double m_l);
    double J3(double q2, double m_B, double m_K, double m_l);
    double J4(double q2, double m_B, double m_K, double m_l);
    double J5(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s);
    double J6s(double q2, double m_B, double m_K, double m_l);
    double J6c(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s);
    double J7(double q2, double m_B, double m_K, double m_l, double f_B, double f_K_par, double m_s);
    double J8(double q2, double m_B, double m_K, double m_l);
    double J9(double q2, double m_B, double m_K, double m_l);

    // Observables

};

#endif // __BXSLLDECAY_H__
