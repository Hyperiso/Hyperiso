#ifndef __BKSTARLLDECAY_H__
#define __BKSTARLLDECAY_H__

#include "DecayParent.h"
#include "Include.h"
#include "ObsQCDProxy.h"
#include "DefaultConfig.h"

enum class FF {A0, A1, A12, V, T1, T2, T23, A2, T3};

struct BKstarllConfig : public DecayConfig {
    enum class FF_Src {BSZ_SR_LAT, BSZ_SR, GRvDV, GKvD_SR_LAT, GKvD_SR, HLMW};
    enum class Power_Corrections_Impl {BFS, BCvDV, KMPW};
    enum class FF_Type {FULL, SOFT};
    enum class B_Charge {B_0, B_PLUS};
    enum class Lepton {E, MU, TAU};

    FF_Src ff_src {FF_Src::BSZ_SR_LAT};
    FF_Type ff_type {FF_Type::FULL};
    Power_Corrections_Impl power_corr_impl {Power_Corrections_Impl::BFS};
    B_Charge charge {B_Charge::B_PLUS};
    Lepton gen {Lepton::MU};

    std::vector<std::pair<double, double>> bins {
        { 1.00,  6.00},
        {14.18, 19.00}
    };
};

struct BKstarllCache {
    std::map<FF, std::array<double, 3>> alpha_ai;
    std::map<FF, double> m_R;
    std::map<WCoef, complex_t> C;
    std::map<WCoef, complex_t> C_bar;

    double G_F, alpha_em;
    double m_l, m_s;
    double mu_b, mu_f;
    double m_c_pole, m_c_mu_b;
    double m_b_PS, m_b_pole, m_b_mu_b;
    double alpha_s_mu_b, alpha_s_mu_f, alpha_s_mb_pole, alpha_s_1_GeV;
    double m_B, m_Ks;
    double zeta_3_A, zeta_3_V, omega_10_A, delta_t_p, delta_t_m;
    double C_F, Nc, beta_0;
    double eta_f;
    double tp, tm, t0;
    double z0;
    double f_B, f_Ks_perp, f_Ks_par;
    double lambda_B_p;
    complex_t lambda_hat_u;
    double a_1_perp, a_2_perp;
    double a_1_par, a_2_par;
    double omega_0;
    double e_q, e_d {-1./3}, e_u {2./3};
    double z_c; // = (mc_pole / mb_PS)^2
    double L_b; // = log(mu_b / mb_PS)
    double Delta_M;
    double kappa;
    double pref_T_perp, pref_T_par;
    complex_t N_0;
    complex_t T_par_m_0;
    double q2_min, q2_max;
    double q2_low, q2_high;

    // Power corrections, hadronic fit
    std::array<complex_t, 3> h_p_fit;
    std::array<complex_t, 3> h_m_fit; 
    std::array<complex_t, 3> h_0_fit;

    // Power corrections, guesstimate
    std::array<double, 6> A_had_err_low_0;
    std::array<double, 6> A_had_err_low_1;
    std::array<double, 8> A_had_err_high;

    // NF Corrections, Van Dyk
    double tp_nf, t0_nf;
    double z0_nf, z_Jpsi_nf, z_psi2S_nf;
    std::array<complex_t, 3> alpha_perp;
    std::array<complex_t, 3> alpha_par;
    std::array<complex_t, 2> alpha_0;

    // NF Corrections, KMPW
    double q2_bar, q2_Jpsi; // q2bar = 1.0 GeV², q2_Jpsi = m_Jpsi^2
    std::array<double, 3> DeltaC9_M_qbar;
    std::array<double, 3> r1_M;
    std::array<double, 3> r2_M;

    // Lookups
    static inline constexpr size_t LOOKUP_SIZE = 50;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_p_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_m_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_par_p_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_par_m_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_p_bar_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_m_bar_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_par_p_bar_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_par_m_bar_lookup;

    std::array<std::vector<double>, 14> J_i_binned;
    std::array<std::vector<double>, 14> J_i_bar_binned; 
};

/**
 * @brief Decay parent for the exclusive B > K* l+ l- decays. Implements the integrated branching ratio and angular observables in several q² bins, as well as the forward-backward asymmetry of the decay. 
 */
class BKstarllDecay : public DecayParentConfigurable<BKstarllConfig> {
public:
    BKstarllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__Kstar_l_l), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar};
        this->max_order = QCDOrder::NNLO;
        // this->load_params();
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

    void set_config_spe(BKstarllConfig config) override {this->cfg = config;}

private:
    BKstarllConfig cfg {};
    BKstarllCache cache;

    static inline constexpr std::array<std::array<double, 6>, 6> P_bar {{
        { 0.5000, 0.0000, 0.0000,  0.0000, 0.0000,  0.0000},
        {-0.1667, 1.0000, 0.0000,  0.0000, 0.0000,  0.0000},
        { 0.0000, 0.0000, 1.0000, -0.1667, 16.000, -2.6667},
        { 0.0000, 0.0000, 0.0000,  0.5000, 0.0000,  8.0000},
        { 0.0000, 0.0000, 1.0000, -0.1667, 4.0000, -0.6667},
        { 0.0000, 0.0000, 0.0000,  0.5000, 0.0000,  2.0000}
    }};

protected:
    // Auxiliary
    complex_t z(double t, double t_p, double t_m);
    double pole(double q2, double m_R);
    complex_t h(double q2, double m_q);
    double phi_Kstar(double u, double a1, double a2);
    complex_t B_0(double s, double m_q);
    complex_t L_1(complex_t x);
    complex_t I_1(double u, double m_q, double q2);
    complex_t F_27_u(double s_hat);
    complex_t F_19_u(double s_hat);
    complex_t F_29_u(double s_hat);
    complex_t A_Seidel(double s);
    complex_t B_Seidel(double s);
    complex_t C_Seidel(double s);

    void fill_wilson_cache();
    void fill_wilson_bar_cache();
    void load_FF_params();

    // Form Factors
    double F_a(FF a, double q2);
    double E_K(double q2);
    double A_2(double q2);
    double T_3(double q2);
    double xi_perp(double q2);
    double xi_par(double q2);
    double f_perp(double q2);
    double f_par(double q2);
    double f_0(double q2);

    // QCDf 
    double F_perp(double s_hat);
    double X_perp(double s_hat);
    double gv_dga_4(double u);
    complex_t F_V(double v, bool bar);

    complex_t Y(double q2);
    complex_t Y_u(double q2);
    double L(double q2);
    complex_t C_perp_0(double q2, double sign, bool bar);
    complex_t C_par_0(double q2, double sign, bool bar);
    complex_t C_perp_f(double q2, double sign, bool bar);
    complex_t C_par_f(double q2, double sign, bool bar);
    complex_t C_perp_nf(double q2, bool bar);
    complex_t C_par_nf(double q2, bool bar);

    complex_t t_perp(double u, double m_q, double q2, double E_Kstar);
    complex_t t_par(double u, double m_q, double q2, double E_Kstar);
    complex_t T_par_p_p_f(double u, double q2, bool bar);
    complex_t T_par_p_m_f(double u, double q2, bool bar);
    complex_t T_perp_p_p_f(double u, double q2, bool bar);
    complex_t T_perp_p_m_f(double u, double q2, bool bar);
    complex_t T_perp_p_nf(double u, double q2, bool bar);
    complex_t T_par_p_nf(double u, double q2, bool bar);
    complex_t T_par_m_nf(double u, double q2, bool bar);

    complex_t inv_lambda_B_m(double q2);

    complex_t I_perp_p(double q2, bool bar); 
    complex_t I_perp_m(double q2, bool bar);
    complex_t I_par_p(double q2, bool bar);
    complex_t I_par_m(double q2, bool bar);

    complex_t I_HSA_1(double q2, bool bar);
    complex_t I_HSA_2(double q2, bool bar);
    complex_t delta_T_perp_WA(double q2);
    complex_t delta_T_perp_HSA(double q2, bool bar);

    complex_t T_perp_p(double q2, bool bar);
    complex_t T_perp_m(double q2, bool bar);
    complex_t T_par_p(double q2, bool bar);
    complex_t T_par_m(double q2, bool bar);
    complex_t Delta_par(double q2);

    complex_t T_perp_p_cached(double q2, bool bar);
    complex_t T_perp_m_cached(double q2, bool bar);
    complex_t T_par_p_cached(double q2, bool bar);
    complex_t T_par_m_cached(double q2, bool bar);

    // Nonfactorisable corrections
    double beta_l(double q2);
    double lambda(double q2);
    complex_t N(double q2, bool bar);

    complex_t delta_A_perp_QCDf(double q2, double sign, bool bar);
    complex_t delta_A_perp_vD(double q2, bool bar);
    complex_t delta_A_perp_K(double q2, bool bar);
    complex_t delta_A_perp(double q2, double sign, bool bar);

    complex_t delta_A_par_QCDf(double q2, double sign, bool bar);
    complex_t delta_A_par_vD(double q2, bool bar);
    complex_t delta_A_par_K(double q2, bool bar);
    complex_t delta_A_par(double q2, double sign, bool bar);

    complex_t delta_A_0_QCDf(double q2, double sign, bool bar);
    complex_t delta_A_0_vD(double q2, bool bar);
    complex_t delta_A_0_K(double q2, bool bar);
    complex_t delta_A_0(double q2, double sign, bool bar);

    // Transversity amplitudes
    complex_t A_perp_low(double q2, double sign, bool bar);
    complex_t A_par_low(double q2, double sign, bool bar);
    complex_t A_0_low(double q2, double sign, bool bar);
    complex_t A_t_low(double q2, bool bar);
    complex_t A_S_low(double q2, bool bar);

    complex_t C7_eff(double q2, bool bar);
    complex_t C9_eff(double q2, bool bar);
    complex_t A_perp_high(double q2, double sign, bool bar);
    complex_t A_par_high(double q2, double sign, bool bar);
    complex_t A_0_high(double q2, double sign, bool bar);
    complex_t A_t_high(double q2, bool bar);
    complex_t A_S_high(double q2, bool bar);

    complex_t interpolate(double q2, complex_t val_low, complex_t val_high);
    complex_t A_perp(double q2, double sign, bool bar);
    complex_t A_par(double q2, double sign, bool bar);
    complex_t A_0(double q2, double sign, bool bar);
    complex_t A_t(double q2, bool bar);
    complex_t A_S(double q2, bool bar);
    
    // Angular coefficients
    double J1s(double q2, bool bar);
    double J1c(double q2, bool bar);
    double J2s(double q2, bool bar);
    double J2c(double q2, bool bar);
    double J3(double q2, bool bar);
    double J4(double q2, bool bar);
    double J5(double q2, bool bar);
    double J6s(double q2, bool bar);
    double J6c(double q2, bool bar);
    double J7(double q2, bool bar);
    double J8(double q2, bool bar);
    double J9(double q2, bool bar);

    void compute_binned_J_i();

    std::vector<ObservableValue> dG_dq2_binned(bool bar);
    double dG_dq2_avg_bin(size_t bin);
    std::vector<ObservableValue> A_FB_binned();
    ObservableValue q0();
    std::vector<ObservableValue> A_CP_binned();
    std::vector<ObservableValue> F_L_binned();
    std::vector<ObservableValue> F_T_binned();
    std::vector<ObservableValue> A_T_1_binned();
    std::vector<ObservableValue> A_T_2_binned();
    std::vector<ObservableValue> A_T_3_binned();
    std::vector<ObservableValue> A_T_4_binned();
    std::vector<ObservableValue> A_T_5_binned();
    std::vector<ObservableValue> A_T_Re_binned();
    std::vector<ObservableValue> A_T_Re_CPV_binned();
    std::vector<ObservableValue> A_Im_binned();
    std::vector<ObservableValue> alpha_K_binned();
    std::vector<ObservableValue> H_T_1_binned();
    std::vector<ObservableValue> H_T_2_binned();
    std::vector<ObservableValue> H_T_3_binned();
    std::vector<ObservableValue> P_2_binned();
    std::vector<ObservableValue> P_3_binned();
    std::vector<ObservableValue> P_6_binned();
    std::vector<ObservableValue> P_8_binned();
    std::vector<ObservableValue> Pp_i_binned(size_t i, bool cpv=false);
    std::vector<ObservableValue> S_i_binned(size_t i, bool cpv=false);
    std::vector<ObservableValue> P_i_CPV_binned(size_t i);

    // Tests
    void test_ff();
    void test_T();
    void test_J();
    void test_binned_obs();
};

#endif // __BXSLLDECAY_H__
