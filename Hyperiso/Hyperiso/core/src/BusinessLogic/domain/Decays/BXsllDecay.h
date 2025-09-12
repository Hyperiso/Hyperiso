#ifndef __BXSLLDECAY_H__
#define __BXSLLDECAY_H__

#include "DecayParent.h"
#include "ObsQCDProxy.h"
#include "Include.h"
#include "Math.h"

/**
 * @brief Decay parent for the inclusive B > X_s l+ l- decays. Implements the integrated branching ratio in both q² \in [1, 6] GeV² and q² > 14.4 GeV², as well as the forward-backward asymmetry of the decay. 
 */
class BXsllDecay : public DecayParent {

private:
    static inline const std::array<double, 6> cc_res_mass      {3.096916, 3.68609,   3.77292, 4.039 , 4.153  , 4.421 };
    static inline const std::array<double, 6> cc_res_br        {5.93e-2 , 7.7e-3 ,   1.1e-5 , 1.4e-5, 1.0e-5 , 1.1e-5};
    static inline const std::array<double, 6> cc_res_width_tot {9.29e-5 , 3.04e-4,   2.73e-2, 8.0e-2, 1.03e-1, 6.2e-2};
    static inline const std::array<double, 6> cc_res_width_had {8.147e-5, 2.9746e-4, 2.36e-2, 5.2e-2, 7.8e-2 , 4.3e-2};

    static inline std::pair<double, double> q2_low_bound {1., 6.};
    static inline std::pair<double, double> q2_high_bound {14.4, 22.};

    static inline constexpr size_t LOOKUP_SIZE = 100;
    static inline constexpr size_t FF_ORDER = 20;
    static inline std::array<scalar_t, LOOKUP_SIZE> F_17_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> F_27_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> F_19_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> F_29_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> C7_new_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> C9_new_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> C10_new_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> CP7_new_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> CP9_new_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> CP10_new_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> C9_eff_LO_lookup;
    static inline std::array<scalar_t, LOOKUP_SIZE> CP9_eff_LO_lookup;

    static inline std::array<scalar_t, 50> delta_brems_lookup_low;
    static inline std::array<scalar_t, 50> delta_brems_lookup_high;

    static inline std::map<WCoef, complex_t> B_FR_wilson_cache;
    static inline std::map<WCoef, complex_t> BPrime_FR_wilson_cache;
    static inline std::map<WCoef, complex_t> BScalar_FR_wilson_cache;
    static inline std::map<WCoef, complex_t> B_FR_wilson_cache_LO;
    static inline std::map<WCoef, complex_t> BPrime_FR_wilson_cache_LO;

protected:
    double alpha_s(double mu);
    double f(double z);
    double h(double z);
    double g_rho(double z);
    double g_lambda(double z);
    double kappa(double f, double h, double alpha_s_mu_b);
    
    double m_hat(double m, double mb);
    double z(double mc_hat);

    double f_7(double s);
    double f_9(double s);
    
    scalar_t Gm1(double t);
    scalar_t G0(double t);
    scalar_t Delta_i_23(double s, double z, double w);
    scalar_t Delta_i_27(double s, double z, double w);
    double tau_22(double s, double w, scalar_t Delta_23, scalar_t Delta_27);
    scalar_t tau_27(double s, double w, scalar_t Delta_23, scalar_t Delta_27);
    scalar_t tau_28(double s, double w, scalar_t Delta_23, scalar_t Delta_27);
    scalar_t tau_29(double s, double w, scalar_t Delta_23, scalar_t Delta_27);
    double tau_77(double s);
    double tau_78(double s);
    double tau_88(double s);
    double tau_89(double s);
    double tau_99(double s);
    double tau_79(double s);
    scalar_t tau_210(double s, double z);
    double tau_710(double s);
    double tau_810(double s);
    double tau_910(double s);
    double sigma(double s);
    double sigma_9(double s);
    double sigma_7(double s, double L_mu);
    scalar_t F(double r);

    scalar_t Sigma_1(double s);
    double Sigma_2(double s);
    scalar_t Sigma_3(double s);
    scalar_t Sigma_7(double s, double z);
    double omega_22(double s, double L_l, double L_b_5_GeV);
    scalar_t omega_27(double s, double L_l, double L_b_5_GeV);
    scalar_t omega_29(double s, double L_l, double L_b_5_GeV);
    scalar_t omega_210(double s, double L_l, double L_b_5_GeV, double z);
    double omega_77(double s, double L_l);
    double omega_79(double s, double L_l);
    double omega_710(double s, double L_l);
    double omega_99(double s, double L_l);
    double omega_910(double s, double L_l);
    double omega_1010(double s, double L_l);

    scalar_t g(double z, double s);
    double breit_wigner(double s, double m_V, double br, double gamma_tot, double gamma_had, double m_b);
    double PV_breit_wigner(double s, double m_V, double br, double gamma_tot, double gamma_had, double m_b, double m_D_hat, double inv_alpha_em);
    double PV_R_cc_cont(double s, double m_D_hat);
    double R_cc(double s, double inv_alpha_em, double m_b);
    double R_cc_cont(double s);
    double PV_R_cc(double s, double inv_alpha_em, double m_b, double m_D_hat);
    scalar_t g_ld(double z, double s, double inv_alpha_em, double m_D_hat, double m_b);
    scalar_t C9_eff_base(double s, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b, QCDOrder order, bool prime);
    scalar_t C9_eff_LO(double s, bool prime);
    
    scalar_t F_17(double s);
    scalar_t F_27(double s);
    scalar_t F_19(double s);
    scalar_t F_29(double s);

    scalar_t C7_new_base(double s, double alpha_s_mu_b, double L_mu, bool prime);
    scalar_t C9_new_base(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b, bool prime);
    scalar_t C10_new_base(double s, double alpha_s_mu_b, bool prime);
    
    scalar_t C7_new(double s, bool prime);
    scalar_t C9_new(double s, bool prime);
    scalar_t C10_new(double s, bool prime);

    double W_7(double s, double alpha_s_mu_b, double L_mu);
    double W_9(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);
    double W_10(double s, double alpha_s_mu_b);
    scalar_t W_27(double s, double alpha_s_mu_b, double L_mu);
    scalar_t W_29(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);
    scalar_t W_210(double s, double alpha_s_mu_b);
    double W_79(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);
    double W_710(double s, double alpha_s_mu_b, double L_mu);
    double W_910(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);

    scalar_t H7(double s, double ml_hat, double alpha_s_mu_b);
    scalar_t H9(double s, double ml_hat, double alpha_s_mu_b);
    scalar_t H10(double s, double ml_hat, double alpha_s_mu_b);
    scalar_t H79(double s, double ml_hat, double alpha_s_mu_b);

    double dB0_ds(double s, double ml_hat, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);

    double delta_mb2(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);
    double delta_mb3(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);
    double delta_mc2(double s, double z, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);
    double delta_bremA(double s);
    double delta_bremB_base(double t, double z, double s_min, double s_max);
    double delta_bremB(double s, double m_b);
    double delta_em(double s, double L_l, double L_b_5_GeV);

    double pref_dB0_ds(double lambda_2, double g_lam_z, double f_z, double m_b, double rho_1, double g_rho_z);
    double pref_delta_mb2(double lambda_2, double m_b, scalar_t V_tb);
    double pref_delta_mb3(double rho_1, double m_b, scalar_t V_tb);
    double pref_delta_mc2(double lambda_2, double m_c, scalar_t V_tb, scalar_t V_ts, scalar_t V_cb, scalar_t V_cs);
    double pref_delta_brems(double alpha_s_mu_b);
    double pref_delta_em(double inv_alpha_em);
    double dB_ds(double s, double ml_hat, double alpha_s_mu_b, double z, double L_l, double L_b_5_GeV,
                    double pref_dB0_ds, double pref_delta_mb2, double pref_delta_mb3, double pref_delta_mc2, 
                    double pref_delta_brems, double pref_delta_em, double m_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat);

    double A_FB_0(double s, double ml_hat, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);

    double delta_A_mb2(double s, double alpha_s_mu_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat, double m_b);
    double delta_A_mc2(double s, double z, double alpha_s_mu_b);
    double delta_A_brem(double s, double z);
    double delta_A_em(double s, double L_l, double L_b_5_GeV, double z);

    double pref_A0_0(double lambda_2, double g_lam_z, double f_z, double m_b);
    double pref_A0_1(double lambda_1, double m_b);
    double A_FB(double s, double ml_hat, double alpha_s_mu_b, double z, double L_l, double L_b_5_GeV, double pref_A0_0, 
                double pref_A0_1, double pref_delta_mb2, double pref_delta_mc2, double pref_delta_brems, double pref_delta_em,
                double m_b, double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat);

    double ckm(scalar_t V_tb, scalar_t V_ts, scalar_t V_cb);
    double pref(double BR_BXclnu, double inv_alpha_em, double f, double kappa, double ckm);

    double BR_B_Xsll(double q2_min, double q2_max, double m_b, double pref, double ml_hat, double alpha_s_mu_b, 
                        double z, double L_l, double L_b_5_GeV, double pref_dB0_ds, double pref_delta_mb2,
                        double pref_delta_mb3, double pref_delta_mc2, double pref_delta_brems, double pref_delta_em,
                        double L_mu, double mc_hat, double inv_alpha_em, double m_D_hat);

public:
    BXsllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParent(matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar};
        this->max_order = QCDOrder::NNLO;
    }

    void build_op_tree() override;
};

#endif // __BXSLLDECAY_H__
