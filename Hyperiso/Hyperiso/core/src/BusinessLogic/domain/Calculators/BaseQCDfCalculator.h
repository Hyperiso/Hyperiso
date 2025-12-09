#ifndef __BASEQCDFCALCULATOR_H__
#define __BASEQCDFCALCULATOR_H__

#include "Include.h"
#include "ObsQCDProxy.h"
#include "ObsParameterProxy.h"

enum class B_FF_Type {FULL, SOFT};

class BaseQCDfCalculator {
protected:
    double mu_b, L_b;
    double m_c_pole, m_b_pole, m_b_PS;
    double m_B, m_X;
    double f_B, f_X_par, f_X_perp {0};
    double z_c;
    double Delta_M;
    double e_q;
    double Lambda_h;
    double alpha_s_mu_f, alpha_s_mu_b;
    double loop_f_mu_f, loop_f_mu_b;
    double a_1_perp {0}, a_2_perp {0}, a_1_par, a_2_par;
    double zeta_3_A {0}, zeta_3_V {0}, omega_10_A {0}, delta_t_p {0}, delta_t_m {0};
    double lambda_B_p;
    double pref_perp {0}, pref_par;
    complex_t n_T_par_m_0, n_T_par_m_0_bar;
    complex_t lambda_hat_u;
    std::map<WCoef, complex_t> C;
    std::map<WCoef, complex_t> C_bar;

    double delta_qu;
    B_FF_Type ff_tp;
    std::string src_block;

    static inline constexpr double e_u = 2. / 3;
    static inline constexpr double e_d = -1. / 3;
    static inline constexpr std::array<std::array<double, 6>, 6> P_bar {{
        { 0.5000, 0.0000, 0.0000,  0.0000, 0.0000,  0.0000},
        {-0.1667, 1.0000, 0.0000,  0.0000, 0.0000,  0.0000},
        { 0.0000, 0.0000, 1.0000, -0.1667, 16.000, -2.6667},
        { 0.0000, 0.0000, 0.0000,  0.5000, 0.0000,  8.0000},
        { 0.0000, 0.0000, 1.0000, -0.1667, 4.0000, -0.6667},
        { 0.0000, 0.0000, 0.0000,  0.5000, 0.0000,  2.0000}
    }};

    static inline const std::map<LhaID, std::string> allowed_decays {
        {{511, 311}, "B_K"},
        {{521, 321}, "B_K"},
        {{511, 313}, "B_Ks"},
        {{521, 323}, "B_Ks"},
        {{531, 333}, "B_phi"}
    };

public:
    BaseQCDfCalculator() = default;
    BaseQCDfCalculator(int B_id, int V_id, double mu_b, const std::map<WCoef, complex_t>& C, B_FF_Type ff_tp);

    void fill_wilson_bar_cache();

    complex_t Y(double q2);

protected:
    double gamma_perp(int n) { return 4. * ObsQCDProxy().get_constants()->C_F * (psi(n + 1) + GAMMA - 1. + 1. / (n + 1)); };
    double gamma_par(int n) { return 4. * ObsQCDProxy().get_constants()->C_F * (psi(n + 2) + GAMMA - .75 - 1. /(2. * (n + 1) * (n + 2))); };

    double E(double q2);
    double phi_X(double u, double a1, double a2);
    double gv_dga_4(double u);
    complex_t F_V(double v, bool bar);
    complex_t Y_u(double q2);
    double L(double q2);
    complex_t inv_lambda_B_m(double q2);

    complex_t C_perp_0(double q2, double sign, bool bar);
    complex_t C_par_0(double q2, double sign, bool bar);
    complex_t C_perp_f(double q2, double sign, bool bar);
    complex_t C_par_f(double q2, double sign, bool bar);
    complex_t C_perp_nf(double q2, bool bar);
    complex_t C_par_nf(double q2, bool bar);

    complex_t t_perp(double u, double m_q, double q2, double E_Kstar);
    complex_t t_par(double u, double m_q, double q2, double E_Kstar);
    complex_t T_par_m_0(bool bar);
    complex_t T_par_p_p_f(double u, double q2, bool bar);
    complex_t T_par_p_m_f(double u, double q2, bool bar);
    complex_t T_perp_p_p_f(double u, double q2, bool bar);
    complex_t T_perp_p_m_f(double u, double q2, bool bar);
    complex_t T_perp_p_nf(double u, double q2, bool bar);
    complex_t T_par_p_nf(double u, double q2, bool bar);
    complex_t T_par_m_nf(double u, double q2, bool bar);

    complex_t I_perp_p(double q2, bool bar); 
    complex_t I_perp_m(double q2, bool bar);
    complex_t I_par_p(double q2, bool bar);
    complex_t I_par_m(double q2, bool bar);
};

#endif // __BASEQCDFCALCULATOR_H__
