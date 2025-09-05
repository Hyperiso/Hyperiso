#ifndef __M0_MIXING_H__
#define __M0_MIXING_H__

#include "DecayParent.h"
#include "General.h"
#include "ObsParameterMutator.h"
#include "ObsQCDProxy.h"

/**
 * @brief Parent for the M0 - M0_bar mixing observables. Currently implements TODO. 
 */
class M0Mixing : public DecayParent {

private:
    std::array<double, 5> Q;
    std::array<double, 5> B;
    std::array<complex_t, 8> C;

    static constexpr std::array<double, 5> N_i {0.66667, -0.41667, 0.08333, 0.50000, 0.16667};
    static constexpr std::array<double, 5> d_i {0.00000,  0.00000, 0.00000, 0.16667, 1.50000};
    static constexpr double J_5 = 5165. / 3174.;
    static constexpr double J_4 = 6719. / 3750.;
    static constexpr double J_3 = 307. / 162.;
    static constexpr double B_t = 6.5;

protected:
    double S_18(double x);
    double S_11(double x);
    double S_1t(double x);
    double dS0_dx(double x);
    double S0_ct(double x_c, double x_t);
    double F_S1(double x, double m_W, double mu_W, double mu_t);

    double alpha(double mu);
    double x_t(double mu, double m_W);
    double c_RGI_B(double alpha_s_mu_b);
    double eta_2B(double alpha_s_mu_W, double x_t, double m_W);

    complex_t lambda(complex_t V_qd, complex_t V_qs);
    double c_RGI_K(double alpha_s_mu_c);
    double eta_tt(double alpha_s_mu_c, double alpha_s_mu_b, double alpha_s_mu_W, double x_t, double m_W);
    complex_t C_sd_SM(double x_t, double x_c, complex_t lambda_c, complex_t lambda_t, double eta_cc, double eta_tt, double eta_ct);

    double Q_i(int i, double B_i, double r_chi, double mf2, bool is_B=false);
    void populate_B(double B_1, double B_2, double B_3, double B_4, double B_5);
    void populate_Q_from_bag(double r_chi, double mf2, bool is_B);
    void populate_Q_from_mels(double Q_1, double Q_2, double Q_3, double Q_4, double Q_5);
    void populate_C(double hadronic_scale, size_t offset=0);
    complex_t M_12_NP(double m_M);
    complex_t M_12_B_SM(double m_B, double eta_2B, double c_RGI_B, int gen);
    complex_t M_12_K_SM(double m_M, double G_F, double m_W, double c_RGI_K, complex_t C_1_sd);

    double phi_q(complex_t M_12);
    double delta_M_B(complex_t M_12);
    double a_fs(complex_t M_12, complex_t G_12, double delta_M, double delta_G);
    double epsilon_K(complex_t M_12, double kappa_e);
    double delta_M_K(complex_t M_12);
    double x_D(complex_t M_12, double tau_D0);

public:
    M0Mixing(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParent(matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::MESON_MIXING};
        this->max_order = QCDOrder::LO;
    }

    void build_op_tree() override;
};

#endif // __M0_MIXING_H__
