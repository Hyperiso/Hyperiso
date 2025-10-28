#ifndef __M0_MIXING_H__
#define __M0_MIXING_H__

#include "DecayParent.h"
#include "General.h"
#include "ObsParameterMutator.h"
#include "ObsQCDProxy.h"

struct M0MixingCache {
    double G_F;
    double m_W, mu_W;
    double x_t, x_c;
    complex_t lambda_t, lambda_c;
    double m_Bd, m_Bs, m_K, m_D;
    double tau_D;
    complex_t G12_s;
    double delta_G_s;
    double kappa_e, eta_cc, eta_ct;
    double alpha_s_mu_W, alpha_s_mu_b, alpha_s_mu_c;

    std::array<double, 5> Q_Bd;
    std::array<double, 5> Q_Bs;
    std::array<double, 5> Q_K;
    std::array<double, 5> Q_D;
    std::array<complex_t, 8> C_Bd;
    std::array<complex_t, 8> C_Bs;
    std::array<complex_t, 8> C_K;
    std::array<complex_t, 8> C_D;
    complex_t C1_Bd_SM;
    complex_t C1_Bs_SM;
};

/**
 * @brief Parent for the M0 - M0_bar mixing observables. Currently implements TODO. 
 */
class M0Mixing : public DecayParentConfigurable<DecayConfig> {
private:
    M0MixingCache cache;

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

    double Q_i(int i, double B_i, double r_chi, double mf2, bool is_B=false);
    void populate_Q_from_bag(std::array<double, 5>& Q, const std::array<double, 5>& B, double r_chi, double mf2, bool is_B);
    void populate_C(std::array<complex_t, 8>& C, double hadronic_scale, size_t offset=0);
    complex_t M_12_NP(const std::array<complex_t, 8>& C, const std::array<double, 5>& Q, double m_M);
    complex_t M_12_B_SM(int gen);
    complex_t M_12_K_SM();

    double phi_q(int gen);
    double delta_M_B(int gen);
    double a_fs();
    double epsilon_K();
    double delta_M_K();
    double x_D();

public:
    M0Mixing(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::M0_Mix), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::MESON_MIXING)};
        this->max_order = QCDOrder::LO;
        // this->load_params();
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
};

#endif // __M0_MIXING_H__
