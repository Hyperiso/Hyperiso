#ifndef BXSDECAY_H
#define BXSDECAY_H

#include "DecayParent.h"
#include "General.h"
#include <array>
#include "DefaultConfig.h"
#include "ObsQCDProxy.h"

struct BXsDecayCache {
    double alpha_em, alpha_em_0;
    double m_b_mb, m_b_kin, r_msbar_1S;
    double m_c, m_s, m_W;
    double ckm_factor;
    double mu_b, mu_W;
    double beta_0;
    double alpha_s_mu_b, alpha_s_upsilon, eta;
    double E0;
    double BR_B__Xc_e_nu_exp;
    double mu_G2, rho_D3, rho_LS3;
    double lambda_2;
    double mu_c;
    double z0, z1;
    double z, delta;
    double m_c_mu_c, m_c_3gev;
    double L_b, L_c;
    const double X_b {-0.168440839756};

    std::map<WCoef, complex_t> C_b_LO;
    std::map<WCoef, complex_t> C_b_NLO;
    std::map<WCoef, complex_t> C_b_NNLO;
    std::map<WCoef, complex_t> C_w;
};

/**
 * @brief Decay parent for the Bq > ll decays. Currently implements both CP-averaged and untagged Bs > mu+ mu- branching ratios and the CP-averaged Bd > mu+ mu- decays. 
 */
class BXsDecay : public DecayParentConfigurable<DecayConfig> {
private:
    BXsDecayCache cache;

protected:
    static constexpr std::array<double, 8> gamma_i7 {-0.85596, 5.1358, -2.1728, -0.6255, -77.432, 19.029, 10.667, -3.5556};
    static constexpr std::array<double, 8> a_i {0.608696, 0.695652, 0.2608696, -0.521739, 0.4086, -0.4230, -0.8994, 0.1456};
    static constexpr std::array<double, 8> d_i {1.4107, -0.8380, -0.4286, -0.0714, -0.6494, -0.0380, -0.0185, -0.0057};

    // Utility
    double gen_P00(const std::array<std::array<double, 8>, 8>& K);
    double gen_P01(const std::array<std::array<double, 8>, 8>& K);

    // P_21 and P_32
    double a(double z);
    double b(double z);
    complex_t G(double t);
    double phi_22(double z, double delta); 
    double phi_27(double z, double delta);  
    double phi_47(double delta); 
    double phi_77(double delta); 
    double phi_78(double delta); 
    double phi_88(double delta); 
    std::array<std::array<double, 8>, 8> phi_1(double delta, double z);
    std::array<double, 8> r_1(double z);
    std::array<std::array<double, 8>, 8> K_1();

    // P_22 (beta_0)
    double F2nf(double z);
    double r22(double z);
    double h22(double z, double delta);
    double h27(double z, double delta);
    double h28(double z, double delta);
    double h88(double delta);
    double h77(double delta);

    std::array<std::array<double, 8>, 8> phi_2_b0(double delta, double z);
    std::array<double, 8> r_hat_2(double z);
    std::array<std::array<double, 8>, 8> K_2_b0();

    // P_22 (rem)
    double F2a(double z);
    double F2na(double z);
    double phi_77_rem(double phi_77_int);
    std::array<std::array<double, 8>, 8> K_2_rem(double z);
    double r2_large_z(double z);
    double dr2_dlogz(double z);
    double r22_large_z(double z);
    double P22_rem();

    double P();

    // Nonperturbative corrections
    double N();

    // Electromagnetic corrections
    double C2_em(double eta);
    double C8_em(double eta);
    complex_t C7_em(double eta);
    double epsilon_em();

    // Branching Ratio
    double C();
    double BR_B_Xs_gamma();

public:
    BXsDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__Xs), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime)};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

};

#endif // __BXSDECAY_H__