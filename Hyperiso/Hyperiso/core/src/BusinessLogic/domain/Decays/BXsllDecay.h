#ifndef __BXSLLDECAY_H__
#define __BXSLLDECAY_H__

#include "DecayParent.h"
#include "ObsQCDProxy.h"
#include "Include.h"
#include "Math.h"

struct kappa_coef {
    size_t i;
    size_t j;
    int l;
    size_t m;
    double value;
};

/**
 * @brief Decay parent for the inclusive B > X_s l+ l- decays. Implements the integrated branching ratio in both q² \in [1, 6] GeV² and q² > 14.4 GeV², as well as the forward-backward asymmetry of the decay. 
 */
class BXsllDecay : public DecayParent {

private:
    static inline const std::array<double, 6> cc_res_mass      {3.096916, 3.68609,   3.77292, 4.039 , 4.153  , 4.421 };
    static inline const std::array<double, 6> cc_res_br        {5.93e-2 , 7.7e-3 ,   1.1e-5 , 1.4e-5, 1.0e-5 , 1.1e-5};
    static inline const std::array<double, 6> cc_res_width_tot {9.29e-5 , 3.04e-4,   2.73e-2, 8.0e-2, 1.03e-1, 6.2e-2};
    static inline const std::array<double, 6> cc_res_width_had {8.147e-5, 2.9746e-4, 2.36e-2, 5.2e-2, 7.8e-2 , 4.3e-2};

    static const std::vector<kappa_coef> kappa_17_low;
    static const std::vector<kappa_coef> kappa_17_high;


protected:
    scalar_t alpha_s(scalar_t mu);
    scalar_t f(scalar_t z);
    scalar_t h(scalar_t z);
    scalar_t kappa(scalar_t f, scalar_t h, scalar_t alpha_s_mu_b);
    
    scalar_t m_hat(scalar_t m, scalar_t mb);
    scalar_t z(scalar_t mc_hat);

    scalar_t f_7(scalar_t s);
    scalar_t f_9(scalar_t s);
    scalar_t tau_77(scalar_t s);
    scalar_t tau_99(scalar_t s);
    scalar_t tau_79(scalar_t s);
    scalar_t tau_710(scalar_t s);
    scalar_t tau_910(scalar_t s);
    scalar_t sigma(scalar_t s);
    scalar_t sigma_9(scalar_t s);
    scalar_t sigma_7(scalar_t s, scalar_t L_mu);

    scalar_t g(scalar_t z, scalar_t s);
    scalar_t R_cc_cont(scalar_t s);
    scalar_t breit_wigner(scalar_t s, scalar_t m_V, scalar_t br, scalar_t gamma_tot, scalar_t gamma_had);
    scalar_t R_cc(scalar_t s, scalar_t inv_alpha_em);
    scalar_t g_ld(scalar_t z, scalar_t s, scalar_t inv_alpha_em, scalar_t m_D_hat);
    scalar_t C9_eff(scalar_t s, scalar_t L_mu, scalar_t mc_hat, scalar_t inv_alpha_em, scalar_t m_D_hat, QCDOrder order);

    scalar_t F_17(scalar_t L_b, scalar_t z, scalar_t s);
    scalar_t F_27(scalar_t L_b, scalar_t z, scalar_t s);
    scalar_t F_19(scalar_t L_b, scalar_t z, scalar_t s);
    scalar_t F_29(scalar_t L_b, scalar_t z, scalar_t s);
    scalar_t C7_new(scalar_t s, scalar_t alpha_s_mu_b, scalar_t C7_eff, scalar_t C8_eff_0, scalar_t L_mu, scalar_t L_c, scalar_t z);
    scalar_t C9_new(scalar_t s, scalar_t alpha_s_mu_b, scalar_t C8_eff_0, scalar_t L_mu, scalar_t L_c, scalar_t z);
    scalar_t C10_new(scalar_t s, scalar_t alpha_s_mu_b);


    scalar_t H7(scalar_t s, scalar_t ml_hat, scalar_t alpha_s_mu_b);
    scalar_t H9(scalar_t s, scalar_t ml_hat, scalar_t alpha_s_mu_b);
    scalar_t H10(scalar_t s, scalar_t ml_hat, scalar_t alpha_s_mu_b);
    scalar_t H79(scalar_t s, scalar_t ml_hat, scalar_t alpha_s_mu_b);

    scalar_t dB_ds(scalar_t s, scalar_t ml_hat, scalar_t alpha_s_mu_w, scalar_t );

public:
    BXsllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParent(matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar};
        this->max_order = QCDOrder::NNLO;
    }

    void build_op_tree() override;

};

#endif // __BXSLLDECAY_H__
