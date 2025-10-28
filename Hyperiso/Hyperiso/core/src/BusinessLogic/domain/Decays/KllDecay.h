#ifndef KLLDECAY_H
#define KLLDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"

struct KllDecayConfig {
    int N_L_sign = 1;
};

struct KllDecayCache {
    double G_F;
    double alpha_em;
    double alpha_em_0;
    double alpha_s_m_Z;
    double sw2;
    double m_mu;
    double m_c_m_c;
    double m_pi;
    double m_rho;
    double m_K;
    double f_K;
    double tau_L;
    double tau_S;
    complex_t lambda_c;
    complex_t lambda_t;
    double lambda;
    double x;
    double r_chi;
    double beta;
    double mu_b;

    double BR_KL_gg_exp;
    double BR_KS_gg_exp;
    double alpha_exp;
    double delta_lambda;

    complex_t C10;
    complex_t CQ1;
    complex_t CQ2;
};

/**
 * @brief Decay parent for the K_L,S > ll decays. Currently implements BR(K_L,S > mu mu) with both possible signs for the K_L > 2 gamma long distance correction.
 */
class KllDecay : public DecayParentConfigurable<KllDecayConfig> {
private:
    KllDecayCache cache;
    KllDecayConfig cfg;

protected:  
    complex_t C_gg(double b);
    complex_t N_L();
    complex_t N_S();
    double P_c();
    double BR_L();
    double BR_S();

public:
    KllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::K__l_l), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::K)};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

    void set_config_spe(KllDecayConfig config) override {this->cfg = config;}
};

#endif // __KLLDECAY_H__