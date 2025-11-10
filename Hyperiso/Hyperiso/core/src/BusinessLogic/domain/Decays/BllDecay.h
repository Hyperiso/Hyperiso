#ifndef BLLDECAY_H
#define BLLDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"

struct BllDecayCache {
    double G_F;
    double alpha_em;
    double m_mu;
    double m_Bd;
    double m_Bs;
    double f_Bd;
    double f_Bs;
    double tau_Bd;
    double tau_Bs;
    complex_t lambda_d;
    complex_t lambda_s;
    double ys;
    double eta_BBS;
    double x_d;
    double x_s;
    double r_d;
    double r_s;
    double beta_d;
    double beta_s;

    complex_t C10_SM;
    complex_t C10;
    complex_t CQ1;
    complex_t CQ2;
    complex_t C10_m;
    complex_t CQ1_m;
    complex_t CQ2_m;
};

/**
 * @brief Decay parent for the Bq > ll decays. Currently implements both CP-averaged and untagged Bs > mu+ mu- branching ratios and the CP-averaged Bd > mu+ mu- decays. 
 */
class BllDecay : public DecayParentConfigurable<DecayConfig> {
private:
    BllDecayCache cache;

protected:
    double BR_avg_Bq_mumu(int q);
    double BR_untag_Bs_mumu();

public:
    BllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__l_l), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime), GroupMapper::to_id(WGroup::BScalar)};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
};

#endif // __BLLDECAY_H__