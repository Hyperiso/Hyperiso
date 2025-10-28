#ifndef BLNUDECAY_H
#define BLNUDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"

struct BlnuDecayCache {
    double G_F;
    double m_tau;
    double m_b;
    double V_ub_2;
    double tau_B;
    double m_B;
    double f_B;
    complex_t C_V;
    complex_t C_S;
};

/**
 * @brief Decay parent for the B > l nu_l decays. Currently implements Bu > tau nu_tau branching ratio and the R_nu_tau ratio of this BR its SM expectation. 
 */
class BlnuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    BlnuDecayCache cache;

protected:
    scalar_t R();
    double BR();

public:
    BlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__l_nu), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::BCC)};
        this->max_order = QCDOrder::LO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

};

#endif // __BLNUDECAY_H__