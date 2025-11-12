#ifndef __DSLNUDECAY_H__
#define __DSLNUDECAY_H__

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"
#include "PlnuCalculator.h"

struct DslnuDecayCache {
    PlnuCalculator calc_mu;
    PlnuCalculator calc_tau;
};

/**
 * @brief Decay parent for the B > l nu_l decays. Currently implements Bu > tau nu_tau branching ratio and the R_nu_tau ratio of this BR its SM expectation. 
 */
class DslnuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    DslnuDecayCache cache;

protected:
    double BR(int gen);

public:
    DslnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::Ds__l_nu), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::CC_cs)};
        this->max_order = QCDOrder::LO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

};

#endif // __DSLNUDECAY_H__
