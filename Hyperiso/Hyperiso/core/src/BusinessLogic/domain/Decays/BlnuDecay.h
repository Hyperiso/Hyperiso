#ifndef BLNUDECAY_H
#define BLNUDECAY_H

#include "DecayParent.h"
#include "DefaultConfig.h"
#include "PlnuCalculator.h"

struct BlnuDecayCache {
    PlnuCalculator calc;
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
    BlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, ObservablePortsConfig& ports) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__l_nu), matching_scale, hadronic_scale, order, ports) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::CC_bu)};
        this->max_order = QCDOrder::LO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

};

#endif // BLNUDECAY_H