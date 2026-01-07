#ifndef __DLNUDECAY_H__
#define __DLNUDECAY_H__

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"
#include "PlnuCalculator.h"

struct DlnuDecayCache {
    PlnuCalculator calc;
};

/**
 * @brief Decay parent for the B > l nu_l decays. Currently implements Bu > tau nu_tau branching ratio and the R_nu_tau ratio of this BR its SM expectation. 
 */
class DlnuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    DlnuDecayCache cache;

protected:
    double BR();

public:
    DlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, ObservablePortsConfig& ports) : DecayParentConfigurable(DecayMapper::to_id(Decays::D__l_nu), matching_scale, hadronic_scale, order, ports) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::CC_cd)};
        this->max_order = QCDOrder::LO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

};

#endif // __DLNUDECAY_H__
