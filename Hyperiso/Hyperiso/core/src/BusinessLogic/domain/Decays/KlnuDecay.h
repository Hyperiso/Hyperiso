#ifndef __KLNU_H__
#define __KLNU_H__

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"
#include "PlnuCalculator.h"

struct KlnuDecayCache {
    double delta_em;
    double f_corr;

    PlnuCalculator calc_K;
    PlnuCalculator calc_pi;
};

/**
 * @brief Decay parent for the B > l nu_l decays. Currently implements Bu > tau nu_tau branching ratio and the R_nu_tau ratio of this BR its SM expectation. 
 */
class KlnuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    KlnuDecayCache cache;

protected:
    double R_mu23();
    double BR_K_BR_pi();

public:
    KlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::K__l_nu), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::CC_su), GroupMapper::to_id(WGroup::CC_du)};
        this->max_order = QCDOrder::LO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

};

#endif // __KLNU_H__
