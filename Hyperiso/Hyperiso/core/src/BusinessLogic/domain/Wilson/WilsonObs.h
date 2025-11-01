#ifndef WILSONDECAY_H
#define WILSONDECAY_H

#include "DecayParent.h"
#include "General.h"

/**
 * @brief "Decay parent" to compute uncertainty on Wilson Coefficients. 
 */
class WilsonDecay : public DecayParentConfigurable<DecayConfig> {

private:
    complex_t get_matching_coef(WGroup group_id, WCoef coef_id);
    complex_t get_running_coef(WGroup group_id, WCoef coef_id);

public:
    WilsonDecay(DecayId id, QCDOrder order, double matching_scale, double hadronic_scale) : DecayParentConfigurable(id, matching_scale, hadronic_scale, order) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime), GroupMapper::to_id(WGroup::BScalar), GroupMapper::to_id(WGroup::CC_bc), GroupMapper::to_id(WGroup::MESON_MIXING)};
        this->max_order = QCDOrder::NNLO;
    }

    void build_op_tree();
    void load_params() override {}

    std::vector<ObservableValue> compute_observable(Observables obs) override {return {};}
    std::vector<ObservableValue> compute_observable(ObservableId obs) override {return {};}
};

#endif // __WILSONDECAY_H__