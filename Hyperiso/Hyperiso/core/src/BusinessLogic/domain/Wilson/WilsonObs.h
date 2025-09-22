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
        this->w_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar, WGroup::BCC, WGroup::MESON_MIXING};
        this->max_order = QCDOrder::NNLO;
    }

    void build_op_tree() override;
};

#endif // __WILSONDECAY_H__