#ifndef __BKSTARLLDECAY_H__
#define __BKSTARLLDECAY_H__

#include "DecayParent.h"
#include "Include.h"

/**
 * @brief Decay parent for the inclusive B > K* l+ l- decays. Implements the integrated branching ratio and angular observables in several q² bins, as well as the forward-backward asymmetry of the decay. 
 */
class BKstarllDecay : public DecayParent {

protected:
    

public:
    BKstarllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParent(matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar};
        this->max_order = QCDOrder::NNLO;
    }

    void build_op_tree() override;

};

#endif // __BXSLLDECAY_H__
