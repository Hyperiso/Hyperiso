#ifndef BLLDECAY_H
#define BLLDECAY_H

#include "DecayParent.h"
#include "General.h"

/**
 * @brief Decay parent for the Bq > ll decays. Currently implements both CP-averaged and untagged Bs > mu+ mu- branching ratios and the CP-averaged Bd > mu+ mu- decays. 
 */
class BllDecay : public DecayParent {

protected:
    scalar_t W1(scalar_t r, scalar_t CQ1, scalar_t CPQ1, bool prime);
    scalar_t W2Q(scalar_t r, scalar_t CQ2, scalar_t CPQ2, bool prime);
    scalar_t W210(scalar_t x, scalar_t C10, scalar_t CP10, bool prime);
    scalar_t ckm(scalar_t V_tb, scalar_t V_tq);
    scalar_t BR_avg_Bq_mumu(scalar_t w1, scalar_t w2q, scalar_t w210, scalar_t ckm, scalar_t b, scalar_t G_F, scalar_t inv_alpha, scalar_t f_Bs, scalar_t m_Bs, scalar_t life_Bs);
    scalar_t A_DG(scalar_t x, scalar_t r, scalar_t w210, scalar_t w1q, scalar_t w2q, scalar_t C10_SM);
    scalar_t BR_untag_Bs_mumu(scalar_t br_avg, scalar_t ys, scalar_t A);

private:
    const QCDOrder max_order = QCDOrder::NNLO;

public:
    BllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<IObsWilsonBuilder<ObsWilsonProxy, WGroup>> wilson_builder) : DecayParent(matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar};
    }

    void build_op_tree() override;
};

#endif // __BLLDECAY_H__