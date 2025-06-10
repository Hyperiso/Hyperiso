#ifndef BLNUDECAY_H
#define BLNUDECAY_H

#include "DecayParent.h"
#include "General.h"

/**
 * @brief Decay parent for the B > l nu_l decays. Currently implements Bu > tau nu_tau branching ratio and the R_nu_tau ratio of this BR its SM expectation. 
 */
class BlnuDecay : public DecayParent {

protected:
    scalar_t R(double m_B, double m_b, double m_tau);
    double ckm(scalar_t V_ub);
    double pref(double G_F, double f_B, double tau_B, double m_B, double m_tau);
    double BR_B_taunu(double pref, double ckm, double R);

public:
    BlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParent(matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::BCC};
        this->max_order = QCDOrder::LO;
    }

    void build_op_tree() override;

};

#endif // __BLNUDECAY_H__