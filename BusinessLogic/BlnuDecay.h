#ifndef __BLNUDECAY_H__
#define __BLNUDECAY_H__

#include "DecayParent.h"
#include "General.h"

/**
 * @brief Decay parent for the B > l nu_l decays. Currently implements Bu > tau nu_tau branching ratio and the R_nu_tau ratio of this BR its SM expectation. 
 */
class BlnuDecay : public DecayParent {

protected:
    complex_t R(double m_B, double m_b, double m_tau);
    double ckm(double V_ub_r, double V_ub_i);
    double pref(double G_F, double f_B, double tau_B, double m_B, double m_tau);
    double BR_B_taunu(double pref, double ckm, double R);

public:
    BlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale) {
        winfo.matching_scale = matching_scale;
        winfo.hadronic_scale = hadronic_scale;
        winfo.model = MemoryManager::GetInstance()->getModel();
        winfo.order = order;
        winfo.basis = BWilsonBasis::STANDARD;
        winfo.wgroups = {WGroup::Blnu};

        build_op_tree();
    }

    void build_op_tree() override;

};

#endif // __BLNUDECAY_H__