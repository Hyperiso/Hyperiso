#ifndef __BDLNUDECAY_H__
#define __BDLNUDECAY_H__

#include "DecayParent.h"
#include "General.h"

/**
 * @brief Decay parent for the B > D l nu_l decays. Currently implements the B > D0 tau nu_tau branching ratio and the xi_Dlnu = BR(B > D0 tau nu_tau) / BR(B > D0 e nu_e) ratio.
 */
class BDlnuDecay : public DecayParent {

protected:
    double ckm(double V_cb_r, double V_cb_i);
    double pref(double G_F, double tau_B, double m_B, double m_D, double G1);
    double t(double rD, double w);
    double G(double rho2, double w);
    double f_1(double rD, double rl, double rho2, double w);
    double f_2(double rD, double rl, double rho2, double w);
    double f_3(double rD, double rl, double rho2, double w);
    double f_4(double rD, double rl, double rho2, double w);
    double w_max(double rD, double rl);
    double I1(double rD, double rl, double rho2, double w_u);
    double I2(double rD, double rl, double rho2, double w_u);
    double I3(double rD, double rl, double rho2, double w_u);
    double I4(double rD, double rl, double rho2, double w_u);
    double W1(double I1);
    double W2(double I2, double Delta2, double rl);
    double W3(double I3, double Delta2, double m_l, double m_bc);
    double W4(double I4, double Delta2, double m_B, double m_bc);
    double W(double W1, double W2, double W3, double W4);
    double BR_B_Dlnu(double pref, double ckm, double W);

public:
    BDlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale) {
        winfo.matching_scale = matching_scale;
        winfo.hadronic_scale = hadronic_scale;
        winfo.model = MemoryManager::GetInstance()->getModel();
        winfo.order = order;
        winfo.basis = BWilsonBasis::STANDARD;
        winfo.wgroups = {WGroup::Blnu};

        max_order = QCDOrder::LO;

        build_op_tree();
    }

    void build_op_tree() override;

};

#endif // __BDLNUDECAY_H__