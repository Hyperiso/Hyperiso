#ifndef BDLNUDECAY_H
#define BDLNUDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "ObsWilsonProxy.h"
#include "ObsQCDProxy.h"
#include "DefaultConfig.h"

struct BDlnuDecayCache {
    double G_F;
    double m_e, m_tau;
    double m_B, m_D;
    double tau_B;
    double V11, rho_D2, Delta;

    complex_t C_V, C_S, C_T;
    bool C_V_flag, C_S_flag, C_T_flag;

    double r_D, r_e, r_tau, r_qp, r_qm;
    double w_e, w_tau;
    double BR_pref;
    double Gamma_p, Gamma_m;
};

/**
 * @brief Decay parent for the B > D l nu_l decays [1309:0301]. Currently implements :
 *     - BR(B > D0 tau nu_tau)
 *     - R_D = Γ(B > D0 tau nu_tau) / Γ(B > D0 e nu_e)
 *     - Forward-Backward asymmetry 
 *     - P_tau = (Γ(s_tau = +1/2) - Γ(s_tau = -1/2)) / Γ
 */
class BDlnuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    BDlnuDecayCache cache;

protected:
    double t(double w);
    double lambda_D(double w);
    double x_l(double rl, double w);
    double phi(double rl, double w);
    double w_max(double rD, double r_l);

    double V_1(double w);
    double S_1(double w);

    double H_V0 (double w);
    double H_Vt (double w);
    double H_S  (double w);
    double H_T  (double w);

    double F_V0_1  (double r_l, double w_m);
    double F_V0_2  (double r_l, double w_m);
    double F_Vt    (double r_l, double w_m);
    double F_S     (double r_l, double w_m);
    double F_T_1   (double r_l, double w_m);
    double F_T_2   (double r_l, double w_m);
    double G_V0_Vt (double r_l, double w_m);
    double G_V0_S  (double r_l, double w_m);
    double G_V0_T  (double r_l, double w_m);
    double G_Vt_S  (double r_l, double w_m);
    double G_Vt_T  (double r_l, double w_m);
    double G_S_T   (double r_l, double w_m);

    double gamma_p(double r_l, double w_m);
    double gamma_m(double r_l, double w_m);

    double BR();
    double A_FB();
    double R_D();
    double P_tau();

public:
    BDlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__D_l_nu), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::CC_bc)};
        this->max_order = QCDOrder::LO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
};

#endif // __BDLNUDECAY_H__