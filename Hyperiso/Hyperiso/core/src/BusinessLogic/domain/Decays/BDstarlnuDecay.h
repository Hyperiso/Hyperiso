#ifndef BDSTARLNUDECAY_H
#define BDSTARLNUDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"
#include "ObsQCDProxy.h"

struct BDstarlnuDecayCache {
    double G_F;
    double m_e, m_tau;
    double m_B, m_D_star;
    double tau_B;
    double h_A1_1, rho_D2, R_11, R_21;

    complex_t C_V1, C_V2, C_A, C_P, C_T;
    bool C_V1_flag, C_V2_flag, C_A_flag, C_P_flag, C_T_flag;

    double r_D, r_e, r_tau, r_qp, r_qm;
    double sqrt_rD, one_m_rD2;
    double w_e, w_tau;
    double BR_pref;
    double Gamma_p, Gamma_m;
};

/**
 * @brief Decay parent for the B > D* l nu_l decays [1309:0301]. Currently implements :
 *     - BR(B > D*0 tau nu_tau) 
 *     - R(D*) = Γ(B > D*0 tau nu_tau) / Γ(B > D*0 e nu_e)
 *     - Forward-Backward asymmetry  
 *     - P_tau = (Γ(s_tau = +1/2) - Γ(s_tau = -1/2)) / Γ
 *     - P_D* = Γ(l_D = 0) / Γ
 */
class BDstarlnuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    BDstarlnuDecayCache cache;

protected:
    // Kinematics
    double t        (double w);
    double lambda_D (double w);
    double x_l      (double rl, double w);
    double phi      (double rl, double w);
    double w_max    (double rl);

    // Form Factors
    double h_A1 (double w);
    double R_1  (double w);
    double R_2  (double w);
    double R_3  (double w);

    // Helicity amplitudes
    double H_Vp (double w);
    double H_Vm (double w);
    double H_V0 (double w);
    double H_Vt (double w);
    double H_S  (double w);
    double H_Tp (double w);
    double H_Tm (double w);
    double H_T0 (double w);

    // Helicity amplitude integrals
    double F_V0_1  (double rl, double w_m);
    double F_V0_2  (double rl, double w_m);
    double F_Vp_1  (double rl, double w_m);
    double F_Vp_2  (double rl, double w_m);
    double F_Vm_1  (double rl, double w_m);
    double F_Vm_2  (double rl, double w_m);
    double F_Vt    (double rl, double w_m);
    double F_S     (double rl, double w_m);
    double F_T0_1  (double rl, double w_m);
    double F_T0_2  (double rl, double w_m);
    double F_Tp_1  (double rl, double w_m);
    double F_Tp_2  (double rl, double w_m);
    double F_Tm_1  (double rl, double w_m);
    double F_Tm_2  (double rl, double w_m);

    double G_Vp_Vm_1 (double rl, double w_m);
    double G_Vp_Vm_2 (double rl, double w_m);
    double G_Vp_Tp   (double rl, double w_m);
    double G_Vp_Tm   (double rl, double w_m);
    double G_Vm_Tp   (double rl, double w_m);
    double G_Vm_Tm   (double rl, double w_m);
    double G_V0_Vt   (double rl, double w_m);
    double G_V0_S    (double rl, double w_m);
    double G_V0_T0   (double rl, double w_m);
    double G_Vt_S    (double rl, double w_m);
    double G_Vt_T0   (double rl, double w_m);
    double G_S_T0    (double rl, double w_m);

    // Intermediate results
    double Gamma_tau_m(double rl, double w_m);
    double Gamma_tau_p(double rl, double w_m);

    // Observables
    double BR();
    double A_FB();
    double R_Dstar();
    double P_tau();
    double P_D();

public:
    BDstarlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__Dstar_l_nu), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::BCC)};
        this->max_order = QCDOrder::LO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

};

#endif // __BDSTARLNUDECAY_H__