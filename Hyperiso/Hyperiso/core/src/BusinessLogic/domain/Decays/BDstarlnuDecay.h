#ifndef BDSTARLNUDECAY_H
#define BDSTARLNUDECAY_H

#include "DecayParent.h"
#include "General.h"

/**
 * @brief Decay parent for the B > D* l nu_l decays [1309:0301]. Currently implements :
 *     - BR(B > D*0 tau nu_tau) 
 *     - R(D*) = Γ(B > D*0 tau nu_tau) / Γ(B > D*0 e nu_e)
 *     - Forward-Backward asymmetry  
 *     - P_tau = (Γ(s_tau = +1/2) - Γ(s_tau = -1/2)) / Γ
 *     - P_D* = Γ(l_D = 0) / Γ
 */
class BDstarlnuDecay : public DecayParent {

protected:
    double ckm  (scalar_t V_cb);
    double pref (double G_F, double tau_B, double m_B, double m_D, double h_A1_1);

    double t        (double rD, double w);
    double lambda_D (double rD, double w);
    double x_l      (double rl, double rD, double w);
    double phi      (double rl, double rD, double w, double rho_D2);
    double w_max    (double rD, double rl);

    double h_A1 (double w, double rho_D2);
    double R_1  (double w, double R_11);
    double R_2  (double w, double R_21);
    double R_3  (double w);

    double H_Vp (double w, double rt_rD, double R_11);
    double H_Vm (double w, double rt_rD, double R_11);
    double H_V0 (double w, double rD, double rt_rD, double R_21);
    double H_Vt (double w, double rD, double rt_rD, double one_m_rD_sq, double R_21);
    double H_S  (double w, double rD, double rt_rD, double one_m_rD_sq, double rqp, double R_21);
    double H_Tp (double w, double rD, double rt_rD, double rqp, double rqm, double R_11);
    double H_Tm (double w, double rD, double rt_rD, double rqp, double rqm, double R_11);
    double H_T0 (double w, double rD, double rt_rD, double one_m_rD_sq, double rqp, double rqm, double R_11);

    double F_V0_1  (double rD, double rt_rD, double rl, double rho_D2, double R_21, double w_m, bool flag);
    double F_V0_2  (double rD, double rt_rD, double rl, double rho_D2, double R_21, double w_m, bool flag);
    double F_Vp_1  (double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag);
    double F_Vp_2  (double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag);
    double F_Vm_1  (double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag);
    double F_Vm_2  (double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag);
    double F_Vt    (double rD, double rt_rD, double one_m_rD_sq, double rl, double rho_D2, double R_21, double w_m, bool flag);
    double F_S     (double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rho_D2, double R_21, double w_m, bool flag);
    double F_T0_1  (double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double F_T0_2  (double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double F_Tp_1  (double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double F_Tp_2  (double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double F_Tm_1  (double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double F_Tm_2  (double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);

    double G_Vp_Vm_1 (double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag);
    double G_Vp_Vm_2 (double rD, double rt_rD, double rl, double rho_D2, double R_11, double w_m, bool flag);
    double G_Vp_Tp   (double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double G_Vp_Tm   (double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double G_Vm_Tp   (double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double G_Vm_Tm   (double rD, double rt_rD, double rl, double rqp, double rqm, double rho_D2, double R_11, double w_m, bool flag);
    double G_V0_Vt   (double rD, double rt_rD, double one_m_rD_sq, double rl, double rho_D2, double R_21, double w_m, bool flag);
    double G_V0_S    (double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rho_D2, double R_21, double w_m, bool flag);
    double G_V0_T0   (double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double R_21, double w_m, bool flag);
    double G_Vt_S    (double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rho_D2, double R_21, double w_m, bool flag);
    double G_Vt_T0   (double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double R_21, double w_m, bool flag);
    double G_S_T0    (double rD, double rt_rD, double one_m_rD_sq, double rl, double rqp, double rqm, double rho_D2, double R_11, double R_21, double w_m, bool flag);

    scalar_t C_V1();
    scalar_t C_V2();
    scalar_t C_A();
    scalar_t C_P();
    scalar_t C_T();

    double c_flag(scalar_t C);

    double Gamma_tau_m(double F_Vp_1, double F_Vm_1, double F_V0_1, double F_Tp_2, double F_Tm_2, double F_T0_2, double G_Vp_Vm_1, double G_T0_V0, double G_Tp_Vp, double G_Tm_Vm, double G_Tp_Vm, double G_Tm_Vp, scalar_t C_V1, scalar_t C_V2, scalar_t C_T);
    double Gamma_tau_p(double F_Vp_2, double F_Vm_2, double F_V0_2, double F_Vt, double F_S, double F_Tp_1, double F_Tm_1, double F_T0_1, double G_Vp_Vm_2, double G_Vt_S, double G_T0_V0, double G_Tp_Vp, double G_Tm_Vm, double G_Tp_Vm, double G_Tm_Vp, scalar_t C_V1, scalar_t C_V2, scalar_t C_A, scalar_t C_P, scalar_t C_T);
    double Gamma_D_0(double F_V0_1, double F_V0_2, double F_Vt, double F_S, double F_T0_1, double F_T0_2, double G_Vt_S, double G_V0_T0, scalar_t C_A, scalar_t C_P, scalar_t C_T);
    double B_theta(double F_Vp_1, double F_Vm_1, double F_Tp_2, double F_Tm_2, double G_V0_Vt, double G_V0_S, double G_Vt_T0, double G_Tp_Vp, double G_Tm_Vm, double G_Tp_Vm, double G_Tm_Vp, double G_T0_S, scalar_t C_V1, scalar_t C_V2, scalar_t C_A, scalar_t C_P, scalar_t C_T);

private:
    const QCDOrder max_order = QCDOrder::LO;

public:
    BDstarlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<IObsWilsonBuilder<ObsWilsonProxy, WGroup>> wilson_builder) : DecayParent(matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::BCLNU};
    }

    void build_op_tree() override;

};

#endif // __BDSTARLNUDECAY_H__