#ifndef BDLNUDECAY_H
#define BDLNUDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "WilsonAdapter.h"
#include "ObsWilsonProxy.h"

/**
 * @brief Decay parent for the B > D l nu_l decays [1309:0301]. Currently implements :
 *     - BR(B > D0 tau nu_tau)
 *     - R_D = Γ(B > D0 tau nu_tau) / Γ(B > D0 e nu_e)
 *     - Forward-Backward asymmetry 
 *     - P_tau = (Γ(s_tau = +1/2) - Γ(s_tau = -1/2)) / Γ
 */
class BDlnuDecay : public DecayParent {
protected:
    double ckm(complex_t V_cb);
    double pref(double G_F, double tau_B, double m_B, double m_D, double V_11);

    double t(double rD, double w);
    double lambda_D(double rD, double w);
    double x_l(double rtau, double rD, double w);
    double phi(double rtau, double rD, double w);
    double w_max(double rD, double rtau);

    double V_1(double w, double rho_D2);
    double S_1(double w, double rho_D2, double Delta);

    double H_V0 (double w, double rD, double rho_D2);
    double H_Vt (double w, double rD, double rho_D2, double Delta);
    double H_S  (double w, double rD, double rqm, double rho_D2, double Delta);
    double H_T  (double w, double rD, double rqp, double rho_D2, double Delta);

    double F_V0_1  (double rD, double rtau, double rho_D2, double w_m, bool flag);
    double F_V0_2  (double rD, double rtau, double rho_D2, double w_m, bool flag);
    double F_Vt    (double rD, double rtau, double rho_D2, double Delta, double w_m, bool flag);
    double F_S     (double rD, double rtau, double rqm, double rho_D2, double Delta, double w_m, bool flag);
    double F_T_1   (double rD, double rtau, double rqp, double rho_D2, double Delta, double w_m, bool flag);
    double F_T_2   (double rD, double rtau, double rqp, double rho_D2, double Delta, double w_m, bool flag);
    double G_V0_Vt (double rD, double rtau, double rho_D2, double Delta, double w_m, bool flag);
    double G_V0_S  (double rD, double rtau, double rqm, double rho_D2, double Delta, double w_m, bool flag);
    double G_V0_T  (double rD, double rtau, double rqp, double rho_D2, double Delta, double w_m, bool flag);
    double G_Vt_S  (double rD, double rtau, double rqm, double rho_D2, double Delta, double w_m, bool flag);
    double G_Vt_T  (double rD, double rtau, double rqp, double rho_D2, double Delta, double w_m, bool flag);
    double G_S_T   (double rD, double rtau, double rqp, double rqm, double rho_D2, double Delta, double w_m, bool flag);

    complex_t C_V();
    complex_t C_S();
    complex_t C_T();

    double c_flag(complex_t C);

    double Gamma_m(double F_V0_1, double F_T_2, double G_V0_T, complex_t C_V, complex_t C_T);
    double Gamma_p(double F_V0_2, double F_Vt, double F_S, double F_T_1, double G_Vt_S, double G_V0_T, complex_t C_V, complex_t C_S, complex_t C_T);
    double Gamma(double gamma_p, double gamma_m);
    double B_theta(double G_V0_Vt, double G_V0_S, double G_Vt_T, double G_S_T, complex_t C_V, complex_t C_S, complex_t C_T);

    double BR_B_Dtaunu(double pref, double ckm, double width);

public:
    BDlnuDecay(QCDOrder order, double matching_scale, double hadronic_scale) {
        order = check_max_order(QCDOrder::LO);
        WilsonAdapter().build({WGroup::BCLNU}, matching_scale, hadronic_scale, order);
        this->w_proxy = std::make_shared<ObsWilsonProxy>();
        build_op_tree();
    }

    void build_op_tree() override;

};

#endif // __BDLNUDECAY_H__