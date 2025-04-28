#ifndef BKSTARDECAY_H
#define BKSTARDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "QCDHelper.h"
#include "Math.h"
#include "ObsParameterMutator.h"

class BKstarDecay : public DecayParent {

protected:
    double alpha_s(double mu);
    double beta_0(double mu);
    double sc(double mb_mu_b, double hadronic_scale);
    double run(double initial_value, double eta, double gamma, double beta);
    double a_n_perp(int n, double a_1_gev, double beta_0, double eta);
    double a_n_par(int n, double a_1_gev, double beta_0, double eta);
    double lambda_B(double lam_B_1_gev, double mu_h, double alpha_s_mu_h);
    double f_Ks_perp(double f_1_gev, double beta_0, double eta);
    complex_t h(double s, double u);
    complex_t g_2(double s);
    double phi_perp(double a1, double a2, double u);
    double gv_dga_4(double a1, double a2, double z3a, double z3v, double w10a, double dtp, double dtm, double u);
    complex_t G(double s, double xbar);
    double F_perp(double a1, double a2); 
    complex_t G_perp(double s, double a1, double a2); 
    complex_t H_perp(double s, double a1par, double a2par, double z3a, double z3v, double w10a, double dtp, double dtm); 
    double X_perp(double a1, double a2, double m_B, double Lambda_h); 
    complex_t G2(double s, double lrb); 
    complex_t G8(double lrb);
    complex_t H2(double s, double a1, double a2);
    double H8(double a1, double a2); 
    complex_t a7c_h(double mu_h, double mu_b, double alpha_s_mu_h, double f_B, double f_Ks_perp, double T1, double m_B, double lambda_B, complex_t h2, complex_t h8); 
    complex_t a7c_b(double alpha_s_mu_b, complex_t g2, complex_t g8); 
    complex_t r1(double mu_0, double mu_b, double F_p);
    complex_t r2(double mu_0, double mu_b);
    complex_t K1(double mb_mb, double m_B, double alpha_s_mu_b, double F_p, complex_t G_p, complex_t X_p, complex_t r1, double mu_b);
    complex_t K2d(double mb_mb, double alpha_s_mu_b, complex_t H_p, complex_t r2, double mu_b);
    complex_t ckm_factor(complex_t Vus, complex_t Vub, complex_t Vcs, complex_t Vcb);
    complex_t K2u(complex_t ckm, complex_t K2d);
    double delta_0(double f_B, double mb_mb, double T1, double f_Ks_perp, double f_Ks_par, double m_Ks, double m_B, double lambda_B, complex_t a7c, complex_t K1, complex_t K2d, complex_t K2u);

private:
    const QCDOrder max_order = QCDOrder::NNLO;

public:
    BKstarDecay(QCDOrder order, double matching_scale, double hadronic_scale) : DecayParent(matching_scale, hadronic_scale, order) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime};
    }

    void enable() {
        DecayParent::enable();
        WilsonAdapter().switchbasis(WGroup::B);
    }

    void build_op_tree() override;

};

#endif // __BKSTARDECAY_H__
