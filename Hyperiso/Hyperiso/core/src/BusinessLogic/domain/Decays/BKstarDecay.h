#ifndef BKSTARDECAY_H
#define BKSTARDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "QCDHelper.h"
#include "Math.h"
#include "ObsParameterMutator.h"
#include "DefaultConfig.h"

struct BKstarDecayCache {
    double m_b_mu_b, m_b_1S;
    double f_Ks_par, f_Ks_perp, f_B;
    double m_B, m_Ks;
    double a_1_perp, a_2_perp;
    double a_1_par, a_2_par;
    double zeta_3_A, zeta_3_V;
    double omega_10_A;
    double delta_t_p, delta_t_m;
    double lambda_B;
    double T1_B_Ks;
    double Lambda_h, mu_0;

    complex_t lambda_hat_u;
    double C_F, Nc, beta_0, n_f;
    double mu_b;
    double alpha_s_mu_b, alpha_s_mu_h;
    double s_c;

    std::unordered_map<WCoef, complex_t> C_b;
    complex_t C2_h, C8_h;
};

class BKstarDecay : public DecayParentConfigurable<DecayConfig> {
private:
    BKstarDecayCache cache;

protected:
    complex_t h(double u);
    double phi_perp(double u);
    double gv_dga_4(double u);
    complex_t G(double xbar);
    complex_t G_perp(); 
    complex_t H_perp(); 
    double X_perp(); 
    complex_t G2(); 
    complex_t G8();
    complex_t H2();
    complex_t K1();
    complex_t K2(int q);
    double delta_0();

public:
    BKstarDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__Kstar), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime};
        this->max_order = QCDOrder::NNLO;
    }

    void enable() {
        DecayParent::enable();
        this->w_proxy->set_basis(WilsonBasis::B_TRADITIONAL);
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
};

#endif // __BKSTARDECAY_H__
