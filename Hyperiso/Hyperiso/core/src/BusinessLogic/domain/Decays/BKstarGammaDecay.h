#ifndef __BKSTARGAMMADECAY_H__
#define __BKSTARGAMMADECAY_H__

#include "DecayParent.h"
#include "Include.h"
#include "ObsQCDProxy.h"
#include "DefaultConfig.h"
#include "BVQCDfCalculator.h"
#include "BVFFCalculator.h"
#include "ObsParameterMutator.h"

struct BKstarGammaConfig : public DecayConfig {
    enum class B_Charge {B_0, B_PLUS};

    BV_FF_Src ff_src {BV_FF_Src::BSZ_SR_LAT};
    B_Charge charge {B_Charge::B_PLUS};
};

struct BKstarGammaCache {
    std::map<WCoef, complex_t> C;
    std::map<WCoef, complex_t> C_trad;
    complex_t C2_h, C8_h;

    BVFFCalculator ff_calculator;
    BVQCDfCalculator qcdf_calculator;

    double alpha_em;
    double m_b_m_b, m_b_mu_b;
    double z;
    double f_Ks_par, f_Ks_perp, f_B;
    double lambda_B;
    double m_B, m_Ks;
    double tau_B;
    complex_t N_prime;
    double mu_b, mu_0, mu_h, L_b;
    complex_t lambda_hat_u;
    double C_F, Nc, n_f;
    double alpha_s_mu_b, alpha_s_mu_h;
};

/**
 * @brief Decay parent for the exclusive B > K* l+ l- decays. Implements the integrated branching ratio and angular observables in several q² bins, as well as the forward-backward asymmetry of the decay. 
 */
class BKstarGammaDecay : public DecayParentConfigurable<BKstarGammaConfig> {
public:
    BKstarGammaDecay(QCDOrder order, double matching_scale, double hadronic_scale, ObservablePortsConfig& ports) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__Kstar_l_l), matching_scale, hadronic_scale, order, ports) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime), GroupMapper::to_id(WGroup::BScalar)};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

    void set_config_spe(BKstarGammaConfig config) override {this->cfg = config;}

private:
    BKstarGammaConfig cfg {};
    BKstarGammaCache cache;

protected:
    // Auxiliary
    void fill_wilson_cache();
    void load_cfg_dependent_params();
    void set_cfg_flags(BKstarGammaConfig::B_Charge charge);

    // Helicity amplitudes
    complex_t H_V(double sign, bool bar);
    complex_t K1();
    complex_t K2(int q);
    
    // Observables
    double BR();
    double delta_0();
};


#endif // __BKSTARGAMMADECAY_H__
