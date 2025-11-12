#ifndef __BKLLDECAY_H__
#define __BKLLDECAY_H__

#include "DecayParent.h"
#include "Include.h"
#include "ObsQCDProxy.h"
#include "DefaultConfig.h"
#include "BPQCDfCalculator.h"
#include "BPFFCalculator.h"

struct BKllConfig : public DecayConfig {
    enum class B_Charge {B_0, B_PLUS};
    enum class Lepton {E, MU, TAU};

    BP_FF_Src ff_src {BP_FF_Src::AS};
    B_FF_Type ff_type {B_FF_Type::FULL};
    B_Charge charge {B_Charge::B_PLUS};
    Lepton gen {Lepton::MU};

    std::vector<std::pair<double, double>> bins {
        { 1.00,  4.00},
        {15.00, 18.00}
    };
};

struct BKllCache {
    std::map<WCoef, complex_t> C;
    BPFFCalculator ff_calculator;
    BPQCDfCalculator qcdf_calculator;

    double G_F, alpha_em;
    double mu_b, L_b;
    double m_l, m_s, m_c_mu_b, m_b_PS, m_b_mu_b;
    double alpha_s_mu_b;
    double m_B, m_K;
    double Delta_M;
    complex_t lambda_hat_u;
    double N_0;
    double q2_min, q2_max;
    double q2_low, q2_high;

    // Power corrections, guesstimate
    std::array<double, 4> A_had_err_low_0;
    std::array<double, 4> A_had_err_low_1;
    std::array<double, 4> A_had_err_high;

    // Lookups
    static inline constexpr size_t LOOKUP_SIZE = 50;
    std::array<scalar_t, LOOKUP_SIZE> T_P_lookup;
    std::array<std::vector<double>, 3> abc_binned; 
};

/**
 * @brief Decay parent for the exclusive B > K l+ l- decays. Implements the integrated branching ratio, forward-backward asymmetry and flat term in several q² bins. 
 */
class BKllDecay : public DecayParentConfigurable<BKllConfig> {
public:
    BKllDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__K_l_l), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime), GroupMapper::to_id(WGroup::BScalar)};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

    void set_config_spe(BKllConfig config) override {this->cfg = config;}

private:
    BKllConfig cfg {};
    BKllCache cache;

protected:
    // Auxiliary
    void fill_wilson_cache();

    // QCDf
    complex_t T_P_cached(double q2);

    // Nonfactorisable corrections
    double beta_l(double q2);
    double lambda(double q2);
    double N(double q2);

    // Helicity amplitudes
    complex_t F_V_low(double q2);
    complex_t F_A_low(double q2);
    complex_t F_P_low(double q2);
    complex_t F_S_low(double q2);

    complex_t C7_eff(double q2);
    complex_t C9_eff(double q2);
    complex_t F_V_high(double q2);
    complex_t F_A_high(double q2);
    complex_t F_P_high(double q2);
    complex_t F_S_high(double q2);

    complex_t interpolate(double q2, complex_t val_low, complex_t val_high);
    complex_t F_V(double q2);
    complex_t F_A(double q2);
    complex_t F_P(double q2);
    complex_t F_S(double q2);
    
    // Angular coefficients
    double a(double q2);
    double b(double q2);
    double c(double q2);

    void compute_binned_abc();

    // Observables
    std::vector<ObservableValue> dG_dq2();
    std::vector<ObservableValue> A_FB();
    std::vector<ObservableValue> F_H();
};


#endif // __BKLLDECAY_H__
