#ifndef __BSPHIDECAY_H__
#define __BSPHIDECAY_H__

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"
#include "ObsQCDProxy.h"
#include "BVFFCalculator.h"
#include "BVQCDfCalculator.h"

struct BsPhiConfig {
    enum class FF_Type {FULL, SOFT};
    enum class Lepton {E, MU, TAU};

    BV_FF_Src ff_src {BV_FF_Src::BSZ_SR_LAT};
    B_FF_Type ff_type {B_FF_Type::FULL};
    Lepton gen {Lepton::MU};

    std::vector<std::pair<double, double>> bins {
        { 1.00,  6.00},
        {16.00, 18.00}
    };
};

struct BsPhiDecayCache {
    std::map<WCoef, complex_t> C;
    BVFFCalculator ff_calculator;
    BVQCDfCalculator qcdf_calculator;

    double G_F, alpha_em;
    double mu_b, L_b;
    double m_l, m_s, m_c_mu_b, m_b_PS, m_b_mu_b;
    double alpha_s_mu_b;
    double m_Bs, m_phi;
    complex_t lambda_hat_u;
    double kappa;
    complex_t N_0;
    double q2_min, q2_max;
    double q2_low, q2_high;
    double ys, phi_s;
    complex_t um, up;

    // Power corrections, guesstimate
    std::array<double, 6> A_had_err_low_0;
    std::array<double, 6> A_had_err_low_1;
    std::array<double, 8> A_had_err_high;

    // Lookups
    static inline constexpr size_t LOOKUP_SIZE = 50;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_p_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_m_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_par_m_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_p_bar_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_m_bar_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_par_m_bar_lookup;

    std::array<std::vector<double>, 15> f_J_i_binned;
};

/**
 * @brief Decay parent for the Bs > phi l+ l- decays. 
 */
class BsPhiDecay : public DecayParentConfigurable<BsPhiConfig> {
private:
    BsPhiDecayCache cache {};
    BsPhiConfig cfg {};

protected:
    // Auxiliary
    void fill_wilson_cache();

    // QCDf 
    complex_t T_perp_p_cached(double q2, bool bar);
    complex_t T_perp_m_cached(double q2, bool bar);
    complex_t T_par_m_cached(double q2, bool bar);

    // Nonfactorisable corrections
    double beta_l(double q2);
    double lambda(double q2);
    complex_t N(double q2, bool bar);

    complex_t delta_A_perp(double q2, double sign, bool bar);
    complex_t delta_A_par(double q2, double sign, bool bar);
    complex_t delta_A_0(double q2, double sign, bool bar);

    // Transversity amplitudes
    complex_t A_perp_low(double q2, double sign, bool bar);
    complex_t A_par_low(double q2, double sign, bool bar);
    complex_t A_0_low(double q2, double sign, bool bar);
    complex_t A_t_low(double q2, bool bar);
    complex_t A_S_low(double q2, bool bar);

    complex_t C7_eff(double q2, bool bar);
    complex_t C9_eff(double q2, bool bar);
    complex_t A_perp_high(double q2, double sign, bool bar);
    complex_t A_par_high(double q2, double sign, bool bar);
    complex_t A_0_high(double q2, double sign, bool bar);
    complex_t A_t_high(double q2, bool bar);
    complex_t A_S_high(double q2, bool bar);

    complex_t interpolate(double q2, complex_t val_low, complex_t val_high);
    complex_t A_perp(double q2, double sign, bool bar);
    complex_t A_par(double q2, double sign, bool bar);
    complex_t A_0(double q2, double sign, bool bar);
    complex_t A_t(double q2, bool bar);
    complex_t A_S(double q2, bool bar);
    
    // Angular coefficients
    double J1s(double q2, bool bar);
    double J1c(double q2, bool bar);
    double J2s(double q2, bool bar);
    double J2c(double q2, bool bar);
    double J3(double q2, bool bar);
    double J4(double q2, bool bar);
    double J5(double q2, bool bar);
    double J6s(double q2, bool bar);
    double J6c(double q2, bool bar);
    double J7(double q2, bool bar);
    double J8(double q2, bool bar);
    double J9(double q2, bool bar);

    double h1s(double q2);
    double h1c(double q2);
    double h2s(double q2);
    double h2c(double q2);
    double h3(double q2);
    double h4(double q2);
    double h5(double q2);
    double h6s(double q2);
    double h6c(double q2);
    double h7(double q2);
    double h8(double q2);
    double h9(double q2);

    double s8(double q2);
    double s9(double q2);

    void compute_binned_J_i();

    std::vector<ObservableValue> dG_dq2_binned(bool bar);
    double dG_dq2_avg_bin(size_t bin);
    std::vector<ObservableValue> F_L();
    std::vector<ObservableValue> A_T_2();
    std::vector<ObservableValue> A_T_Re_CPV();
    std::vector<ObservableValue> A_T_Im_CPV();
    std::vector<ObservableValue> Pp_4();
    std::vector<ObservableValue> Pp_6();
    std::vector<ObservableValue> S_i(int i);
    std::vector<ObservableValue> A_i(int i);
    std::vector<ObservableValue> A_FB_CPV();
    std::vector<ObservableValue> P_2_CPV();
    std::vector<ObservableValue> P_3_CPV();
    std::vector<ObservableValue> Pp_5_CPV();
    std::vector<ObservableValue> Pp_8_CPV();
    std::vector<ObservableValue> Q_8_m();
    std::vector<ObservableValue> Q_8_p();
    std::vector<ObservableValue> Q_9();

public:
    BsPhiDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__l_nu), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {WGroup::B, WGroup::BPrime, WGroup::BScalar};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

    void set_config_spe(BsPhiConfig config) override {this->cfg = config;}
};


#endif // __BSPHIDECAY_H__
