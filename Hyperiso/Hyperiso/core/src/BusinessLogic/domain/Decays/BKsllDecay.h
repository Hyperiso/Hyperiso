#ifndef BKSTARLLDECAY_H
#define BKSTARLLDECAY_H

#include <thread>
#include <algorithm>

#include "DecayParent.h"
#include "Include.h"
#include "ObsQCDProxy.h"
#include "DefaultConfig.h"
#include "BVQCDfCalculator.h"
#include "BVFFCalculator.h"

struct BKstarllConfig : public DecayConfig {
    enum class Power_Corrections_Impl {BFS, BCvDV, KMPW};
    enum class B_Charge {B_0, B_PLUS};
    enum class Lepton {E, MU, TAU};

    BV_FF_Src ff_src {BV_FF_Src::BSZ_SR_LAT};
    B_FF_Type ff_type {B_FF_Type::FULL};
    Power_Corrections_Impl power_corr_impl {Power_Corrections_Impl::BFS};
    B_Charge charge {B_Charge::B_PLUS};
    Lepton gen {Lepton::MU};

    size_t n_threads {1};
};

struct BKstarllCache {
    std::map<WCoef, complex_t> C;
    BVFFCalculator ff_calculator;
    BVQCDfCalculator qcdf_calculator;

    double G_F, alpha_em;
    double mu_b, L_b;
    double m_l, m_s, m_c_mu_b, m_b_PS, m_b_mu_b;
    double Delta_M;
    double alpha_s_mu_b;
    double m_B, m_Ks;
    double life_B;
    complex_t lambda_hat_u;
    double kappa;
    complex_t N_0;
    double q2_min, q2_max;
    double q2_low, q2_high;
    
    double q2_lookup_min;

    // Power corrections, hadronic fit
    std::array<complex_t, 3> h_p_fit {};
    std::array<complex_t, 3> h_m_fit {}; 
    std::array<complex_t, 3> h_0_fit {};

    // Power corrections, guesstimate
    std::array<complex_t, 6> a_k_low {};
    std::array<complex_t, 6> b_k_low {};
    std::array<complex_t, 8> a_k_high {};
    std::array<complex_t, 6> phi_k_low {};
    std::array<complex_t, 6> theta_k_low {};
    std::array<complex_t, 8> phi_k_high {};

    // NF Corrections, Van Dyk
    double tp_nf, t0_nf;
    double z0_nf, z_Jpsi_nf, z_psi2S_nf;
    std::array<complex_t, 3> alpha_perp;
    std::array<complex_t, 3> alpha_par;
    std::array<complex_t, 2> alpha_0;

    // NF Corrections, KMPW
    double q2_bar, q2_Jpsi; // q2bar = 1.0 GeV², q2_Jpsi = m_Jpsi^2
    std::array<double, 3> DeltaC9_M_qbar;
    std::array<double, 3> r1_M;
    std::array<double, 3> r2_M;

    // Lookups
    static inline constexpr size_t LOOKUP_SIZE = 50;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_p_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_m_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_par_m_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_p_bar_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_perp_m_bar_lookup;
    std::array<scalar_t, LOOKUP_SIZE> T_par_m_bar_lookup;

    std::array<std::vector<double>, 15> J_i_binned;
    std::array<std::vector<double>, 15> J_i_bar_binned; 
};

/**
 * @brief Decay parent for the exclusive B > K* l+ l- decays. Implements the integrated branching ratio and angular observables in several q² bins, as well as the forward-backward asymmetry of the decay. 
 */
class BKstarllDecay : public DecayParentConfigurable<BKstarllConfig> {
public:
    BKstarllDecay(QCDOrder order, double matching_scale, double hadronic_scale, ObservablePortsConfig& ports) : DecayParentConfigurable(DecayMapper::to_id(Decays::B__Kstar_l_l), matching_scale, hadronic_scale, order, ports) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime), GroupMapper::to_id(WGroup::BScalar)};
        this->max_order = QCDOrder::NNLO;
        this->binned = true;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

    void set_config_spe(BKstarllConfig config) override {this->cfg = config;}
    void set_n_threads(size_t n_threads) override;

private:
    BKstarllConfig cfg {};
    BKstarllCache cache;

protected:
    // Auxiliary
    void fill_wilson_cache();
    void load_cfg_dependent_params();
    void set_lepton_gen_and_charge(BKstarllConfig::Lepton gen, BKstarllConfig::B_Charge charge);

    // QCDf
    complex_t T_perp_p_cached(double q2, bool bar);
    complex_t T_perp_m_cached(double q2, bool bar);
    complex_t T_par_m_cached(double q2, bool bar);

    // Nonfactorisable corrections
    double beta_l(double q2);
    double lambda(double q2);
    complex_t N(double q2, bool bar);

    complex_t delta_A_perp_QCDf(double q2, double sign, bool bar);
    complex_t delta_A_perp_vD(double q2, bool bar);
    complex_t delta_A_perp_K(double q2, bool bar);
    complex_t delta_A_perp(double q2, double sign, bool bar);

    complex_t delta_A_par_QCDf(double q2, double sign, bool bar);
    complex_t delta_A_par_vD(double q2, bool bar);
    complex_t delta_A_par_K(double q2, bool bar);
    complex_t delta_A_par(double q2, double sign, bool bar);

    complex_t delta_A_0_QCDf(double q2, double sign, bool bar);
    complex_t delta_A_0_vD(double q2, bool bar);
    complex_t delta_A_0_K(double q2, bool bar);
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

    void compute_binned_J_i();

    std::vector<ObservableValue> dBR_dq2_binned(bool bar, Observables id, bool br=true);
    double dG_dq2_avg_bin(size_t bin);
    std::vector<ObservableValue> A_FB_binned(Observables id, bool cpv);
    ObservableValue q0(Observables id);
    std::vector<ObservableValue> A_CP_binned(Observables id);
    std::vector<ObservableValue> F_L_binned(Observables id);
    std::vector<ObservableValue> F_T_binned(Observables id);
    std::vector<ObservableValue> A_T_1_binned(Observables id);
    std::vector<ObservableValue> A_T_2_binned(Observables id);
    std::vector<ObservableValue> A_T_3_binned(Observables id);
    std::vector<ObservableValue> A_T_4_binned(Observables id);
    std::vector<ObservableValue> A_T_5_binned(Observables id);
    std::vector<ObservableValue> A_T_Re_binned(Observables id);
    std::vector<ObservableValue> A_T_Re_CPV_binned(Observables id);
    std::vector<ObservableValue> A_Im_binned(Observables id);
    std::vector<ObservableValue> alpha_K_binned(Observables id);
    std::vector<ObservableValue> H_T_1_binned(Observables id);
    std::vector<ObservableValue> H_T_2_binned(Observables id);
    std::vector<ObservableValue> H_T_3_binned(Observables id);
    std::vector<ObservableValue> P_2_binned(Observables id);
    std::vector<ObservableValue> P_3_binned(Observables id);
    std::vector<ObservableValue> P_6_binned(Observables id);
    std::vector<ObservableValue> P_8_binned(Observables id);
    std::vector<ObservableValue> Pp_i_binned(size_t i, bool cpv, Observables id);
    std::vector<ObservableValue> S_i_binned(size_t i, bool cpv, Observables id);
    std::vector<ObservableValue> P_i_CPV_binned(size_t i, Observables id);
    std::vector<ObservableValue> Rm1_BKstar(Observables id, BKstarllConfig::B_Charge charge);

    // Tests
    void test_ff();
    void test_T();
    void test_J();
    void test_binned_obs();
};

#endif // BKSTARLLDECAY_H