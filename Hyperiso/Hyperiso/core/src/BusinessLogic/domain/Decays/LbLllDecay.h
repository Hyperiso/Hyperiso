#ifndef LBLLLDECAY_H
#define LBLLLDECAY_H

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"
#include "ObsQCDProxy.h"
#include "LbLFFCalculator.h"

struct LbLllConfig {
    enum class Lepton {E, MU, TAU};

    LbL_FF_Src ff_src {LbL_FF_Src::DM};
    Lepton gen {Lepton::MU};
};

struct LbLllDecayCache {
    LbLllDecayCache(std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_sm) : ff_calculator(iobspp_sm) {}
    std::map<WCoef, complex_t> C;
    LbLFFCalculator ff_calculator;

    double m_l, m_b_mu_b;
    double m_Lb, m_L;
    double life_L;
    complex_t N_0;
    double alpha_L;
    double q2_min, q2_max;

    std::array<std::vector<double>, 6> K_i_binned;
    std::array<std::vector<double>, 6> K_i_bar_binned;
};

/**
 * @brief Decay parent for the Lambda_b > Lambda l+ l- decays. 
 */
class LbLllDecay : public DecayParentConfigurable<LbLllConfig> {
private:
    LbLllDecayCache cache;
    LbLllConfig cfg {};

protected:
    // Auxiliary
    void fill_wilson_cache();
    void set_cfg_flags(LbLllConfig::Lepton gen);
    void load_cfg_dep_params();

    // Kinematics
    double beta_l(double q2);
    double lambda(double q2);
    double s_p(double q2);
    double s_m(double q2);
    complex_t N(double q2, bool bar);

    // Transversity amplitudes
    complex_t A_perp_1(double q2, double sign, bool bar);
    complex_t A_par_1(double q2, double sign, bool bar);
    complex_t A_perp_0(double q2, double sign, bool bar);
    complex_t A_par_0(double q2, double sign, bool bar);

    // Angular coefficients
    double K1ss(double q2, bool bar);
    double K1cc(double q2, bool bar);
    double K1c(double q2, bool bar);
    double K2ss(double q2, bool bar);
    double K2cc(double q2, bool bar);
    double K2c(double q2, bool bar);

    void compute_binned_K_i();

    std::vector<ObservableValue> dBR_dq2_binned(Observables oid);
    double dG_dq2_avg_bin(size_t bin);
    std::vector<ObservableValue> A_FB_l(Observables oid);
    std::vector<ObservableValue> A_FB_h(Observables oid);
    std::vector<ObservableValue> A_FB_lh(Observables oid);
    std::vector<ObservableValue> F_L(Observables oid);
    std::vector<ObservableValue> F_T(Observables oid);

public:
    LbLllDecay(QCDOrder order, double matching_scale, double hadronic_scale, ObservablePortsConfig& ports) :  DecayParentConfigurable(DecayMapper::to_id(Decays::Lambda_b__Lambda_l_l), matching_scale, hadronic_scale, order, ports), cache(ports.iobspp_sm) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime)};
        this->max_order = QCDOrder::NNLO;
        this->binned = true;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;

    void set_config_spe(LbLllConfig config) override {this->cfg = config;}
};

#endif // __LBLLLDECAY_H__
