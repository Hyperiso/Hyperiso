#ifndef BKNUNU_DECAY_H
#define BKNUNU_DECAY_H

#include <array>
#include <memory>

#include "BPFFCalculator.h"
#include "DecayParent.h"
#include "DefaultConfig.h"

struct BKnunuConfig : public DecayConfig {
    // Keep the default aligned with BKllConfig so B -> K X channels use the
    // same form-factor convention unless explicitly overridden.
    BP_FF_Src ff_src {BP_FF_Src::AS};
};

struct BKnunuDecayCache {
    double G_F{};
    double alpha_em{};
    complex_t lambda_t{};

    double m_Bp{};
    double m_B0{};
    double m_Kp{};
    double m_K0{};
    double tau_Bp{};
    double tau_B0{};

    double m_tau{};
    double f_K{};
    double f_B{};
    complex_t Vus{};
    complex_t Vub{};

    std::array<complex_t, 9> C_L{};
    std::array<complex_t, 9> C_R{};

    std::shared_ptr<BPFFCalculator> ff_charged;
    std::shared_ptr<BPFFCalculator> ff_neutral;
};

/**
 * @brief B -> K nu anti-nu total branching fractions.
 *
 * The B -> K form factors are evaluated through BPFFCalculator so their
 * configured nuisance parameters and correlations are propagated by the
 * statistical layer.
 */
class BKnunuDecay : public DecayParentConfigurable<BKnunuConfig> {
private:
    BKnunuConfig cfg {};
    BKnunuDecayCache cache;

    double loop_br(bool charged);
    double charged_tree_br();
    double coefficient_sum_plus() const;

public:
    BKnunuDecay(QCDOrder order, double matching_scale, double hadronic_scale,
                ObservablePortsConfig& ports)
        : DecayParentConfigurable(DecayMapper::to_id(Decays::B__K_nu_nu),
                                  matching_scale, hadronic_scale, order, ports)
    {
        this->w_config.groups = {GroupMapper::to_id(WGroup::BNuNu)};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    void set_config_spe(BKnunuConfig config) override { this->cfg = config; }
    std::any get_config() const override { return cfg; }
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
};

#endif
