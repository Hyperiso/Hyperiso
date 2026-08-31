#ifndef BKNUNU_DECAY_H
#define BKNUNU_DECAY_H

#include <array>

#include "DecayParent.h"
#include "DefaultConfig.h"

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
};

/**
 * @brief B -> K nu anti-nu total branching fractions.
 *
 * Implements the WET normalization and central B->K form-factor fit of
 * arXiv:2301.06990, including the charged-mode on-shell-tau weak-annihilation
 * contribution.  Form-factor nuisance parameters are deliberately postponed
 * to the follow-up assets/bindings patch; this C++ core patch uses the paper's
 * central fit values.
 */
class BKnunuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    BKnunuDecayCache cache;

    double loop_br(bool charged);
    double charged_tree_br();
    double f_plus(double q2, double m_B, double m_K) const;
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
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
};

#endif
