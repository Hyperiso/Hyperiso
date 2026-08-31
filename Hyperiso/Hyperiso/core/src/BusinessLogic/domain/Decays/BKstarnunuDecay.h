#ifndef BKSTARNUNU_DECAY_H
#define BKSTARNUNU_DECAY_H

#include <array>
#include <memory>
#include <utility>

#include "BVFFCalculator.h"
#include "DecayParent.h"
#include "DefaultConfig.h"

struct BKstarnunuDecayCache {
    double G_F{};
    double alpha_em{};
    complex_t lambda_t{};

    double m_Bp{};
    double m_B0{};
    double m_Kstarp{};
    double m_Kstar0{};
    double tau_Bp{};
    double tau_B0{};

    double m_tau{};
    double f_Kstar{};
    double f_B{};
    complex_t Vus{};
    complex_t Vub{};

    std::array<complex_t, 9> C_L{};
    std::array<complex_t, 9> C_R{};

    std::shared_ptr<BVFFCalculator> ff_charged;
    std::shared_ptr<BVFFCalculator> ff_neutral;
};

/** @brief B -> K* nu anti-nu total branching fractions. */
class BKstarnunuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    BKstarnunuDecayCache cache;

    double loop_br(bool charged);
    double charged_tree_br();
    std::pair<double, double> coefficient_sums() const;

public:
    BKstarnunuDecay(QCDOrder order, double matching_scale, double hadronic_scale,
                    ObservablePortsConfig& ports)
        : DecayParentConfigurable(DecayMapper::to_id(Decays::B__Kstar_nu_nu),
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
