#ifndef __KPINUNUDECAY_H__
#define __KPINUNUDECAY_H__

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"

#include <array>

struct KPinunuDecayCache {
    double alpha_s_m_Z;
    double sw2;
    double m_c_m_c;
    complex_t lambda_c;
    complex_t lambda_t;
    double lambda;

    double kappa_L;
    double kappa_p;
    double delta_em;

    std::array<complex_t, 3> C_L{};
    std::array<complex_t, 3> C_R{};
};

/**
 * @brief Rare K -> pi nu anti-nu decay parent.
 *
 * The short-distance top contribution is flavour resolved in the KNuNu Wilson
 * group.  The charm contribution P_c is SM-like and diagonal/universal in the
 * current implementation.
 */
class KPinunuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    KPinunuDecayCache cache;

protected:
    double P_c();
    double BR_L();
    double BR_p();

public:
    KPinunuDecay(QCDOrder order, double matching_scale, double hadronic_scale, ObservablePortsConfig& ports)
        : DecayParentConfigurable(DecayMapper::to_id(Decays::K__pi_nu_nu), matching_scale, hadronic_scale, order, ports) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::KNuNu)};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
};

#endif // __KPINUNUDECAY_H__
