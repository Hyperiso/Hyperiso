#ifndef __KPINUNUDECAY_H__
#define __KPINUNUDECAY_H__

#include "DecayParent.h"
#include "General.h"
#include "DefaultConfig.h"

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

    complex_t CL;
};

/**
 * @brief Decay parent for the K_L,S > ll decays. Currently implements BR(K_L,S > mu mu) with both possible signs for the K_L > 2 gamma long distance correction.
 */
class KPinunuDecay : public DecayParentConfigurable<DecayConfig> {
private:
    KPinunuDecayCache cache;

protected:  
    double P_c();
    double BR_L();
    double BR_p();

public:
    KPinunuDecay(QCDOrder order, double matching_scale, double hadronic_scale, std::shared_ptr<ObsWilsonBuilder>& wilson_builder) : DecayParentConfigurable(DecayMapper::to_id(Decays::K__pi_nu_nu), matching_scale, hadronic_scale, order, wilson_builder) {
        this->w_config.groups = {GroupMapper::to_id(WGroup::K)};
        this->max_order = QCDOrder::NNLO;
    }

    void load_params() override;
    std::vector<ObservableValue> compute_observable(Observables obs) override;
    std::vector<ObservableValue> compute_observable(ObservableId obs) override;
};

#endif // __KPINUNUDECAY_H__
