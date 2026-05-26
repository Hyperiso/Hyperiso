#include "ObservableInterfaceProxy.h"

ObservableInterfaceProxy::ObservableInterfaceProxy(
    std::shared_ptr<ObservableInterface> obs,
    std::vector<ParamId> p_specs,
    std::vector<ParamId> eta_specs)
    : p_specs_(std::move(p_specs)), eta_specs_(std::move(eta_specs)) {
    oi_ = obs;
}

ObservableInterfaceProxy::ObservableInterfaceProxy(std::shared_ptr<ObservableInterface> obs, std::shared_ptr<IStatParamOptimizerProxy> spop) { oi_ = obs; spop_ = spop;}
std::size_t ObservableInterfaceProxy::n_observables() const { return oi_->get_current_observables().size(); }

std::unordered_set<ParamId> ObservableInterfaceProxy::get_obs_deps(ObservableId id) {
    return oi_->get_all_ops_deps(id);
}

std::vector<BinnedObservableId> ObservableInterfaceProxy::get_obs_ids() {
    return oi_->get_current_observables();
}

std::map<ObservableId, std::vector<ObservableValue>> ObservableInterfaceProxy::predict_optimized(
    const std::map<ParamId, double>& p,
    const std::map<ParamId, double>& eta)
{
    bool has_nonzero_fit_param = false;
    double max_abs_p = 0.0;
    for (const auto& [pid, val] : p) {
        max_abs_p = std::max(max_abs_p, std::abs(val));
        if (std::abs(val) > 1e-12) {
            has_nonzero_fit_param = true;
        }
    }

    for (auto p_elem : p) {
        const auto& s = p_elem.first;
        spop_->set_value(s.block, s.code, p_elem.second);
    }
    for (auto eta_elem : eta) {
        const auto& s = eta_elem.first;
        spop_->set_value(s.block, s.code, eta_elem.second);
    }
    spop_->commit();

    auto pred = oi_->compute_all();

    return pred;
}

void ObservableInterfaceProxy::compute_observables() const {
    oi_->compute_all();
};
