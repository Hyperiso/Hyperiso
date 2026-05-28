#include "Observable.h"

scalar_t Observable::get_exp_val(std::pair<double, double> bins, std::string exp) const {
    return (*iobspp_obs)(ParamId(ParameterType::OBSERVABLE, "FOBS_" + exp, BinnedObservableId(this->id, bins).flha()), DataType::VALUE);
}

scalar_t Observable::get_exp_uncertainty(std::pair<double, double> bins, std::string exp, UncertaintyType u_type) const {
    return (*iobspp_obs)(ParamId(ParameterType::OBSERVABLE, "FOBS_" + exp, BinnedObservableId(this->id, bins).flha()), UncertaintyTypeMapper::d_type(u_type));
}

std::vector<ObservableValue> Observable::compute() const {
    return decay_parent->compute_observable(id);
}

void Observable::add_dependence(const ParamId &param_name) {
    dependences.emplace(param_name);
}

void Observable::add_dependences(const std::unordered_set<ParamId> &param_names) {
    LOG_DEBUG("Adding parameter list to compound");
    std::set_union(dependences.begin(), dependences.end(),
                   param_names.begin(), param_names.end(),
                   std::inserter(dependences, dependences.begin()));
}

const std::unordered_set<ParamId> &Observable::get_dependences() const {
    return dependences;
}
