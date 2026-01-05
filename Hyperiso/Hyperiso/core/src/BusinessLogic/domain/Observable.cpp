#include "Observable.h"

scalar_t Observable::get_exp_val() const {
    ObsParameterProxy opp {ParameterType::OBSERVABLE};
    return opp("FOBS", ObservableMapper::flha_of(this->id).value(), DataType::VALUE);
}

scalar_t Observable::get_exp_uncertainty(UncertaintyType u_type) const {
    ObsParameterProxy opp {ParameterType::OBSERVABLE};
    return opp("FOBS", ObservableMapper::flha_of(this->id).value(), UncertaintyTypeMapper::d_type(u_type));
}

std::vector<ObservableValue> Observable::compute() const {
    decay_parent->enable();
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
