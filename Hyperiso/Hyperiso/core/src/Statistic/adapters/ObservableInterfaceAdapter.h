#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include "ports/IModel.h"
#include "ObservableInterface.h" // provided by you


// Adapter to your ObservableInterface. It maps (p, η) vectors to your internal parameters,
// calls evaluate_all(), and returns the observables in a fixed order.


// struct ParamSpec {
// std::string block; int code; int type; // ParameterType as int to avoid header deps here
// };


class ObservableInterfaceAdapter final : public IModel {
public:
// obs_ids : the ordered list of ObservableId to report
// p_specs : how to set each component of p into Parameters
// eta_specs : how to set each component of eta into Parameters
ObservableInterfaceAdapter(ObservableInterface& oi,
std::vector<ObservableId> obs_ids,
std::vector<ParamId> p_specs,
std::vector<ParamId> eta_specs)
: oi_(oi), obs_ids_(std::move(obs_ids)), p_specs_(std::move(p_specs)), eta_specs_(std::move(eta_specs)) {}


std::size_t n_observables() const override { return obs_ids_.size(); }


Vec predict(const Vec& p, const Vec& eta) override {
if (p.size()!=p_specs_.size() || eta.size()!=eta_specs_.size())
throw std::invalid_argument("(p,eta) vector sizes do not match specs");


// set p
for (std::size_t i=0;i<p.size();++i) {
const auto& s = p_specs_[i];
oi_.set_param(s.block, s.code, p[i], s.type.value_or(ParameterType::SM)); //TODO check value_or
}
// set eta
for (std::size_t i=0;i<eta.size();++i) {
const auto& s = eta_specs_[i];
oi_.set_param(s.block, s.code, eta[i], s.type.value_or(ParameterType::SM)); //TODO check value_or
}


// compute all current observables
auto all = oi_.compute_all();
Vec out; out.reserve(obs_ids_.size());
for (ObservableId oid : obs_ids_) {
auto it = all.find(oid);
if (it==all.end()) throw std::runtime_error("ObservableId missing in compute_all()");
out.push_back(it->second.central_value);
}
return out;
}


private:
ObservableInterface& oi_;
std::vector<ObservableId> obs_ids_;
std::vector<ParamId> p_specs_;
std::vector<ParamId> eta_specs_;
};