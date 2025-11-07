#pragma once
#include <vector>
#include <stdexcept>
#include "ports/IModel.h"
#include "ObservableInterface.h" // your API


struct ParamSpec { std::string block; LhaID code; ParameterType type; };


class ObservableInterfaceAdapterObs final : public IModel {
public:
    ObservableInterfaceAdapterObs(ObservableInterface& oi,
    std::vector<Observables> obs_ids,
    std::vector<ParamSpec> p_specs,
    std::vector<ParamSpec> eta_specs)
    : oi_(oi), obs_ids_(std::move(obs_ids)), p_specs_(std::move(p_specs)), eta_specs_(std::move(eta_specs)) {}


    std::size_t n_observables() const override { return obs_ids_.size(); }


    Vec predict(const Vec& p, const Vec& eta) const override {
        if (p.size()!=p_specs_.size() || eta.size()!=eta_specs_.size())
        throw std::invalid_argument("(p,eta) sizes do not match specs");
        for (std::size_t i=0;i<p.size();++i) {
            const auto& s = p_specs_[i];
            oi_.set_param(s.block, s.code, p[i], static_cast<ParameterType>(s.type));
        }
        for (std::size_t i=0;i<eta.size();++i) {
            const auto& s = eta_specs_[i];
            oi_.set_param(s.block, s.code, eta[i], static_cast<ParameterType>(s.type));
        }
        auto all = oi_.compute_all();
        Vec out; out.reserve(obs_ids_.size());
        for (auto oid : obs_ids_) {
            out.push_back(oi_.compute_observable(oid).front().value); // or from compute_all()
        }
        return out;
    }
private:
    ObservableInterface& oi_;
    std::vector<Observables> obs_ids_;
    std::vector<ParamSpec> p_specs_, eta_specs_;
};