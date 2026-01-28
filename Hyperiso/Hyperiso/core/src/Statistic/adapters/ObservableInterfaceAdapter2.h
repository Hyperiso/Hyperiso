#pragma once
#include <vector>
#include <stdexcept>
#include "ports/IModel.h"
#include "ObservableInterface.h"
#include "StatParamOptimizerProxy.h"


class ObservableInterfaceAdapterObs final : public IModel {
public:
    ObservableInterfaceAdapterObs(
    std::shared_ptr<ObservableInterface> obs,
    std::map<ObservableId, QCDOrder> obs_ids,
    std::vector<ParamId> p_specs,
    std::vector<ParamId> eta_specs)
    : p_specs_(std::move(p_specs)), eta_specs_(std::move(eta_specs)) {
        oi_ = obs;

        oi_->add_observables(obs_ids, true);

        for (auto elem : obs_ids) {
            obs_ids_.push_back(elem.first);
        }
    }

    ObservableInterfaceAdapterObs(std::shared_ptr<ObservableInterface> obs) { oi_ = obs; }
    std::size_t n_observables() const override { return obs_ids_.size(); }

    void add_observables(std::map<ObservableId, QCDOrder> obs_ids) override {
        oi_->add_observables(obs_ids, true);

        for (auto elem : obs_ids) {
            obs_ids_.push_back(elem.first);
        }
    }

    std::unordered_set<ParamId> get_obs_deps(ObservableId id) override {
        return oi_->get_all_ops_deps(id);
    }

    Vec predict(const Vec& p, const Vec& eta) override {
        auto obs = oi_->get_current_observables();
        oi_->enable_obs();

        if (p.size()!=p_specs_.size() || eta.size()!=eta_specs_.size())
        throw std::invalid_argument("(p,eta) sizes do not match specs");
        for (std::size_t i=0;i<p.size();++i) {
            const auto& s = p_specs_[i];
            oi_->set_param(s.block, s.code, p[i], s.type.value_or(ParameterType::SM)); //TODO check value_or
        }
        for (std::size_t i=0;i<eta.size();++i) {
            const auto& s = eta_specs_[i];
            oi_->set_param(s.block, s.code, eta[i], s.type.value_or(ParameterType::SM)); //TODO check value_or
        }
        // auto all = oi_.compute_all();
        Vec out; out.reserve(obs_ids_.size());
        for (auto oid : obs_ids_) {
            out.push_back(oi_->compute_observable(oid).front().value); // or from compute_all()
        }
        return out;
    }

    std::map<ObservableId, double> predict(const std::map<ParamId, double>& p, const std::map<ParamId, double>& eta) override {
        auto obs = oi_->get_current_observables();
        oi_->enable_obs();

        // if (p.size()!=p_specs_.size() || eta.size()!=eta_specs_.size())
        // throw std::invalid_argument("(p,eta) sizes do not match specs");
        for (auto p_elem : p) {
            const auto& s = p_elem.first;
            oi_->set_param(s.block, s.code, p_elem.second, s.type.value_or(ParameterType::SM)); //TODO check value_or
        }
        for (auto eta_elem : eta) {
            const auto& s = eta_elem.first;
            oi_->set_param(s.block, s.code, eta_elem.second, s.type.value_or(ParameterType::SM)); //TODO check value_or
        }
        // auto all = oi_.compute_all();
        std::map<ObservableId, double> out;
        for (auto oid : obs_ids_) {
            out[oid] = oi_->compute_observable(oid).front().value; // or from compute_all()
        }
        return out;
    }

    std::map<ObservableId, double> predict_optimized(const std::map<ParamId, double>& p, const std::map<ParamId, double>& eta) override {
        StatParamOptimizerProxy spop = StatParamOptimizerProxy();

        auto obs = oi_->get_current_observables();
        oi_->enable_obs();

        // if (p.size()!=p_specs_.size() || eta.size()!=eta_specs_.size())
        // throw std::invalid_argument("(p,eta) sizes do not match specs");
        for (auto p_elem : p) {
            const auto& s = p_elem.first;
            spop.set_value(s.block, s.code, p_elem.second);
            // oi_->set_param(s.block, s.code, p_elem.second, s.type.value_or(ParameterType::SM)); //TODO check value_or
        }
        for (auto eta_elem : eta) {
            const auto& s = eta_elem.first;
            spop.set_value(s.block, s.code, eta_elem.second);
            // oi_->set_param(s.block, s.code, eta_elem.second, s.type.value_or(ParameterType::SM)); //TODO check value_or
        }
        spop.commit();
        auto v1 = oi_->compute_observable(obs_ids_[0]).front().value;
        auto v2 = oi_->compute_observable(obs_ids_[0]).front().value;


        double a = oi_->compute_observable(obs_ids_[0]).front().value;
        double b = oi_->compute_observable(obs_ids_[1]).front().value;
        // recompute le premier après avoir calculé le second
        double a2 = oi_->compute_observable(obs_ids_[0]).front().value;

        // auto all = oi_.compute_all();
        std::map<ObservableId, double> out;
        for (auto oid : obs_ids_) {
            out[oid] = oi_->compute_observable(oid).front().value; // or from compute_all()
        }
        return out;
    }
private:
    std::shared_ptr<ObservableInterface> oi_;
    std::vector<ObservableId> obs_ids_;
    std::vector<ParamId> p_specs_, eta_specs_;
};