#pragma once
#include <vector>
#include <stdexcept>
#include "ports/IModel.h"
#include "ObservableInterface.h"
#include "StatParamOptimizerProxy.h"
#include "StatParameterProxy.h"

class ObservableInterfaceAdapterObs final : public IModel {
public:
    ObservableInterfaceAdapterObs(
    std::shared_ptr<ObservableInterface> obs,
    std::vector<ParamId> p_specs,
    std::vector<ParamId> eta_specs)
    : p_specs_(std::move(p_specs)), eta_specs_(std::move(eta_specs)) {
        oi_ = obs;
    }

    ObservableInterfaceAdapterObs(std::shared_ptr<ObservableInterface> obs) { oi_ = obs; }
    std::size_t n_observables() const override { return oi_->get_current_observables().size(); }

    // void add_observables(std::map<ObservableId, QCDOrder> obs_ids) override {
    //     oi_->add_observables(obs_ids, true);

    //     for (auto elem : obs_ids) {
    //         obs_ids_.push_back(elem.first);
    //     }
    // }

    std::unordered_set<ParamId> get_obs_deps(ObservableId id) override {
        return oi_->get_all_ops_deps(id);
    }

    std::vector<BinnedObservableId> get_obs_ids() override {
        return oi_->get_current_observables();
    }

    // Vec predict(const Vec& p, const Vec& eta) override {
    //     auto obs = oi_->get_current_observables();
    //     oi_->enable_obs();

    //     if (p.size()!=p_specs_.size() || eta.size()!=eta_specs_.size())
    //     throw std::invalid_argument("(p,eta) sizes do not match specs");
    //     for (std::size_t i=0;i<p.size();++i) {
    //         const auto& s = p_specs_[i];
    //         oi_->set_param(s.block, s.code, p[i], s.type.value_or(ParameterType::SM)); //TODO check value_or
    //     }
    //     for (std::size_t i=0;i<eta.size();++i) {
    //         const auto& s = eta_specs_[i];
    //         oi_->set_param(s.block, s.code, eta[i], s.type.value_or(ParameterType::SM)); //TODO check value_or
    //     }
    //     // auto all = oi_.compute_all();
    //     Vec out; out.reserve(obs_ids_.size());
    //     for (auto oid : obs_ids_) {
    //         out.push_back(oi_->compute_observable(oid).front().value); // or from compute_all()
    //     }
    //     return out;
    // }

    // std::map<ObservableId, double> predict(const std::map<ParamId, double>& p, const std::map<ParamId, double>& eta) override {
    //     auto obs = oi_->get_current_observables();
    //     oi_->enable_obs();

    //     // if (p.size()!=p_specs_.size() || eta.size()!=eta_specs_.size())
    //     // throw std::invalid_argument("(p,eta) sizes do not match specs");
    //     for (auto p_elem : p) {
    //         const auto& s = p_elem.first;
    //         oi_->set_param(s.block, s.code, p_elem.second, s.type.value_or(ParameterType::SM)); //TODO check value_or
    //     }
    //     for (auto eta_elem : eta) {
    //         const auto& s = eta_elem.first;
    //         oi_->set_param(s.block, s.code, eta_elem.second, s.type.value_or(ParameterType::SM)); //TODO check value_or
    //     }
    //     // auto all = oi_.compute_all();
    //     std::map<ObservableId, double> out;
    //     for (auto oid : obs_ids_) {
    //         out[oid] = oi_->compute_observable(oid).front().value; // or from compute_all()
    //     }
    //     return out;
    // }

    // std::map<ObservableId, std::vector<ObservableValue>> predict_optimized(const std::map<ParamId, double>& p, const std::map<ParamId, double>& eta) override {
    //     StatParamOptimizerProxy spop = StatParamOptimizerProxy();
    //     // auto obs = oi_->get_current_observables();

    //     // if (p.size()!=p_specs_.size() || eta.size()!=eta_specs_.size())
    //     // throw std::invalid_argument("(p,eta) sizes do not match specs");
    //     for (auto p_elem : p) {
    //         const auto& s = p_elem.first;
    //         spop.set_value(s.block, s.code, p_elem.second);
    //         // oi_->set_param(s.block, s.code, p_elem.second, s.type.value_or(ParameterType::SM)); //TODO check value_or
    //     }
    //     for (auto eta_elem : eta) {
    //         const auto& s = eta_elem.first;
    //         spop.set_value(s.block, s.code, eta_elem.second);
    //         // oi_->set_param(s.block, s.code, eta_elem.second, s.type.value_or(ParameterType::SM)); //TODO check value_or
    //     }
    //     spop.commit();
        
    //     // oi_->enable_obs();
    //     // auto v1 = oi_->compute_observable(obs_ids_[0]).front().value;
    //     // auto v2 = oi_->compute_observable(obs_ids_[0]).front().value;


    //     // double a = oi_->compute_observable(obs_ids_[0]).front().value;
    //     // double b = oi_->compute_observable(obs_ids_[1]).front().value;
    //     // // recompute le premier après avoir calculé le second
    //     // double a2 = oi_->compute_observable(obs_ids_[0]).front().value;

    //     // oi_->compute_all();
    //     // std::map<ObservableId, double> out;
    //     // for (auto oid : obs_ids_) {
    //     //     out[oid] = oi_->compute_observable(oid).front().value; // or from compute_all()
    //     // }
    //     return oi_->compute_all();
    // }

    std::map<ObservableId, std::vector<ObservableValue>> predict_optimized(
        const std::map<ParamId, double>& p,
        const std::map<ParamId, double>& eta) override
    {
        StatParamOptimizerProxy spop = StatParamOptimizerProxy();

        for (auto p_elem : p) {
            const auto& s = p_elem.first;
            spop.set_value(s.block, s.code, p_elem.second);
        }
        for (auto eta_elem : eta) {
            const auto& s = eta_elem.first;
            spop.set_value(s.block, s.code, eta_elem.second);
        }
        spop.commit();

        static std::size_t dbg_call = 0;
        if (dbg_call < 20) {
            StatParameterProxy spp;

            std::cout << "[PDET] call=" << dbg_call << "\n";
            for (const auto& [pid, val] : p) {
                auto p_read = spp.get_param(pid);
                std::cout << "[PDET] requested " << pid
                        << " -> " << val
                        << " | readback same = "
                        << (p_read ? p_read->get_val() : scalar_t(0.0))
                        << "\n";

                auto parts = pid.code.get_parts();
                if (parts.size() >= 4) {
                    LhaID same_sm (parts[0], parts[1], parts[2], 0);
                    LhaID same_bsm(parts[0], parts[1], parts[2], 1);
                    LhaID same_tot(parts[0], parts[1], parts[2], 2);

                    try {
                        auto p_sm  = spp.get_param(pid.block, same_sm);
                        auto p_bsm = spp.get_param(pid.block, same_bsm);
                        auto p_tot = spp.get_param(pid.block, same_tot);

                        std::cout << "[PDET] same triplet in block " << pid.block << "\n";
                        std::cout << "       SM  = " << (p_sm  ? p_sm->get_val()  : scalar_t(0.0)) << "\n";
                        std::cout << "       BSM = " << (p_bsm ? p_bsm->get_val() : scalar_t(0.0)) << "\n";
                        std::cout << "       TOT = " << (p_tot ? p_tot->get_val() : scalar_t(0.0)) << "\n";
                    } catch (...) {
                        std::cout << "[PDET] could not read SM/BSM/TOT triplet for " << pid << "\n";
                    }
                }
            }
        }
        ++dbg_call;

        return oi_->compute_all();
    }

    void compute_observables() const {
        oi_->compute_all();
    };

private:
    std::shared_ptr<ObservableInterface> oi_;
    std::vector<ParamId> p_specs_, eta_specs_;
};