#pragma once
#include <vector>
#include <stdexcept>
#include "ports/IModel.h"
#include "ObservableInterface.h"
#include "StatParamOptimizerProxy.h"
#include "StatParameterProxy.h"
#include "IStatParamOptimizerProxy.h"

class ObservableInterfaceAdapterObs final : public IModel {
public:
    ObservableInterfaceAdapterObs(
    std::shared_ptr<ObservableInterface> obs,
    std::vector<ParamId> p_specs,
    std::vector<ParamId> eta_specs)
    : p_specs_(std::move(p_specs)), eta_specs_(std::move(eta_specs)) {
        oi_ = obs;
    }

    ObservableInterfaceAdapterObs(std::shared_ptr<ObservableInterface> obs, std::shared_ptr<IStatParamOptimizerProxy> spop) { oi_ = obs; spop_ = spop;}
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

        // static std::size_t dbg_call = 0;
        // if (has_nonzero_fit_param && dbg_call < 40) {
        //     StatParameterProxy spp(ParameterType::WILSON);

        //     std::cout << "[PDET] call=" << dbg_call << "\n";

        //     for (const auto& [pid, val] : p) {
        //         auto p_read = spp.get_param(pid);
        //         std::cout << "[PDET] requested " << pid
        //                 << " -> " << val
        //                 << " | readback same = "
        //                 << (p_read ? p_read->get_val() : scalar_t(0.0))
        //                 << "\n";

        //         auto parts = pid.code.get_parts();
        //         if (parts.size() >= 4) {
        //             LhaID same_sm (parts[0], parts[1], parts[2], 0);
        //             LhaID same_bsm(parts[0], parts[1], parts[2], 1);
        //             LhaID same_tot(parts[0], parts[1], parts[2], 2);

        //             // Si on est sur un bloc intermédiaire, on lit seulement la composante locale
        //             const std::string suffix_bsm = "__BSM_INTERMEDIATE";
        //             const std::string suffix_sm  = "__SM_INTERMEDIATE";

        //             std::string final_block = pid.block;
        //             bool is_bsm_intermediate = false;
        //             bool is_sm_intermediate  = false;

        //             if (final_block.size() >= suffix_bsm.size() &&
        //                 final_block.compare(final_block.size() - suffix_bsm.size(),
        //                                     suffix_bsm.size(),
        //                                     suffix_bsm) == 0) {
        //                 final_block.erase(final_block.size() - suffix_bsm.size());
        //                 is_bsm_intermediate = true;
        //             } else if (final_block.size() >= suffix_sm.size() &&
        //                     final_block.compare(final_block.size() - suffix_sm.size(),
        //                                         suffix_sm.size(),
        //                                         suffix_sm) == 0) {
        //                 final_block.erase(final_block.size() - suffix_sm.size());
        //                 is_sm_intermediate = true;
        //             }

        //             if (is_bsm_intermediate || is_sm_intermediate) {
        //                 std::cout << "[PDET] intermediate block " << pid.block << "\n";
        //                 std::cout << "       local = "
        //                         << (p_read ? p_read->get_val() : scalar_t(0.0)) << "\n";
        //             }

        //             // Toujours lire le triplet dans le bloc final hadronique
        //             try {
        //                 auto f_sm  = spp.get_param(final_block, same_sm);
        //                 auto f_bsm = spp.get_param(final_block, same_bsm);
        //                 auto f_tot = spp.get_param(final_block, same_tot);

        //                 std::cout << "[PDET] final triplet in block " << final_block << "\n";
        //                 std::cout << "       SM  = " << (f_sm  ? f_sm->get_val()  : scalar_t(0.0)) << "\n";
        //                 std::cout << "       BSM = " << (f_bsm ? f_bsm->get_val() : scalar_t(0.0)) << "\n";
        //                 std::cout << "       TOT = " << (f_tot ? f_tot->get_val() : scalar_t(0.0)) << "\n";
        //             } catch (...) {
        //                 std::cout << "[PDET] could not read final triplet in block "
        //                         << final_block << "\n";
        //             }
        //         }
        //     }
        // }

        auto pred = oi_->compute_all();

        // if (has_nonzero_fit_param && dbg_call < 40) {
        //     std::size_t shown = 0;
        //     for (const auto& [oid, vals] : pred) {
        //         for (const auto& ov : vals) {
        //             std::cout << "[PREDDBG] obs " << ObservableMapper::str(oid)
        //                     << " = " << ov.value;
        //             if (ov.bin.has_value()) {
        //                 std::cout << " bin=[" << ov.bin->first << "," << ov.bin->second << "]";
        //             }
        //             std::cout << "\n";
        //             if (++shown >= 6) break;
        //         }
        //         if (shown >= 6) break;
        //     }
        // }

        // if (has_nonzero_fit_param) {
        //     ++dbg_call;
        // }
        return pred;
    }

    void compute_observables() const {
        oi_->compute_all();
    };

private:
    std::shared_ptr<ObservableInterface> oi_;
    std::shared_ptr<IStatParamOptimizerProxy> spop_;
    std::vector<ParamId> p_specs_, eta_specs_;
};