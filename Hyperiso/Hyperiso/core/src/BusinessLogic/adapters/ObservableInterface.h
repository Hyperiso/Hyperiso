#ifndef OBSERVABLE_INTERFACE_H
#define OBSERVABLE_INTERFACE_H

#include <map>
#include <memory>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <cmath>

#include "Include.h"
#include "ModelEvaluator.h"
#include "Decays.h"
#include "ObsManager.h"

class ObservableInterface {
private:
    std::shared_ptr<ObsManager> manager;

public:
    ObservableInterface();

    ObservableInterface& add_observable(Observables obs, QCDOrder order, bool add_dependencies=false) {  
        manager->add_obs(obs, order, add_dependencies);
        return *this;
    }

    void add_observables(std::map<Observables, QCDOrder> obss, bool add_dependencies=false) {  
        for (auto &[k, v] : obss) {
            add_observable(k, v, add_dependencies);
        }
    }

    void add_observables(Decays decay, QCDOrder order, bool add_dependencies=false) {  
        for (auto &obs : DecayMapper::get_observables(decay)) {
            add_observable(obs, order, add_dependencies);
        }
    }

    void add_observable_parameter(Observables obs, ParamId pid) {
        manager->add_obs_dep(obs, pid);
    }

    void add_observable_parameters(Observables obs, std::unordered_set<ParamId> pids) {
        manager->add_obs_deps(obs, pids);
    }

    scalar_t compute_observable(Observables obs) const {
        return manager->evaluate(obs);
    }

    scalar_t compute_uncertainty(Observables obs, UncertaintyType u_type=UncertaintyType::COMBINED) const {
        return manager->get_uncertainty(obs);
    }

    std::unordered_map<ParamId, scalar_t> compute_leading_uncertainties(Observables obs, size_t n, UncertaintyType u_type=UncertaintyType::COMBINED) const {
        return manager->get_leading_uncertainties(obs, n);
    }

    std::unordered_map<Observables, scalar_t> compute_all_uncertainties() const {
        return manager->get_all_uncertainties();
    }

    double compute_chi2() const {
        return manager->get_chi2();
    }

    void remove_observable(Observables id) {
        manager->remove_obs(id);
    }

    void remove_observables(std::unordered_set<Observables> ids) {
        for (Observables id : ids) {
            remove_observable(id);
        }
    };

    void remove_observables(Decays dec) {
        for (Observables id : DecayMapper::get_observables(dec)) {
            remove_observable(id);
        }
    };

    scalar_t get_exp_value(Observables id) {
        return manager->get_obs(id)->get_exp_val();
    };

    scalar_t get_exp_uncertainty(Observables id, UncertaintyType u_type=UncertaintyType::COMBINED) {
        return manager->get_obs(id)->get_exp_uncertainty(u_type);
    }

    std::unordered_set<Observables> get_current_observables() {
        return manager->get_current_obss();
    }

    std::unordered_map<Observables, Estimate> compute_all() {
        return manager->evaluate_all();
    }

    std::unordered_map<Observables, Estimate> get_all_exp() {
        std::unordered_map<Observables, Estimate> all_exp;
        for (Observables id : get_current_observables()) {
            all_exp.emplace(
                id, 
                Estimate {
                    get_exp_value(id),
                    get_exp_uncertainty(id, UncertaintyType::STAT),
                    get_exp_uncertainty(id, UncertaintyType::SYST), 
                }
            );
        }
        return all_exp;
    }

    void set_param(const std::string& block, int code, double value, ParameterType type) {
        Parameters::GetInstance(type)->setBlockValue(block, code, value);
    }

    double get_param(const std::string& block, int code, ParameterType type) {
        return Parameters::GetInstance(type)->operator()(block, code);
    }

    double get_observable_evaluations(Observables obs) {
        return manager->get_obs_evals(obs);
    }

    void update_gradient(Observables obs) {
        manager->update_gradient(obs);
    }

    template <typename EnumType>
    void set_config_flag(Decays decay_id, EnumType e) {
        manager->set_config_flag(decay_id, e);
    }

    template <typename EnumType>
    void set_config_flag(Observables obs_id, EnumType e) {
        manager->set_config_flag(DecayMapper::get_decay(obs_id), e);
    }
};

#endif