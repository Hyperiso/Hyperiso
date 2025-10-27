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

    void add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr);

    ObservableInterface& add_observable(Observables obs, QCDOrder order, bool add_dependencies=false) {
        manager->add_obs(obs, order, add_dependencies);
        return *this;
    }

    ObservableInterface& add_observable(ObservableId obs, QCDOrder order, bool add_dependencies=false) {  
        manager->add_obs(obs, order, add_dependencies);
        return *this;
    }

    void add_observables(std::map<Observables, QCDOrder> obss, bool add_dependencies=false) {  
        for (auto &[k, v] : obss) {
            add_observable(k, v, add_dependencies);
        }
    }

    void add_observables(std::map<ObservableId, QCDOrder> obss, bool add_dependencies=false) {  
        for (auto &[k, v] : obss) {
            add_observable(k, v, add_dependencies);
        }
    }

    void add_observables(Decays decay, QCDOrder order, bool add_dependencies=false) {  
        for (auto &obs : DecayMapper::get_observables(decay)) {
            add_observable(obs, order, add_dependencies);
        }
    }

    // TODO : Overload get_observables(DecayId)
    // void add_observables(DecayId decay, QCDOrder order, bool add_dependencies=false) {  
    //     for (auto &obs : DecayMapper::get_observables(decay)) {
    //         add_observable(obs, order, add_dependencies);
    //     }
    // }

    void add_observable_parameter(Observables obs, ParamId pid) {
        manager->add_obs_dep(obs, pid);
    }

    void add_observable_parameter(ObservableId obs, ParamId pid) {
        manager->add_obs_dep(obs, pid);
    }

    void add_observable_parameters(Observables obs, std::unordered_set<ParamId> pids) {
        manager->add_obs_deps(obs, pids);
    }

    void add_observable_parameters(ObservableId obs, std::unordered_set<ParamId> pids) {
        manager->add_obs_deps(obs, pids);
    }

    std::vector<ObservableValue> compute_observable(Observables obs) const {
        return manager->evaluate(obs);
    }

    std::vector<ObservableValue> compute_observable(ObservableId obs) const {
        return manager->evaluate(obs);
    }

    scalar_t compute_uncertainty(Observables obs, UncertaintyType u_type=UncertaintyType::COMBINED) const {
        return manager->get_uncertainty(obs);
    }

    scalar_t compute_uncertainty(ObservableId obs, UncertaintyType u_type=UncertaintyType::COMBINED) const {
        return manager->get_uncertainty(obs);
    }

    std::unordered_map<ParamId, scalar_t> compute_leading_uncertainties(Observables obs, size_t n, UncertaintyType u_type=UncertaintyType::COMBINED) const {
        return manager->get_leading_uncertainties(obs, n);
    }

    std::unordered_map<ParamId, scalar_t> compute_leading_uncertainties(ObservableId obs, size_t n, UncertaintyType u_type=UncertaintyType::COMBINED) const {
        return manager->get_leading_uncertainties(obs, n);
    }

    std::unordered_map<ObservableId, scalar_t> compute_all_uncertainties() const {
        return manager->get_all_uncertainties();
    }

    double compute_chi2() const {
        return manager->get_chi2();
    }

    void remove_observable(Observables id) {
        manager->remove_obs(id);
    }

    void remove_observable(ObservableId id) {
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

    // TODO : same as above
    // void remove_observables(DecayId dec) {
    //     for (Observables id : DecayMapper::get_observables(dec)) {
    //         remove_observable(id);
    //     }
    // };

    scalar_t get_exp_value(Observables id) {
        return manager->get_obs(id)->get_exp_val();
    };

    scalar_t get_exp_value(ObservableId id) {
        return manager->get_obs(id)->get_exp_val();
    };

    scalar_t get_exp_uncertainty(Observables id, UncertaintyType u_type=UncertaintyType::COMBINED) {
        return manager->get_obs(id)->get_exp_uncertainty(u_type);
    }

    scalar_t get_exp_uncertainty(ObservableId id, UncertaintyType u_type=UncertaintyType::COMBINED) {
        return manager->get_obs(id)->get_exp_uncertainty(u_type);
    }

    std::unordered_set<ObservableId> get_current_observables() {
        return manager->get_current_obss();
    }

    std::unordered_map<ObservableId, Estimate> compute_all() {
        return manager->evaluate_all();
    }

    std::unordered_map<ObservableId, Estimate> get_all_exp() {
        std::unordered_map<ObservableId, Estimate> all_exp;
        for (ObservableId id : get_current_observables()) {
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

    void update_gradient(Observables obs) {
        manager->update_gradient(obs);
    }

    void update_gradient(ObservableId obs) {
        manager->update_gradient(obs);
    }
};

#endif