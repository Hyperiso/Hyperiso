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
#include "Decays.h"
#include "ObsManager.h"
#include "ObsUseMarty.h"
#include "WilsonFreezer.h"

class ObservableInterface {
private:
    std::shared_ptr<ObsManager> manager;
    std::shared_ptr<ObservablePortsConfig> ports;

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

    std::unordered_map<ObservableId, std::vector<ObservableValue>> compute_all() {
        return manager->evaluate_all();
    }


    void set_param(const std::string& block, LhaID code, double value, ParameterType type) {
        Parameters::GetInstance(type)->setBlockValue(block, code, value);
    }

    double get_param(const std::string& block, LhaID code, ParameterType type) {
        return Parameters::GetInstance(type)->operator()(block, code);
    }

    std::unordered_set<ParamId> get_all_ops_deps(ObservableId id) {
        return manager->get_all_ops_deps(id);
    }

    std::unordered_set<ParamId> get_all_ops_deps(Observables id) {
        return manager->get_all_ops_deps(ObservableMapper::to_id(id));
    }

    void reload_params() {
        manager->reload_params();
    }

    void enable_obs() {
        manager->enable_obs();
    }
    void set_decay_config(Decays dec, std::any config) {
        manager->set_decay_config(dec, config);  
    };

    ObservablePortsConfig& get_ports() {return this->manager->get_ports();}
};

#endif