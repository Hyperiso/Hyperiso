#ifndef OBSERVABLE_INTERFACE_H
#define OBSERVABLE_INTERFACE_H

#include <map>
#include <memory>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <cmath>

#include "General.h"
#include "ModelEvaluator.h"
#include "Decays.h"
#include "ObsManager.h"

class ObservableInterface {
private:
    std::shared_ptr<ObsManager> manager;

    void disable_decays() const {
        manager->disable_decays();
    }

public:
    ObservableInterface();

    ObservableInterface& add_observable(Observables obs, QCDOrder order, bool add_dependencies=false) {  
        // disable_decays();
        manager->add_obs(obs, order, add_dependencies);
        return *this;
    }

    void add_observables(std::map<Observables, QCDOrder> obss, bool add_dependencies=false) {  
        for (auto &[k, v] : obss) {
            add_observable(k, v, add_dependencies);
        }
    }

    void add_observable_parameter(Observables obs, ParamId pid) {
        // disable_decays();
        manager->add_obs_dep(obs, pid);
    }

    void add_observable_parameters(Observables obs, std::unordered_set<ParamId> pids) {
        // disable_decays();
        manager->add_obs_deps(obs, pids);
    }

    scalar_t compute_observable(Observables obs) const {
        // disable_decays();
        return manager->evaluate(obs);
    }

    std::unordered_map<Observables, scalar_t> compute_all_observables() const {
        // disable_decays();
        return manager->evaluate_all();
    }

    scalar_t compute_uncertainty(Observables obs) const {
        // disable_decays();
        return manager->get_uncertainty(obs);
    }

    std::unordered_map<ParamId, scalar_t> compute_leading_uncertainties(Observables obs, size_t n) const {
        // disable_decays();
        return manager->get_leading_uncertainties(obs, n);
    }

    std::unordered_map<Observables, scalar_t> compute_all_uncertainties() const {
        return manager->get_all_uncertainties();
    }

    double compute_chi2() const {
        return manager->get_chi2();
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

    std::unordered_set<Observables> get_current_obss() {
        return manager->get_current_obss();
    }

    void update_gradient(Observables obs) {
        // disable_decays();
        manager->update_gradient(obs);
    }

};

#endif