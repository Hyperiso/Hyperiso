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

public:
    ObservableInterface() : manager(ObsManager::GetInstance()) {}

    ObservableInterface add_observable(Observables obs, QCDOrder order) {  
        manager->add_obs(obs, order);
        return *this;
    }

    void add_observables(std::map<Observables, QCDOrder> obss) {  
        for (auto &[k, v] : obss) {
            add_observable(k, v);
        }
    }

    void add_observable_parameter(Observables obs, ParamId pid) {
        manager->add_obs_dep(obs, pid);
    }

    void add_observable_parameters(Observables obs, std::vector<ParamId> pids) {
        manager->add_obs_deps(obs, pids);
    }

    double compute_observable(Observables obs) const {
        return manager->evaluate(obs);
    }

    std::map<Observables, double> compute_all_observables() const {
        return manager->evaluate_all();
    }

    double compute_uncertainty(Observables obs) const {
        return manager->get_uncertainty(obs);
    }

    std::map<ParamId, double> compute_leading_uncertainties(Observables obs, size_t n) const {
        return manager->get_leading_uncertainties(obs, n);
    }

    std::map<Observables, double> compute_all_uncertainties() const {
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

};
