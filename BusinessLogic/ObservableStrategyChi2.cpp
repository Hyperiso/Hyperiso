#include "ObservableStrategyChi2.h"


Observable::Observable(const std::string& name, std::unique_ptr<ObservableStrategy> strategy, const std::vector<std::string>& relevant_params)
    : name(name), strategy(std::move(strategy)), relevant_parameters(relevant_params), value(0.0) {}

void Observable::calculate(const std::map<std::string, Nuisance>& params) {
    std::vector<Nuisance> relevant_params;
    for (const auto& param_name : relevant_parameters) {
        relevant_params.push_back(params.at(param_name));
    }
    value = strategy->calculate(relevant_params);
}

double SpecificObservable::calculate(const std::vector<Nuisance>& params) {
    double result = 0.0;
    for (const auto& param : params) {
        if (param.name == "alphas_MZ") {
            result += param.central_value;
        }
        if (param.name == "mass_b") {
            result += param.central_value;
        }
        if (param.name == "AI_BKstargamma") {
            result += param.central_value;
        }
        if (param.name == "BR_BXsgamma") {
            result += param.central_value;
        }
    }
    return result;
}

std::unique_ptr<Observable> ObservableFactory::createObservable(const std::string& type, const std::string& name, const std::vector<std::string>& relevant_params) {
    if (type == "SpecificObservable") {
        return std::make_unique<Observable>(name, std::make_unique<SpecificObservable>(), relevant_params);
    }
    return nullptr;
}