#include "ObservableStrategyChi2.h"


TheoObservable::TheoObservable(Observables id, std::unique_ptr<ObservableStrategy> strategy, int model, int order, double scale, int wilson_basis)
    : id(id), strategy(std::move(strategy)), model(model), order(order), scale(scale), wilson_basis(wilson_basis), value(0.0) {}

void TheoObservable::calculate(const std::map<std::string, Nuisance>& params) {
    std::vector<Nuisance> relevant_params_objects;
    for (const auto& param : params) {
        relevant_params_objects.push_back(param.second);
    }
    value = strategy->calculate(relevant_params_objects);
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

// std::unique_ptr<TheoObservable> ObservableFactory::createObservable(Observables id, int model, int order, double scale, int wilson_basis) {
//     return std::make_unique<TheoObservable>(id, std::make_unique<SpecificObservable>(), model, order, scale, wilson_basis);
// }