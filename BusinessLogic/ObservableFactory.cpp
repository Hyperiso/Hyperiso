#include "ObservableFactory.h"

std::unique_ptr<Observable> ObservableFactory::createObservable(Observables id, int model, int order, double scale, int wilson_basis) {
    return std::make_unique<Observable>(id, scale, order, model, wilson_basis);
}

std::unique_ptr<Nuisance> ObservableFactory::createNuisance(const std::string& name, double central_value, double deviation, double sys_error) {
    return std::make_unique<Nuisance>(name, central_value, deviation, sys_error);
}
