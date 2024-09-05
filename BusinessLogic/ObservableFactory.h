#ifndef OBSERVABLEFACTORY_H
#define OBSERVABLEFACTORY_H

#include <memory>
#include "Observable.h"
#include "Observables.h"
#include "Nuisance.h"

class ObservableFactory {
public:
    static std::unique_ptr<Observable> createObservable(Observables id, int model, int order, double scale, int wilson_basis);
    static std::unique_ptr<Nuisance> createNuisance(const std::string& name, double central_value, double deviation, double sys_error);
};

#endif // OBSERVABLEFACTORY_H
