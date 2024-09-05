#ifndef OBSERVABLESTRATEGY_H
#define OBSERVABLESTRATEGY_H

#include <vector>
#include <string>
#include <memory>
#include <map>
#include "Nuisance.h"
#include "Observables.h"


class ObservableStrategyChi2 {
public:
    virtual ~ObservableStrategyChi2() = default;
    virtual double calculate(const std::vector<Nuisance>& params) = 0;
};


class ObservableStrategy {
public:
    virtual ~ObservableStrategy() = default;
    virtual double calculate(const std::vector<Nuisance>& params) = 0;
};

class SpecificObservable : public ObservableStrategy {
public:
    double calculate(const std::vector<Nuisance>& params) override;
};




class TheoObservable {
public:
    Observables id;
    double value;
    std::unique_ptr<ObservableStrategy> strategy;
    int model;
    int order;
    double scale;
    int wilson_basis;

    TheoObservable(Observables id, std::unique_ptr<ObservableStrategy> strategy, int model, int order, double scale, int wilson_basis);

    void calculate(const std::map<std::string, Nuisance>& params);
};


// class ObservableFactory {
// public:
//     static std::unique_ptr<TheoObservable> createObservable(Observables id, int model, int order, double scale, int wilson_basis);
// };

#endif // OBSERVABLESTRATEGY_H
