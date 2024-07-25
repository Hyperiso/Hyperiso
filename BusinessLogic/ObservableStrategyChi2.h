#ifndef OBSERVABLESTRATEGY_H
#define OBSERVABLESTRATEGY_H

#include <vector>
#include <string>
#include <memory>
#include <map>
#include "Nuisance.h"

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

class Observable {
public:
    std::string name;
    double value;
    std::unique_ptr<ObservableStrategy> strategy;
    std::vector<std::string> relevant_parameters;

    Observable(const std::string& name, std::unique_ptr<ObservableStrategy> strategy, const std::vector<std::string>& relevant_params);

    void calculate(const std::map<std::string, Nuisance>& params);
};

class ObservableFactory {
public:
    static std::unique_ptr<Observable> createObservable(const std::string& type, const std::string& name, const std::vector<std::string>& relevant_params);
};

#endif // OBSERVABLESTRATEGY_H
