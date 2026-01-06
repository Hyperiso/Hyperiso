#ifndef OBSERVABELEVALUE_H
#define OBSERVABELEVALUE_H

#include "Include.h"
#include <optional>
#include <utility>

struct ObservableValue {
    ObservableId id;
    double value;
    std::optional<std::pair<double, double>> bin;

    ObservableValue(ObservableId id, double value) : id(id), value(value) {};
    ObservableValue(ObservableId id, double value, std::pair<double, double> bin) : ObservableValue(id, value) { this->bin = bin; };
};

#endif // OBSERVABELEVALUE_H
