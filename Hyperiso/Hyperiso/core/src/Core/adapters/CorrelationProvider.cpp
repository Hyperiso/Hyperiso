#include "CorrelationProvider.h"

double CorrelationProvider::operator()(const ParamId &pid_1, const ParamId &pid_2, CorrelationType type) {
    return 0.0;
}

double CorrelationProvider::operator()(const Observables &pid_1, const Observables &pid_2, CorrelationType type) {
    return 0.0;
}
