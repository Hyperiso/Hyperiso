#include "CorrelationProvider.h"

double CorrelationProvider::operator()(const ParamId &pid_1, const ParamId &pid_2, CorrelationType type) {
    // TODO, return 0 if correlation not present in repo
    return (double)(pid_1 == pid_2);
}

double CorrelationProvider::operator()(const Observables &pid_1, const Observables &pid_2, CorrelationType type) {
    // TODO
    return (double)(pid_1 == pid_2);
}
