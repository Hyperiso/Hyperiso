#include "CorrelationProxy.h"

double CorrelationProxy::operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationProvider::CorrelationType type) {
    return this->cp(pid_1, pid_2, type);
}
double CorrelationProxy::operator()(const Observables& pid_1, const Observables& pid_2, CorrelationProvider::CorrelationType type) {
    return this->cp(pid_1, pid_2, type);
}
