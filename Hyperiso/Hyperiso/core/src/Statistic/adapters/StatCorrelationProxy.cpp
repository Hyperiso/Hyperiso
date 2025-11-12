#include "StatCorrelationProxy.h"

double StatCorrelationProxy::operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationProvider::CorrelationType type) {
    // if (!this->cp.exists(pid_1, pid_2, type)) {
    //     return 0;
    // }
    return this->cp(pid_1, pid_2, type);
}
double StatCorrelationProxy::operator()(const Observables& pid_1, const Observables& pid_2, CorrelationProvider::CorrelationType type) {
    // if (!this->cp.exists(pid_1, pid_2, type)) {
    //     return 0;
    // }
    return this->cp(pid_1, pid_2, type);
}

double StatCorrelationProxy::operator()(const ObservableId& pid_1, const ObservableId& pid_2, CorrelationProvider::CorrelationType type) {
    // if (!this->cp.exists(pid_1, pid_2, type)) {
    //     return 0;
    // }
    return this->cp(pid_1, pid_2, type);
}

