#include "StatCorrelationProxy.h"

/**
 * @copydoc IStatCorrelationProxy::operator()(const ParamId&, const ParamId&, Type)
 */
double StatCorrelationProxy::operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationProvider::CorrelationType type) {
    return this->cp(pid_1, pid_2, type);
}

/**
 * @copydoc IStatCorrelationProxy::operator()(const Observables&, const Observables&, Type)
 */
double StatCorrelationProxy::operator()(const Observables& pid_1, const Observables& pid_2, CorrelationProvider::CorrelationType type) {
    return this->cp(pid_1, pid_2, type);
}

/**
 * @copydoc IStatCorrelationProxy::operator()(const BinnedObservableId&, const BinnedObservableId&, Type)
 */
double StatCorrelationProxy::operator()(const BinnedObservableId& pid_1, const BinnedObservableId& pid_2, CorrelationProvider::CorrelationType type) {
    return this->cp(pid_1, pid_2, type);
}

