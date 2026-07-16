#include "CorrelationProxy.h"

/**
 * @copydoc IObsCorrelationProxy::operator()(const ParamId&, const ParamId&, Type)
 */
double CorrelationProxy::operator()(const ParamId& pid_1,
                                    const ParamId& pid_2,
                                    CorrelationProvider::CorrelationType type) {
    return this->cp(pid_1, pid_2, type);
}

/**
 * @copydoc IObsCorrelationProxy::operator()(const ExperimentObs&, const ExperimentObs&, Type)
 */
double CorrelationProxy::operator()(const ExperimentObs& oid_1,
                                    const ExperimentObs& oid_2,
                                    CorrelationProvider::CorrelationType type) {
    return this->cp(oid_1, oid_2, type);
}

/**
 * @copydoc IObsCorrelationProxy::operator()(const std::string&, const Observables&, const Observables&, Type)
 */
double CorrelationProxy::operator()(const std::string& experiment,
                                    const Observables& oid_1,
                                    const Observables& oid_2,
                                    CorrelationProvider::CorrelationType type) {
    return this->cp(experiment, oid_1, oid_2, type);
}

/**
 * @copydoc IObsCorrelationProxy::operator()(const std::string&, const ObservableId&, const ObservableId&, Type)
 */
double CorrelationProxy::operator()(const std::string& experiment,
                                    const ObservableId& oid_1,
                                    const ObservableId& oid_2,
                                    CorrelationProvider::CorrelationType type) {
    return this->cp(experiment, oid_1, oid_2, type);
}

/**
 * @copydoc IObsCorrelationProxy::operator()(const std::string&, const BinnedObservableId&, const BinnedObservableId&, Type)
 */
double CorrelationProxy::operator()(const std::string& experiment,
                                    const BinnedObservableId& oid_1,
                                    const BinnedObservableId& oid_2,
                                    CorrelationProvider::CorrelationType type) {
    return this->cp(experiment, oid_1, oid_2, type);
}

/**
 * @copydoc IObsCorrelationProxy::operator()(const std::string&, const BinnedObservableId&, const std::string&, const BinnedObservableId&, Type)
 */
double CorrelationProxy::operator()(const std::string& exp_1,
                                    const BinnedObservableId& oid_1,
                                    const std::string& exp_2,
                                    const BinnedObservableId& oid_2,
                                    CorrelationProvider::CorrelationType type) {
    return this->cp(exp_1, oid_1, exp_2, oid_2, type);
}