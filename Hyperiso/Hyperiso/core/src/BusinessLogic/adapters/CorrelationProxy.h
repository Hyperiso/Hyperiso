#ifndef OBS_CORRELATION_PROXY_H
#define OBS_CORRELATION_PROXY_H

#include "IObsCorrelationProxy.h"
#include "ParameterProvider.h"
#include "CorrelationProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"
#include "ExperimentObs.h"
#include "BinnedObservableId.h"

/**
 * @file CorrelationProxy.h
 * @brief Proxy “user-facing” to retrieve correlations between parameters or experiment-scoped observables.
 *
 * This class is a thin wrapper around @ref CorrelationProvider, exposed through
 * a simple call-operator API so it can be used like a function.
 *
 * It supports multiple identifier types:
 *  - @ref ParamId            : generic typed parameter identifier,
 *  - @ref ExperimentObs      : fully explicit experiment-scoped observable identifier,
 *  - @ref Observables        : observable enum, together with an experiment name,
 *  - @ref ObservableId       : internal observable identifier, together with an experiment name,
 *  - @ref BinnedObservableId : binned observable identifier, together with an experiment name.
 *
 * Observable correlations are experiment-scoped. Therefore, when querying
 * correlations between observables, one must either:
 *  - pass the experiment explicitly, or
 *  - use fully explicit @ref ExperimentObs objects.
 *
 * The @ref CorrelationProvider::CorrelationType selects which quantity is requested:
 *  - statistical correlation,
 *  - systematic correlation,
 *  - combined correlation.
 *
 * Typical usage:
 * @code
 *   CorrelationProxy corr;
 *
 *   double rho_param =
 *       corr(pid1, pid2, CorrelationProvider::CorrelationType::STAT);
 *
 *   double rho_obs =
 *       corr("LHCb", obs1, obs2, CorrelationProvider::CorrelationType::COMBINED);
 *
 *   double rho_explicit =
 *       corr(ExperimentObs{"LHCb", b1}, ExperimentObs{"Belle", b2},
 *            CorrelationProvider::CorrelationType::STAT);
 * @endcode
 *
 * @see CorrelationProvider
 * @see IObsCorrelationProxy
 */
class CorrelationProxy : public IObsCorrelationProxy<
    ParamId,
    ExperimentObs,
    Observables,
    ObservableId,
    BinnedObservableId,
    CorrelationProvider::CorrelationType> {
public:
    /// Alias for the correlation selector type.
    using Type = CorrelationProvider::CorrelationType;

    /**
     * @brief Returns the requested correlation quantity for two parameters.
     *
     * @param pid_1 First parameter identifier.
     * @param pid_2 Second parameter identifier.
     * @param type  Correlation quantity requested.
     * @return Requested correlation quantity.
     */
    double operator()(const ParamId& pid_1, const ParamId& pid_2, Type type) override;

    /**
     * @brief Returns the requested correlation quantity for two fully explicit
     *        experiment-scoped observables.
     *
     * @param oid_1 First experiment-scoped observable.
     * @param oid_2 Second experiment-scoped observable.
     * @param type  Correlation quantity requested.
     * @return Requested correlation quantity.
     */
    double operator()(const ExperimentObs& oid_1, const ExperimentObs& oid_2, Type type) override;

    /**
     * @brief Returns the requested correlation quantity for two observables
     *        identified by enum in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First observable (enum).
     * @param oid_2 Second observable (enum).
     * @param type  Correlation quantity requested.
     * @return Requested correlation quantity.
     */
    double operator()(const std::string& experiment,
                      const Observables& oid_1,
                      const Observables& oid_2,
                      Type type) override;

    /**
     * @brief Returns the requested correlation quantity for two internal
     *        observable identifiers in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First observable identifier.
     * @param oid_2 Second observable identifier.
     * @param type  Correlation quantity requested.
     * @return Requested correlation quantity.
     */
    double operator()(const std::string& experiment,
                      const ObservableId& oid_1,
                      const ObservableId& oid_2,
                      Type type) override;

    /**
     * @brief Returns the requested correlation quantity for two binned
     *        observables in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First binned observable identifier.
     * @param oid_2 Second binned observable identifier.
     * @param type  Correlation quantity requested.
     * @return Requested correlation quantity.
     */
    double operator()(const std::string& experiment,
                      const BinnedObservableId& oid_1,
                      const BinnedObservableId& oid_2,
                      Type type) override;

    /**
     * @brief Returns the requested correlation quantity for two binned
     *        observables coming from possibly different experiments.
     *
     * @param exp_1 Experiment associated with the first observable.
     * @param oid_1 First binned observable identifier.
     * @param exp_2 Experiment associated with the second observable.
     * @param oid_2 Second binned observable identifier.
     * @param type  Correlation quantity requested.
     * @return Requested correlation quantity.
     */
    double operator()(const std::string& exp_1,
                      const BinnedObservableId& oid_1,
                      const std::string& exp_2,
                      const BinnedObservableId& oid_2,
                      Type type) override;

private:
    /// Underlying provider that actually computes / looks up correlations.
    CorrelationProvider cp;
};

#endif // OBS_CORRELATION_PROXY_H