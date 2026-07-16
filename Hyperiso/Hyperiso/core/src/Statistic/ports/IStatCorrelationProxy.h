#ifndef ISTAT_CORRELATION_PROXY_H
#define ISTAT_CORRELATION_PROXY_H

#include "CorrelationProvider.h"
#include "Include.h"
#include "ExperimentObs.h"

/**
 * @file IStatCorrelationProxy.h
 * @brief Statistical-layer port for read-only access to correlation coefficients.
 *
 * This interface defines the contract used by the statistics module to query
 * correlations between:
 * - generic parameters identified by @ref ParamId,
 * - fully explicit experiment-scoped observables identified by @ref ExperimentObs,
 * - public observable enums identified by @ref Observables together with an experiment name,
 * - observables identified by @ref ObservableId together with an experiment name,
 * - binned observables identified by @ref BinnedObservableId together with an experiment name.
 *
 * The exact quantity returned is selected through
 * @ref CorrelationProvider::CorrelationType:
 * - @ref CorrelationProvider::CorrelationType::STAT     for statistical correlation,
 * - @ref CorrelationProvider::CorrelationType::SYST     for systematic correlation,
 * - @ref CorrelationProvider::CorrelationType::COMBINED for combined correlation.
 *
 * Implementations are expected to be read-only adapters over the internal
 * correlation repository / provider layer.
 *
 * @see StatCorrelationProxy
 * @see CorrelationProvider
 * @see CorrelationRepository
 */

/**
 * @struct IStatCorrelationProxy
 * @brief Abstract read-only interface for correlation queries in the statistics module.
 *
 * This port provides a uniform callable API so higher-level statistical code
 * can request correlation coefficients without depending directly on the
 * storage/backend layer.
 *
 * Observable correlations are experiment-scoped: the experiment name is part
 * of the lookup key unless a fully explicit @ref ExperimentObs object is used.
 */
struct IStatCorrelationProxy {
    /// Alias for the correlation selector type used by the backend provider.
    using Type = CorrelationProvider::CorrelationType;

    /**
     * @brief Virtual destructor.
     */
    virtual ~IStatCorrelationProxy() = default;

    /**
     * @brief Returns the correlation coefficient between two parameters.
     *
     * Both parameters are identified by fully typed @ref ParamId objects.
     *
     * @param lhs  First parameter identifier.
     * @param rhs  Second parameter identifier.
     * @param type Kind of correlation requested (statistical, systematic, combined).
     * @return Requested correlation coefficient.
     */
    virtual double operator()(const ParamId& lhs, const ParamId& rhs, Type type) = 0;

    /**
     * @brief Returns the correlation coefficient between two fully explicit
     *        experiment-scoped observables.
     *
     * This overload is the most explicit observable interface. Each observable
     * is identified by an @ref ExperimentObs, i.e. by the pair
     * `(experiment, BinnedObservableId)`.
     *
     * @param lhs  First experiment-scoped observable.
     * @param rhs  Second experiment-scoped observable.
     * @param type Kind of correlation requested (statistical, systematic, combined).
     * @return Requested correlation coefficient.
     */
    virtual double operator()(const ExperimentObs& lhs, const ExperimentObs& rhs, Type type) = 0;

    /**
     * @brief Returns the correlation coefficient between two observables
     *        identified by enum in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param lhs        First observable.
     * @param rhs        Second observable.
     * @param type       Kind of correlation requested (statistical, systematic, combined).
     * @return Requested correlation coefficient.
     */
    virtual double operator()(const std::string& experiment,
                              const Observables& lhs,
                              const Observables& rhs,
                              Type type) = 0;

    /**
     * @brief Returns the correlation coefficient between two observables
     *        identified by ObservableId in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param lhs        First observable identifier.
     * @param rhs        Second observable identifier.
     * @param type       Kind of correlation requested (statistical, systematic, combined).
     * @return Requested correlation coefficient.
     */
    virtual double operator()(const std::string& experiment,
                              const ObservableId& lhs,
                              const ObservableId& rhs,
                              Type type) = 0;

    /**
     * @brief Returns the correlation coefficient between two binned observables
     *        in the same experiment.
     *
     * This overload is intended for observable entries carrying an explicit
     * bin definition through @ref BinnedObservableId.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param lhs        First binned observable identifier.
     * @param rhs        Second binned observable identifier.
     * @param type       Kind of correlation requested (statistical, systematic, combined).
     * @return Requested correlation coefficient.
     */
    virtual double operator()(const std::string& experiment,
                              const BinnedObservableId& lhs,
                              const BinnedObservableId& rhs,
                              Type type) = 0;

    /**
     * @brief Returns the correlation coefficient between two binned observables
     *        coming from possibly different experiments.
     *
     * @param exp_lhs Experiment associated with the first observable.
     * @param lhs     First binned observable identifier.
     * @param exp_rhs Experiment associated with the second observable.
     * @param rhs     Second binned observable identifier.
     * @param type    Kind of correlation requested (statistical, systematic, combined).
     * @return Requested correlation coefficient.
     */
    virtual double operator()(const std::string& exp_lhs,
                              const BinnedObservableId& lhs,
                              const std::string& exp_rhs,
                              const BinnedObservableId& rhs,
                              Type type) = 0;
};

#endif // ISTAT_CORRELATION_PROXY_H