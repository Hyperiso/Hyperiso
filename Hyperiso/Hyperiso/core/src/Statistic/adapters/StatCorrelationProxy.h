#ifndef STAT_CORRELATION_PROXY_H
#define STAT_CORRELATION_PROXY_H

#include "IStatCorrelationProxy.h"
#include "CorrelationProvider.h"
#include "Include.h"
#include "BinnedObservableId.h"
#include "ExperimentObs.h"

/**
 * @file StatCorrelationProxy.h
 * @brief Concrete statistical proxy forwarding correlation queries to CorrelationProvider.
 *
 * This class is a thin read-only adapter used by the statistics layer to access
 * correlation coefficients stored in the backend correlation repository.
 *
 * It supports correlation queries between:
 * - parameters identified by @ref ParamId,
 * - fully explicit experiment-scoped observables identified by @ref ExperimentObs,
 * - observables identified by @ref Observables together with an experiment name,
 * - observables identified by @ref ObservableId together with an experiment name,
 * - binned observables identified by @ref BinnedObservableId together with an experiment name.
 *
 * The proxy does not implement any correlation logic by itself; it simply
 * delegates all requests to an internal @ref CorrelationProvider instance.
 *
 * Typical usage:
 * @code
 *   StatCorrelationProxy corr;
 *
 *   double rho_stat_param = corr(p1, p2, CorrelationProvider::CorrelationType::STAT);
 *
 *   double rho_stat_obs =
 *       corr("LHCb", obs1, obs2, CorrelationProvider::CorrelationType::STAT);
 *
 *   double rho_tot_cross =
 *       corr("LHCb", obs1, "Belle", obs2, CorrelationProvider::CorrelationType::COMBINED);
 * @endcode
 *
 * @see IStatCorrelationProxy
 * @see CorrelationProvider
 * @see CorrelationRepository
 */

/**
 * @class StatCorrelationProxy
 * @brief Concrete implementation of @ref IStatCorrelationProxy for the statistics module.
 *
 * This proxy centralizes correlation access behind a small stable interface,
 * allowing the statistics layer to remain decoupled from the underlying
 * correlation storage / provider implementation.
 *
 * Observable correlations are experiment-scoped: the experiment name is part
 * of the lookup key unless a fully explicit @ref ExperimentObs object is used.
 */
class StatCorrelationProxy : public IStatCorrelationProxy {
public:
    /// Alias for the correlation selector type.
    using Type = CorrelationProvider::CorrelationType;

    /**
     * @brief Returns the correlation coefficient between two parameters.
     *
     * Delegates directly to the internal @ref CorrelationProvider.
     *
     * @param lhs  First parameter identifier.
     * @param rhs  Second parameter identifier.
     * @param type Kind of correlation requested.
     * @return Requested correlation coefficient.
     */
    double operator()(const ParamId& lhs, const ParamId& rhs, Type type) override;

    /**
     * @brief Returns the correlation coefficient between two fully explicit
     *        experiment-scoped observables.
     *
     * Delegates directly to the internal @ref CorrelationProvider.
     *
     * @param lhs  First experiment-scoped observable.
     * @param rhs  Second experiment-scoped observable.
     * @param type Kind of correlation requested.
     * @return Requested correlation coefficient.
     */
    double operator()(const ExperimentObs& lhs, const ExperimentObs& rhs, Type type) override;

    /**
     * @brief Returns the correlation coefficient between two observables
     *        identified by enum in the same experiment.
     *
     * Delegates directly to the internal @ref CorrelationProvider.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param lhs        First observable.
     * @param rhs        Second observable.
     * @param type       Kind of correlation requested.
     * @return Requested correlation coefficient.
     */
    double operator()(const std::string& experiment,
                      const Observables& lhs,
                      const Observables& rhs,
                      Type type) override;

    /**
     * @brief Returns the correlation coefficient between two observables
     *        identified by ObservableId in the same experiment.
     *
     * Delegates directly to the internal @ref CorrelationProvider.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param lhs        First observable identifier.
     * @param rhs        Second observable identifier.
     * @param type       Kind of correlation requested.
     * @return Requested correlation coefficient.
     */
    double operator()(const std::string& experiment,
                      const ObservableId& lhs,
                      const ObservableId& rhs,
                      Type type) override;

    /**
     * @brief Returns the correlation coefficient between two binned observables
     *        in the same experiment.
     *
     * Delegates directly to the internal @ref CorrelationProvider.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param lhs        First binned observable identifier.
     * @param rhs        Second binned observable identifier.
     * @param type       Kind of correlation requested.
     * @return Requested correlation coefficient.
     */
    double operator()(const std::string& experiment,
                      const BinnedObservableId& lhs,
                      const BinnedObservableId& rhs,
                      Type type) override;

    /**
     * @brief Returns the correlation coefficient between two binned observables
     *        coming from possibly different experiments.
     *
     * Delegates directly to the internal @ref CorrelationProvider.
     *
     * @param exp_lhs Experiment associated with the first observable.
     * @param lhs     First binned observable identifier.
     * @param exp_rhs Experiment associated with the second observable.
     * @param rhs     Second binned observable identifier.
     * @param type    Kind of correlation requested.
     * @return Requested correlation coefficient.
     */
    double operator()(const std::string& exp_lhs,
                      const BinnedObservableId& lhs,
                      const std::string& exp_rhs,
                      const BinnedObservableId& rhs,
                      Type type) override;

private:
    CorrelationProvider cp;  ///< Underlying provider used to query the correlation repository.
};

#endif // STAT_CORRELATION_PROXY_H