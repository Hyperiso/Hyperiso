#ifndef STAT_CORRELATION_PROXY_H
#define STAT_CORRELATION_PROXY_H

#include "IStatCorrelationProxy.h"
#include "CorrelationProvider.h"
#include "Include.h"
#include "BinnedObservableId.h"

/**
 * @file StatCorrelationProxy.h
 * @brief Concrete statistical proxy forwarding correlation queries to CorrelationProvider.
 *
 * This class is a thin read-only adapter used by the statistics layer to access
 * correlation coefficients stored in the backend correlation repository.
 *
 * It supports correlation queries between:
 * - parameters identified by @ref ParamId,
 * - observables identified by @ref Observables,
 * - binned observables identified by @ref BinnedObservableId.
 *
 * The proxy does not implement any correlation logic by itself; it simply
 * delegates all requests to an internal @ref CorrelationProvider instance.
 *
 * Typical usage:
 * @code
 *   StatCorrelationProxy corr;
 *
 *   double rho_stat = corr(obs1, obs2, CorrelationProvider::CorrelationType::STAT);
 *   double rho_syst = corr(obs1, obs2, CorrelationProvider::CorrelationType::SYST);
 *   double rho_tot  = corr(obs1, obs2, CorrelationProvider::CorrelationType::COMBINED);
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
 */
class StatCorrelationProxy  : public IStatCorrelationProxy {
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
    double operator()(const ParamId&, const ParamId&, Type) override;

    /**
     * @brief Returns the correlation coefficient between two observables.
     *
     * Delegates directly to the internal @ref CorrelationProvider.
     *
     * @param lhs  First observable.
     * @param rhs  Second observable.
     * @param type Kind of correlation requested.
     * @return Requested correlation coefficient.
     */
    double operator()(const Observables&, const Observables&, Type) override;

    /**
     * @brief Returns the correlation coefficient between two binned observables.
     *
     * Delegates directly to the internal @ref CorrelationProvider.
     *
     * @param lhs  First binned observable identifier.
     * @param rhs  Second binned observable identifier.
     * @param type Kind of correlation requested.
     * @return Requested correlation coefficient.
     */
    double operator()(const BinnedObservableId&, const BinnedObservableId&, Type) override;

private:
    CorrelationProvider cp;     /// Underlying provider used to query the correlation repository.
};


#endif 