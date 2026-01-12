#ifndef OBS_CORRELATION_PROXY_H
#define OBS_CORRELATION_PROXY_H

#include "IObsCorrelationProxy.h"
#include "ParameterProvider.h"
#include "CorrelationProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

/**
 * @file CorrelationProxy.h
 * @brief Proxy “user-facing” to retrieve correlations/covariances between observables.
 *
 * This class is a thin wrapper around @ref CorrelationProvider, exposed through
 * a simple call-operator API so it can be used like a function.
 *
 * It supports multiple identifier types:
 *  - @ref ParamId     : generic typed parameter identifier (block + code + type)
 *  - @ref Observables : observable enum (high-level public API)
 *  - @ref ObservableId: internal observable identifier (when available)
 *
 * The @ref CorrelationProvider::CorrelationType selects which correlation-like object
 * is requested (e.g. correlation coefficient, covariance, etc. depending on the provider).
 *
 * Typical usage:
 * @code
 *   CorrelationProxy corr;
 *   double rho = corr(Observables::BR_B_to_Xs_gamma, Observables::ACP_B_to_Xs_gamma,
 *                     CorrelationProvider::CorrelationType::CORRELATION);
 * @endcode
 *
 * @see CorrelationProvider
 * @see IObsParameterProxy
 */
class CorrelationProxy : public IObsCorrelationProxy<ParamId, Observables, ObservableId, CorrelationProvider::CorrelationType> {
public:
    /**
     * @brief Returns correlation-like quantity for two parameters identified by ParamId.
     * @param pid_1 First parameter id.
     * @param pid_2 Second parameter id.
     * @param type  Correlation quantity requested (provider-defined).
     * @return Requested correlation quantity.
     */
    double operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationProvider::CorrelationType type);

    /**
     * @brief Returns correlation-like quantity for two observables (public enum IDs).
     * @param pid_1 First observable (enum).
     * @param pid_2 Second observable (enum).
     * @param type  Correlation quantity requested (provider-defined).
     * @return Requested correlation quantity.
     */

    double operator()(const Observables& pid_1, const Observables& pid_2, CorrelationProvider::CorrelationType type);

    /**
     * @brief Returns correlation-like quantity for two internal observable identifiers.
     * @param pid_1 First observable id.
     * @param pid_2 Second observable id.
     * @param type  Correlation quantity requested (provider-defined).
     * @return Requested correlation quantity.
     */
    double operator()(const ObservableId& pid_1, const ObservableId& pid_2, CorrelationProvider::CorrelationType type);

    
private:
    /// Underlying provider that actually computes/looks up correlations.
    CorrelationProvider cp;
};

#endif 