#ifndef IOBS_CORRELATION_PROXY_H
#define IOBS_CORRELATION_PROXY_H

/**
 * @file IObsCorrelationProxy.h
 * @brief Interface for retrieving correlation-like quantities between observables/parameters.
 *
 * Observable fits and combinations often require correlation information (or covariance).
 * This interface provides a small, user-friendly API to query correlation-like quantities
 * between two entities identified by various id types.
 *
 * Template parameters:
 *  - T : generic parameter identifier type (e.g. @ref ParamId)
 *  - V : public observable identifier type (e.g. @ref Observables enum)
 *  - U : internal observable identifier type (e.g. @ref ObservableId)
 *  - W : correlation quantity selector type (e.g. CorrelationProvider::CorrelationType)
 *
 * Implementations are expected to forward these queries to a correlation backend
 * (e.g. @ref CorrelationProvider), possibly applying model-specific logic.
 *
 * @see CorrelationProxy
 * @see CorrelationProvider
 */
template<typename T, typename V, typename U, typename W>
class IObsCorrelationProxy {
public:
    /**
     * @brief Query correlation-like quantity for two parameters identified by T.
     *
     * @param pid_1 First id.
     * @param pid_2 Second id.
     * @param type  Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const T& pid_1, const T& pid_2, W type) = 0;

    /**
     * @brief Query correlation-like quantity for two observables (public ids).
     *
     * @param pid_1 First observable id (public enum).
     * @param pid_2 Second observable id (public enum).
     * @param type  Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const V& pid_1, const V& pid_2, W type) = 0;

    /**
     * @brief Query correlation-like quantity for two observables (internal ids).
     *
     * @param pid_1 First internal observable id.
     * @param pid_2 Second internal observable id.
     * @param type  Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const U& pid_1, const U& pid_2, W type) = 0;
};

#endif