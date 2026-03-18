#ifndef IOBS_CORRELATION_PROXY_H
#define IOBS_CORRELATION_PROXY_H

#include <string>

/**
 * @file IObsCorrelationProxy.h
 * @brief Interface for retrieving correlation-like quantities between parameters and observables.
 *
 * Observable fits and combinations often require correlation information
 * (or covariance-like quantities). This interface provides a small,
 * user-friendly API to query such quantities between two entities identified
 * by various id types.
 *
 * Template parameters:
 *  - T : generic parameter identifier type (e.g. @ref ParamId)
 *  - X : fully explicit experiment-scoped observable identifier
 *        (e.g. @ref ExperimentObs)
 *  - V : public observable identifier type (e.g. @ref Observables enum)
 *  - U : internal observable identifier type (e.g. @ref ObservableId)
 *  - B : binned observable identifier type (e.g. @ref BinnedObservableId)
 *  - W : correlation quantity selector type
 *        (e.g. CorrelationProvider::CorrelationType)
 *
 * Implementations are expected to forward these queries to a correlation
 * backend (e.g. @ref CorrelationProvider), possibly applying model-specific
 * logic.
 *
 * Observable correlations are experiment-scoped: the experiment name is part
 * of the lookup key unless a fully explicit identifier of type @p X is used.
 *
 * This interface is independent from @ref IObsParameterProxy and therefore
 * does not affect proxies such as @ref ObsParameterProxy.
 *
 * @see CorrelationProxy
 * @see CorrelationProvider
 */
template<typename T, typename X, typename V, typename U, typename B, typename W>
class IObsCorrelationProxy {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IObsCorrelationProxy() = default;

    /**
     * @brief Query a correlation-like quantity for two parameters identified by @p T.
     *
     * @param pid_1 First parameter identifier.
     * @param pid_2 Second parameter identifier.
     * @param type  Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const T& pid_1, const T& pid_2, W type) = 0;

    /**
     * @brief Query a correlation-like quantity for two fully explicit
     *        experiment-scoped observables identified by @p X.
     *
     * @param pid_1 First experiment-scoped observable identifier.
     * @param pid_2 Second experiment-scoped observable identifier.
     * @param type  Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const X& pid_1, const X& pid_2, W type) = 0;

    /**
     * @brief Query a correlation-like quantity for two public observables
     *        in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param pid_1      First observable identifier (public enum).
     * @param pid_2      Second observable identifier (public enum).
     * @param type       Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const std::string& experiment,
                              const V& pid_1,
                              const V& pid_2,
                              W type) = 0;

    /**
     * @brief Query a correlation-like quantity for two internal observables
     *        in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param pid_1      First internal observable identifier.
     * @param pid_2      Second internal observable identifier.
     * @param type       Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const std::string& experiment,
                              const U& pid_1,
                              const U& pid_2,
                              W type) = 0;

    /**
     * @brief Query a correlation-like quantity for two binned observables
     *        in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param pid_1      First binned observable identifier.
     * @param pid_2      Second binned observable identifier.
     * @param type       Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const std::string& experiment,
                              const B& pid_1,
                              const B& pid_2,
                              W type) = 0;

    /**
     * @brief Query a correlation-like quantity for two binned observables
     *        coming from possibly different experiments.
     *
     * @param exp_1 First experiment name.
     * @param pid_1 First binned observable identifier.
     * @param exp_2 Second experiment name.
     * @param pid_2 Second binned observable identifier.
     * @param type  Which correlation-like quantity is requested.
     * @return Requested correlation-like quantity (provider-defined semantics).
     */
    virtual double operator()(const std::string& exp_1,
                              const B& pid_1,
                              const std::string& exp_2,
                              const B& pid_2,
                              W type) = 0;
};

#endif // IOBS_CORRELATION_PROXY_H