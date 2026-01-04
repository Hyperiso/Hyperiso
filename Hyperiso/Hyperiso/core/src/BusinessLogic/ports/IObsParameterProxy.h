#ifndef IOBS_PARAMETER_PROXY_H
#define IOBS_PARAMETER_PROXY_H

#include <concepts>

#include "Math.h"
#include "ParameterProvider.h"

/**
 * @file IObsParameterProxy.h
 * @brief CRTP base class for “proxy” accessors used by the observable layer.
 *
 * This class provides a uniform call syntax `proxy(args...)` while delegating
 * the actual implementation to the derived class via CRTP (Curiously Recurring
 * Template Pattern).
 *
 * Template parameters:
 *  - T: derived class type implementing `operator()(Args...)`
 *  - V: identifier/value helper type (kept to match existing design; may be unused here)
 *
 * The `requires HasCallableOperator<T, Args...>` constraint ensures at compile time
 * that the derived type provides an appropriate call operator.
 *
 * Note:
 * - The return type is currently `double` in this interface. Derived implementations
 *   may still return `scalar_t` (implicitly convertible) if desired.
 *
 * @see ObsParameterProxy
 * @see CorrelationProxy
 */
template<typename T, typename V>
class IObsParameterProxy {
public:
    /**
     * @brief Generic forwarding call operator.
     *
     * Forwards all arguments to the derived implementation:
     * `static_cast<T*>(this)->operator()(args...)`.
     */
    template<typename... Args>
    requires HasCallableOperator<T, Args...>
    double operator()(Args&&... args) {
        return static_cast<T*>(this)->operator()(std::forward<Args>(args)...);
}
};

#endif