#ifndef IOBS_PARAMETER_PROXY_H
#define IOBS_PARAMETER_PROXY_H

#include <memory>

#include "Math.h"
#include "ParameterProvider.h"

/**
 * @file IObsParameterProxy.h
 * @brief Interface for read-only access to parameters used by observable computations.
 *
 * The observable layer (a.k.a. business/phenomenology layer) needs to *read* many kinds
 * of parameters produced by the core:
 *  - SM and BSM inputs (masses, couplings, ...)
 *  - Wilson coefficients (matching/hadronic blocks, various contributions/orders)
 *  - Flavor / decay / observable blocks computed upstream
 *
 * This interface abstracts the retrieval mechanism, so observable code does not
 * depend on a concrete storage/backend implementation.
 *
 * Implementations typically wrap a @ref ParameterProvider (or a proxy thereof),
 * and may enforce business rules (allowed parameter types, default behavior
 * for missing Wilson parameters, etc.).
 *
 * Template parameters:
 *  - T : typed parameter identifier type (e.g. @ref ParamId)
 *  - V : “data selector” type (e.g. @ref DataType for VALUE/ERROR/...)
 *  - U : block identifier type (e.g. std::string)
 *  - W : in-block code type (e.g. @ref LhaID)
 *
 * Two access styles are supported:
 *  1) Fully typed access via (T pid, V d_type) — recommended.
 *  2) Low-level access via (U block, W id, V d_type) — useful for legacy calls.
 *
 * Implementations may also expose the underlying @ref Parameter object when needed.
 *
 * @see ObsParameterProxy
 * @see ParameterProvider
 * @see ParamId
 * @see LhaID
 */
template<typename T, typename V, typename U, typename W>
class IObsParameterProxy {
public:
    /**
     * @brief Retrieve a parameter value using a typed identifier.
     *
     * @param pid     Typed parameter identifier (usually includes type + block + code).
     * @param d_type  Which value component to retrieve (central value, error, etc.).
     * @return The requested value as @ref scalar_t.
     *
     * @note Implementations may:
     *  - warn on missing type information,
     *  - restrict allowed parameter namespaces (SM/BSM/WILSON/...),
     *  - return a default-constructed @ref scalar_t when a parameter is absent
     *    (commonly done for optional Wilson coefficients).
     */
    virtual scalar_t operator()(const T& pid, V d_type) = 0;

    /**
     * @brief Retrieve a parameter value using (block, code) addressing.
     *
     * @param block   Block name.
     * @param id      In-block identifier (LHA-like code).
     * @param d_type  Which value component to retrieve (central value, error, etc.).
     * @return The requested value as @ref scalar_t.
     *
     * @note Implementations may apply special missing-parameter behavior
     * (e.g. Wilson blocks may be sparse depending on group/backend).
     */
    virtual scalar_t operator()(const U& block, const W& id, V d_type) const = 0;

    /**
     * @brief Retrieve the underlying Parameter object for a given id.
     *
     * This is useful when downstream code needs metadata beyond the scalar value
     * (mode, uncertainty model, provenance, etc.).
     *
     * @param pid Typed parameter identifier.
     * @return Shared pointer to the underlying @ref Parameter.
     */
    virtual std::shared_ptr<Parameter> get_parameter(const T& pid) const = 0;
};

#endif // IOBS_PARAMETER_PROXY_H