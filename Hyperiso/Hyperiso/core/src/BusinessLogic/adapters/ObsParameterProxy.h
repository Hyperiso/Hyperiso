#ifndef OBS_PARAMETER_PROXY_H
#define OBS_PARAMETER_PROXY_H

#include "IObsParameterProxy.h"
#include "ParameterProvider.h"
#include "Include.h"

/**
 * @file ObsParameterProxy.h
 * @brief Read-only proxy to access parameters needed by observable computations.
 *
 * This proxy is meant to be used by the observable/business layer to query:
 *  - SM parameters
 *  - BSM parameters
 *  - Wilson coefficients
 *  - Flavor/decay/observable blocks
 *
 * The proxy enforces an allow-list of parameter types (see @ref ALLOWED).
 *
 * Special behavior for Wilson parameters:
 *  - If the parameter does not exist, returns a default-constructed @ref scalar_t
 *    instead of raising (this is convenient when some Wilson blocks are optional,
 *    group-dependent, or not produced for a given backend).
 *
 * Two access styles are provided:
 *  - by @ref ParamId (typed access, recommended)
 *  - by (block, LhaID) (legacy/low-level access)
 *
 * Typical usage:
 * @code
 *   ObsParameterProxy sm(ParameterType::SM);
 *   double mt = sm("MASS", 6); // depending on your LhaID overloads
 *
 *   ObsParameterProxy wil(ParameterType::WILSON);
 *   scalar_t c7 = wil("B_MATCHING", LhaID(...));
 * @endcode
 *
 * @see ParameterProvider
 * @see IObsParameterProxy
 */
class ObsParameterProxy : public IObsParameterProxy<ParamId, DataType, std::string, LhaID> {
public:
    /**
     * @brief Construct a proxy limited to a given parameter type.
     *
     * If @p type is not in @ref ALLOWED, this triggers a ValueError.
     *
     * @param type The parameter type namespace to read from (default: SM).
     */
    ObsParameterProxy(ParameterType type = ParameterType::SM);

    /**
     * @brief Access a parameter by fully-typed ParamId.
     *
     * If pid.type is missing, a warning is emitted (untyped ParamId usage).
     * If pid.type is not allowed, this triggers a ValueError.
     *
     * For Wilson parameters: returns default @ref scalar_t if parameter does not exist.
     *
     * @param pid    Typed parameter identifier.
     * @param d_type Which value to return (central/value, error, etc.).
     * @return The requested scalar value (or default if missing Wilson parameter).
     */
    scalar_t operator()(const ParamId& pid, DataType d_type=DataType::VALUE) override;

    /**
     * @brief Access a parameter by block name and LHA-like id.
     *
     * For Wilson parameters: if missing, returns default @ref scalar_t.
     * For other types: delegates directly to the provider.
     *
     * @param block  Parameter block name.
     * @param id     Parameter code inside the block.
     * @param d_type Which value to return (central/value, error, etc.).
     * @return The requested scalar value (or default if missing Wilson parameter).
     */
    scalar_t operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE) const override;
    
    /**
     * @brief Returns the underlying Parameter object (shared_ptr) for a given id.
     * @param pid Typed parameter identifier.
     * @return Shared pointer to the parameter (may be null depending on provider policy).
     */
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const override;
private:
    /// Generic provider (can access typed ParamId as-is).
    ParameterProvider pp;

    /// Provider “bound” to a specific ParameterType (set in constructor).
    ParameterProvider pp_with_type;

    /// Allow-list of parameter types that the observable layer is allowed to read.
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::BSM, ParameterType::WILSON, ParameterType::FLAVOR, ParameterType::DECAY, ParameterType::OBSERVABLE};
};

#endif 