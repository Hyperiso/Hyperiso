#ifndef PARAMETER_PROXY_H
#define PARAMETER_PROXY_H

#include "IParameterProxy.h"
#include "ParameterProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

/**
 * @file ParameterProxy.h
 * @brief Concrete IParameterProxy backed by ParameterProvider / Parameters.
 *
 * This header defines @ref ParameterProxy, an implementation of
 * @ref IParameterProxy that delegates reads/existence/scale queries to
 * @ref ParameterProvider (and therefore to the underlying @ref Parameters system).
 *
 * Special behavior:
 * - When the proxy is configured with ParameterType::WILSON, reads are made
 *   tolerant: if the parameter does not exist, a default-constructed scalar_t
 *   is returned instead of hard failing.
 *
 * @see IParameterProxy
 * @see ParameterProvider
 * @see Parameters
 */

/**
 * @class ParameterProxy
 * @ingroup DataProvidersModule
 * @brief Parameter access proxy bound to a specific ParameterType.
 *
 * This class wraps a @ref ParameterProvider instance configured for one of:
 * - ParameterType::SM
 * - ParameterType::BSM
 * - ParameterType::WILSON
 *
 * Access policy:
 * - SM/BSM: direct provider access (errors propagate according to provider behavior).
 * - WILSON: returns 0/default value for missing entries (common when coefficients are optional).
 */
class ParameterProxy : public IParameterProxy<std::string, LhaID> {
public:
    /**
     * @brief Constructs a proxy bound to a specific ParameterType.
     *
     * @param type Parameter type to access (SM, BSM, WILSON).
     * @throws via LOG_ERROR if @p type is not allowed.
     */
    ParameterProxy(ParameterType type);

    /**
     * @brief Retrieves a parameter value from (block, id).
     *
     * - For WILSON type: returns a default scalar_t if the entry does not exist.
     * - Otherwise: forwards directly to ParameterProvider.
     */
    scalar_t operator()(const std::string& block, const LhaID& id) const override;
    
    /**
     * @brief Checks whether a parameter exists.
     */
    bool exist(const std::string& block, const LhaID& id) const override;

    /**
     * @brief Returns the block scale.
     */
    double get_scale(const std::string& block) const override;
    
private:
    /// Underlying provider used for actual data access.
    ParameterProvider pp;

     /**
     * @brief Parameter types supported by this proxy.
     *
     * This proxy is intended for physical access layers, thus restricting to
     * SM/BSM/WILSON.
     */
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::BSM, ParameterType::WILSON};
};

#endif 