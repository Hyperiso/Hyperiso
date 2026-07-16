#ifndef MARTYPARAMETERPROXY_H
#define MARTYPARAMETERPROXY_H

#include "IMartyParameterProxy.h"
#include "ParameterProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

/**
 * @file MartyParameterProxy.h
 * @brief Concrete parameter proxy that reads values from Hyperiso Parameters.
 *
 * This header defines ::MartyParameterProxy, a bridge between the
 * Hyperiso ::Parameters system and the MARTY code-generation layer.
 */

/**
 * @class MartyParameterProxy
 * @ingroup CodeGenerationModule
 * @brief Implementation of IMartyParameterProxy using ParameterProvider.
 *
 * MartyParameterProxy:
 *  - checks that the requested ::ParameterType is allowed (SM or BSM),
 *  - optionally forces SM-only access if the global model is SM,
 *  - retrieves values through a ::ParameterProvider instance.
 */
class MartyParameterProxy : public IMartyParameterProxy<std::string, LhaID> {
public:
    /**
     * @brief Constructs a proxy for a given parameter type.
     *
     * If the global model (queried through ModelAPI / Hyperiso) is SM,
     * the proxy will always use ::ParameterType::SM internally, regardless
     * of the type passed in.
     *
     * @param type Desired parameter type (SM or BSM).
     *
     * @throws LOG_ERROR (ValueError) if @p type is not allowed.
     */
    MartyParameterProxy(ParameterType type);

    /**
     * @brief Retrieves the value of a parameter.
     *
     * @param block Name of the parameter block.
     * @param id    LHA identifier inside the block.
     * @return The parameter value as a scalar.
     */
    scalar_t operator()(const std::string& block, const LhaID& id) const override;
    
private:
    ParameterProvider pp;  ///< Underlying parameter provider.

    /// Allowed types for MARTY access (currently SM and BSM).
    static inline const std::unordered_set<ParameterType> ALLOWED {ParameterType::SM, ParameterType::BSM};
};

#endif // MARTYPARAMETERPROXY_H
