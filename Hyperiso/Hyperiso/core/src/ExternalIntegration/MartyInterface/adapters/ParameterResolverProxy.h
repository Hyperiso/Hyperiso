#ifndef PARAMETER_RESOLVER_PROXY_H
#define PARAMETER_RESOLVER_PROXY_H

#include <unordered_map>
#include <string>
#include <memory>
#include <iomanip>

#include "IParameterResolver.h"
#include "IMappingDatabasePort.h"
#include "Logger.h"

/**
 * @file ParameterResolverProxy.h
 * @brief Declares a proxy implementation of IParameterResolver.
 *
 * This proxy combines two mapping databases (SM and model) to resolve
 * parameters, reporting errors when a BSM parameter is used in SM-only
 * calculations or when mappings are missing.
 */

/**
 * @class ParameterResolverProxy
 * @ingroup MappingModule
 * @brief Resolves parameters by querying SM and model mapping databases.
 *
 * Resolution strategy:
 *  - First look up the parameter in the SM database.
 *  - If not found, try the model database.
 *  - If still not found, emit a mapping error.
 *  - If a BSM parameter is found while @p modelIsSM is true, emit a logic error.
 */
class ParameterResolverProxy final : public IParameterResolver {
public:
    /**
     * @brief Constructs a resolver from SM and model databases.
     *
     * @param sm    Shared pointer to the SM mapping database port.
     * @param model Shared pointer to the model-specific mapping database port.
     */
    ParameterResolverProxy(std::shared_ptr<IMappingDatabasePort> sm,
                           std::shared_ptr<IMappingDatabasePort> model);

    /// @copydoc IParameterResolver::resolve()
    std::unordered_map<std::string, ResolvedParam>
    resolve(const std::vector<Extractor::Parameter>& params,
            bool modelIsSM) const override;
    
    /// @copydoc IParameterResolver::clone()
    std::unique_ptr<IParameterResolver> clone() const override;

private:
    std::shared_ptr<IMappingDatabasePort> sm_;      ///< SM mapping database.
    std::shared_ptr<IMappingDatabasePort> model_;   ///< Model mapping database.
};

#endif // PARAMETER_RESOLVER_PROXY_H
