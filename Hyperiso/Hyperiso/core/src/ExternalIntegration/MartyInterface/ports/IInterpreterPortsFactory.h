#ifndef IINTERPRETER_PORTS_FACTORY_H
#define IINTERPRETER_PORTS_FACTORY_H

#include <memory>
#include <string>

#include "IParameterResolver.h"

/**
 * @file IInterpreterPortsFactory.h
 * @brief Declares a factory interface for building parameter resolvers.
 *
 * The factory provides a way to construct ::IParameterResolver instances
 * for a given model and a pair of JSON mapping files.
 */

/**
 * @class IInterpreterPortsFactory
 * @ingroup MappingModule
 * @brief Factory interface for creating parameter resolver objects.
 *
 * Implementations (e.g. ::DefaultInterpreterPortsFactory) encapsulate
 * the construction of mapping databases and the associated
 * ::IParameterResolver used by ::Interpreter.
 */
class IInterpreterPortsFactory {
public:
    virtual ~IInterpreterPortsFactory() = default;

    /**
     * @brief Creates a new ::IParameterResolver instance.
     *
     * @param modelName     Logical name of the model (e.g. `"SM"`, `"THDM"`).
     * @param modelJsonPath Path to the JSON mapping for @p modelName.
     * @param smJsonPath    Path to the JSON mapping for the SM.
     * @return A fresh, owning pointer to an ::IParameterResolver.
     */
    virtual std::unique_ptr<IParameterResolver>
    makeResolver(const std::string& modelName,
                 const std::string& modelJsonPath,
                 const std::string& smJsonPath) const = 0;
};

#endif 
