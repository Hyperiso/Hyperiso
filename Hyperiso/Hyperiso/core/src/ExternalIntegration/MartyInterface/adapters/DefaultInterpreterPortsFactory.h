#ifndef DEFAULT_INTERPRETER_PORTS_FACTORY_H
#define DEFAULT_INTERPRETER_PORTS_FACTORY_H

#include <memory>

#include "IInterpreterPortsFactory.h"

/**
 * @file DefaultInterpreterPortsFactory.h
 * @brief Declares the default factory for interpreter ports.
 *
 * This factory wires together JSON loaders, mapping adapters and database
 * proxies into a concrete IParameterResolver used by Interpreter.
 */

class IMappingAdapterFactory;
class IParamMappingSource;

/**
 * @class DefaultInterpreterPortsFactory
 * @ingroup MappingModule
 * @brief Default implementation of IInterpreterPortsFactory.
 *
 * This factory:
 *  - constructs a mapping adapter using a provided or default factory,
 *  - loads mapping information through the existing JSON loading pipeline,
 *  - builds MappingDatabaseProxy objects for the model and the SM,
 *  - validates that model/BSM mappings do not collide with SM mapping keys,
 *  - returns a ParameterResolverProxy that uses these databases.
 */
class DefaultInterpreterPortsFactory : public IInterpreterPortsFactory {
public:
    explicit DefaultInterpreterPortsFactory(
        std::shared_ptr<IMappingAdapterFactory> adapterFactory = nullptr,
        std::shared_ptr<IParamMappingSource>    loader         = nullptr
    );

    /**
     * @brief Creates a new IParameterResolver for the given model.
     *
     * @param modelName     Model name.
     * @param modelJsonPath Path to the model or user BSM mapping JSON file.
     * @param smJsonPath    Path to the read-only SM mapping JSON file.
     */
    std::unique_ptr<IParameterResolver>
    makeResolver(const std::string& modelName,
                 const std::string& modelJsonPath,
                 const std::string& smJsonPath) const override;

private:
    std::shared_ptr<IMappingAdapterFactory> adapterFactory_;
    std::shared_ptr<IParamMappingSource>    loader_;
};

#endif
