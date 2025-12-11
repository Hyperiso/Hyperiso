#ifndef DEFAULT_INTERPRETER_PORTS_FACTORY_H
#define DEFAULT_INTERPRETER_PORTS_FACTORY_H

#include <memory>

#include "IInterpreterPortsFactory.h"

/**
 * @file DefaultInterpreterPortsFactory.h
 * @brief Declares the default factory for interpreter ports.
 *
 * This header defines ::DefaultInterpreterPortsFactory, which wires
 * together JSON loaders, mapping adapters and database proxies into
 * a concrete ::IParameterResolver used by ::Interpreter.
 */


class IMappingAdapterFactory;
class IParamMappingSource;

/**
 * @class DefaultInterpreterPortsFactory
 * @ingroup MappingModule
 * @brief Default implementation of IInterpreterPortsFactory.
 *
 * This factory:
 *  - constructs a mapping adapter (using a provided or default factory),
 *  - loads parameter mapping information (e.g. from JSON files),
 *  - builds ::MappingDatabaseProxy objects for the model and the SM,
 *  - and returns a ::ParameterResolverProxy that uses these databases.
 */
class DefaultInterpreterPortsFactory : public IInterpreterPortsFactory {
public:
    /**
     * @brief Constructs a factory with optional dependencies.
     *
     * If @p adapterFactory or @p loader are null, default implementations
     * will be created internally (e.g. ::DefaultMappingAdapterFactory and
     * ::JsonParamMappingAdapter).
     *
     * @param adapterFactory Optional factory for mapping adapters.
     * @param loader         Optional mapping source loader.
     */
    explicit DefaultInterpreterPortsFactory(
        std::shared_ptr<IMappingAdapterFactory> adapterFactory = nullptr,
        std::shared_ptr<IParamMappingSource>    loader         = nullptr
    );

    /**
     * @brief Creates a new ::IParameterResolver for the given model.
     *
     * @param modelName     Model name (e.g. `"THDM"`, `"MSSM"`).
     * @param modelJsonPath Path to the model’s mapping JSON file.
     * @param smJsonPath    Path to the SM mapping JSON file.
     * @return A unique pointer to a ::ParameterResolverProxy instance.
     */
    std::unique_ptr<IParameterResolver>
    makeResolver(const std::string& modelName,
                 const std::string& modelJsonPath,
                 const std::string& smJsonPath) const override;

private:
    std::shared_ptr<IMappingAdapterFactory> adapterFactory_;    ///< Factory used to build mapping adapters.
    std::shared_ptr<IParamMappingSource>    loader_;            ///< Loader used to fetch mapping data.
};

#endif
