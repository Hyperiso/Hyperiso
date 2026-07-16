#ifndef IMAPPINGADAPTERFACTORY_H
#define IMAPPINGADAPTERFACTORY_H

#include <memory>
#include <string>

#include "IMappingDatabaseAdapter.h"
#include "IParamMappingSource.h"

/**
 * @file IMappingAdapterFactory.h
 * @brief Interface for factories creating mapping database adapters.
 *
 * This abstraction allows higher-level code to request a fully configured
 * IMappingDatabaseAdapter without knowing how it is constructed internally
 * (file format, loader type, etc.).
 */

/**
 * @class IMappingAdapterFactory
 * @brief Factory interface for building IMappingDatabaseAdapter instances.
 *
 * A concrete factory typically:
 *   - chooses or creates an IParamMappingSource (e.g. JsonParamMappingAdapter),
 *   - uses it to populate a concrete MappingDatabase,
 *   - wraps that database inside a IMappingDatabaseAdapter implementation.
 *
 * Example implementation:
 *   - DefaultMappingAdapterFactory.
 */
class IMappingAdapterFactory {
public:
    virtual ~IMappingAdapterFactory() = default;

    /**
     * @brief Creates a mapping database adapter from an underlying JSON-like file.
     *
     * @param instanceName  Human-readable name for the mapping database instance.
     * @param jsonFilePath  Path to the JSON (or JSON-like) file containing
     *                      parameter mapping definitions.
     * @param loader        Optional mapping source responsible for reading
     *                      and interpreting the file. If null, the concrete
     *                      factory is expected to instantiate a sensible
     *                      default (e.g. JsonParamMappingAdapter + JSONParser).
     *
     * @return Shared pointer to a fully-initialized IMappingDatabaseAdapter.
     *
     * @throws std::runtime_error if the underlying database or adapter cannot
     *         be created or initialized.
     */

    virtual std::shared_ptr<IMappingDatabaseAdapter> create(
        const std::string& instanceName,
        const std::string& jsonFilePath,
        std::shared_ptr<IParamMappingSource> loader
    ) const = 0;
};

#endif // IMAPPINGADAPTERFACTORY_H
