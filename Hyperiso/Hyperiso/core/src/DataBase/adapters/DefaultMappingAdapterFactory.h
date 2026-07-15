#ifndef DEFAULTMAPPINGADAPTERFACTORY_H
#define DEFAULTMAPPINGADAPTERFACTORY_H

#include "IMappingAdapterFactory.h"
#include "DefaultMappingDatabaseAdapter.h"
#include "JsonParamMappingAdapter.h"
#include "JsonParser.h"
#include "MappingDatabase.h"

/**
 * @file DefaultMappingAdapterFactory.h
 * @brief Factory for creating default mapping database adapters.
 *
 * This factory:
 *   - optionally constructs a JsonParamMappingAdapter with a JSONParser,
 *   - builds a MappingDatabase using the given instance name and JSON file,
 *   - wraps it in a DefaultMappingDatabaseAdapter.
 */

/**
 * @class DefaultMappingAdapterFactory
 * @brief Default implementation of IMappingAdapterFactory.
 *
 * This factory is meant as a convenient, ready-to-use adapter builder for
 * the common case:
 *   - parameter mapping stored in a JSON file,
 *   - loaded via JsonParamMappingAdapter + JSONParser,
 *   - stored in MappingDatabase,
 *   - exposed via DefaultMappingDatabaseAdapter.
 */
class DefaultMappingAdapterFactory : public IMappingAdapterFactory {
public:
    /**
     * @brief Creates a mapping database adapter from a JSON file.
     *
     * @param instanceName  Logical name for the mapping database.
     * @param jsonFilePath  Path to the JSON (or JSON-like) mapping file.
     * @param loader        Optional parameter loader; if null, a
     *                      JsonParamMappingAdapter with a JSONParser will
     *                      be constructed.
     *
     * @return Shared pointer to an IMappingDatabaseAdapter wrapping
     *         a MappingDatabase instance.
     *
     * @throws std::runtime_error if the MappingDatabase cannot be created.
     */
    std::shared_ptr<IMappingDatabaseAdapter> create(
        const std::string& instanceName,
        const std::string& jsonFilePath,
        std::shared_ptr<IParamMappingSource> loader
    ) const override
    {
        if (!loader) {
            loader = std::make_shared<JsonParamMappingAdapter>(std::make_shared<JSONParser>());
        }

        auto db = std::make_shared<MappingDatabase>(MappingDatabase(instanceName, jsonFilePath, loader));
        if (!db) throw std::runtime_error("DefaultMappingAdapterFactory: database is unavailable");

        return std::make_shared<DefaultMappingDatabaseAdapter>(instanceName, db);
    }
};

#endif
