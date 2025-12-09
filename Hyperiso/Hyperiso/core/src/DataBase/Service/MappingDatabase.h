#ifndef MAPPING_DATABASE_H
#define MAPPING_DATABASE_H

#include <unordered_map>
#include <string>
#include <memory>
#include "IParamMappingSource.h" 

/**
 * @file MappingDatabase.h
 * @brief In-memory store of parameter mappings loaded from an external file.
 *
 * This class encapsulates the logic needed to load a map of
 * (string name → InterpretedParam) from some external description,
 * using an IParamMappingSource implementation (e.g. JSON-based).
 */

/**
 * @class MappingDatabase
 * @brief Lightweight container for InterpretedParam mappings.
 *
 * A MappingDatabase is essentially a cache of the parameter mapping
 * file. It is responsible for:
 *   - invoking an IParamMappingSource to parse the mapping file,
 *   - storing the resulting associations in an in-memory map,
 *   - providing read-only access to the resulting map.
 *
 * It does *not* deal with higher-level concepts such as database
 * instances or multiple mapping sets; that logic is handled by
 * IMappingDatabaseAdapter implementations and factories.
 */
class MappingDatabase {
public:
    /**
     * @brief Constructs a MappingDatabase and immediately loads mappings.
     *
     * If @p loader is null, a default JSON-based loader (JsonParamMappingAdapter
     * + JSONParser) is instantiated. Any error during loading is reported to
     * stderr and results in an empty mapping.
     *
     * @param instanceName  Human-readable name of this mapping instance
     *                      (used for logging only).
     * @param jsonFilePath  Path to the mapping file to be loaded.
     * @param loader        Source object responsible for reading and
     *                      interpreting the mapping file.
     */
    MappingDatabase(const std::string& instanceName, const std::string& jsonFilePath,
        std::shared_ptr<IParamMappingSource> loader);
    
    /**
     * @brief Returns the current parameter mapping.
     *
     * The returned map is a copy of the internal storage to keep the
     * class logically immutable from the caller’s perspective.
     *
     * @return Map of parameter name → InterpretedParam.
     */
    std::unordered_map<std::string, InterpretedParam> getParams() const;

private:
    /**
     * @brief Helper that loads mappings from a file using the given loader.
     *
     * On failure, an error is printed to stderr and the internal map
     * is cleared.
     *
     * @param jsonFilePath  Path to the mapping file.
     * @param loader        Loader used to parse the file content.
     */
    void load(const std::string& jsonFilePath,
              const std::shared_ptr<IParamMappingSource>& loader);

private:
    /// In-memory map from logical parameter names to interpreted identifiers.
    std::unordered_map<std::string, InterpretedParam> paramsMap;
};

#endif // MAPPING_DATABASE_H
