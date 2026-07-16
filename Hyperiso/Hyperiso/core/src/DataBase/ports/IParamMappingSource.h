#ifndef IPARAMMAPPINGSOURCE_H
#define IPARAMMAPPINGSOURCE_H

#include <memory>
#include <string>
#include <unordered_map>

#include "Include.h"
#include "IParser.h"
#include "InterpretedParam.h"

/**
 * @file IParamMappingSource.h
 * @brief Interface for loading parameter mappings from an external source.
 *
 * Implementations of this interface know how to read some external
 * representation (JSON, YAML, database, etc.) and convert it into a
 * map of symbolic parameter names to InterpretedParam descriptors.
 */

/**
 * @class IParamMappingSource
 * @brief Abstract source of InterpretedParam mappings.
 *
 * A mapping source is responsible for:
 *   - locating and reading an input file (or other medium),
 *   - interpreting its structure,
 *   - returning a mapping from string names to InterpretedParam objects.
 *
 * Example implementations:
 *   - JsonParamMappingAdapter: loads from a JSON/YAML file via IParser.
 */
class IParamMappingSource {
public:
    virtual ~IParamMappingSource() = default;

    /**
     * @brief Loads all parameter mappings from a given file.
     *
     * The exact file format (JSON, YAML, etc.) is implementation-specific.
     *
     * @param filePath Path to the mapping file to load.
     * @return Map from parameter name (string) to InterpretedParam.
     *
     * @throws std::runtime_error or derived exceptions if the file cannot
     *         be read or parsed correctly.
     */
    virtual std::unordered_map<std::string, InterpretedParam>
    loadFromFile(const std::string& filePath) const = 0;
};

#endif // IPARAMMAPPINGSOURCE_H
