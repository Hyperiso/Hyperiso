#ifndef JSONPARAMMAPPINGADAPTER_H
#define JSONPARAMMAPPINGADAPTER_H

#include <variant>
#include <iostream>

#include "IParamMappingSource.h"
#include "DBNode.h"

/**
 * @file JsonParamMappingAdapter.h
 * @brief Adapter to load parameter mappings from a JSON/YAML-like Node tree.
 *
 * This adapter uses an IParser implementation (typically JSONParser) to read
 * a mapping file, and exposes it as a collection of InterpretedParam entries
 * keyed by parameter name.
 */

/**
 * @class JsonParamMappingAdapter
 * @brief IParamMappingSource implementation backed by a JSON-like file.
 *
 * Expected JSON (or YAML) schema:
 * @code
 * {
 *   "paramName1": {
 *     "block":   "<blockName>",
 *     "pdgCode": "<LHA ID string or integer-like>"
 *   },
 *   "paramName2": { ... }
 * }
 * @endcode
 *
 * For each top-level key, the adapter constructs an InterpretedParam by:
 *   - reading "block"  as a string,
 *   - reading "pdgCode" and converting it to an LhaID.
 */
class JsonParamMappingAdapter : public IParamMappingSource {
public:
    /**
     * @brief Constructs the adapter with a given parser instance.
     *
     * @param parser Concrete parser (e.g. JSONParser) used to load the file
     *               into a Node tree.
     */
    explicit JsonParamMappingAdapter(std::shared_ptr<IParser> parser)
        : parser_(std::move(parser)) {}

    /**
     * @brief Loads parameter mappings from a file.
     *
     * @param filePath Path to the JSON/YAML mapping file.
     * @return Map from parameter name to InterpretedParam.
     *
     * @throws std::runtime_error if the parser is null, the root Node is null,
     *         or the file cannot be parsed.
     */
    std::unordered_map<std::string, InterpretedParam>
    loadFromFile(const std::string& filePath) const override;

private:
    /**
     * @brief Converts a scalar Node::Value to a string.
     *
     * Supported variants:
     *   - BlockName -> its string alias,
     *   - int       -> std::to_string(i),
     *   - double    -> std::to_string(d),
     *   - bool      -> "true" / "false".
     *
     * Any non-scalar variant triggers a runtime error.
     *
     * @param v Node::Value to convert.
     * @return String representation of the value.
     */
    static std::string asString(const Node::Value& v) {
        if (auto s = std::get_if<BlockName>(&v)) return *s;
        if (auto i = std::get_if<int>(&v))      return std::to_string(*i);
        if (auto d = std::get_if<double>(&v))   return std::to_string(*d);
        if (auto b = std::get_if<bool>(&v))     return *b ? "true" : "false";
        throw std::runtime_error("JSON mapping: valeur non scalaire convertissable en string");
    }

private:
    /// Underlying parser used to load the JSON/YAML mapping file.
    std::shared_ptr<IParser> parser_;
};

#endif // JSONPARAMMAPPINGADAPTER_H
