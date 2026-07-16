#ifndef NUISANCEREADER_H
#define NUISANCEREADER_H

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_set>
#include <initializer_list>

#include "DBNode.h"
#include "DBNodeProviderFactory.h"
#include "INuisancePathsProvider.h"
#include "INuisanceReader.h"
#include "NuisanceSpec.h"

namespace fs = std::filesystem;

/**
 * @file NuisanceReader.h
 * @brief Concrete reader for nuisance-parameter definition files.
 *
 * @ref NuisanceReader loads nuisance specifications from DBNode-compatible
 * files, typically a default JSON registry plus an optional user YAML override.
 *
 * @see INuisanceReader
 * @see INuisancePathsProvider
 * @see NuisanceSpec
 */

/**
 * @class NuisanceReader
 * @brief Loads and parses nuisance-parameter specifications.
 *
 * This class reads nuisance entries from the paths supplied by an
 * @ref INuisancePathsProvider. The public @ref load function merges the default
 * registry with the user registry, where user entries override defaults with the
 * same @ref ParamId.
 *
 * Supported entry fields include:
 * - `block`: LHA block name,
 * - `code`: LHA code,
 * - `min_val` or `min`: lower bound,
 * - `max_val` or `max`: upper bound,
 * - `distribution` or `marginal`: marginal type.
 *
 * @throws std::invalid_argument if constructed with a null path provider.
 * @throws std::runtime_error for missing required fields, invalid bounds,
 *         unknown distributions, missing files, or malformed input nodes.
 */
class NuisanceReader : public INuisanceReader {
public:
    /**
     * @brief Constructs a nuisance reader with an external path provider.
     *
     * @param paths_provider Provider used to resolve default and user nuisance
     *        configuration files.
     *
     * @throws std::invalid_argument if @p paths_provider is null.
     */
    explicit NuisanceReader(std::shared_ptr<INuisancePathsProvider> paths_provider);

    /**
     * @copydoc INuisanceReader::load_default
     */
    NuisanceRegistry load_default() const override;

    /**
     * @copydoc INuisanceReader::load_user()
     */
    NuisanceRegistry load_user() const override;

    /**
     * @copydoc INuisanceReader::load_user(const fs::path&) const
     */
    NuisanceRegistry load_user(const fs::path& user_path) const override;

    /**
     * @brief Loads the merged nuisance registry.
     *
     * The default registry is loaded first. User-defined entries are then merged
     * on top and override default specifications with identical parameter ids.
     *
     * @return Merged nuisance registry.
     *
     * @throws std::runtime_error if either source cannot be loaded or parsed.
     */
    NuisanceRegistry load() const;

    /**
     * @brief Returns the configured default nuisance-file path.
     *
     * @return Path supplied by @ref INuisancePathsProvider::default_nuisances_path.
     */
    fs::path default_path() const;

    /**
     * @brief Returns the configured user nuisance-file path.
     *
     * @return Path supplied by @ref INuisancePathsProvider::user_nuisances_path.
     */
    fs::path user_path() const;

private:
    std::shared_ptr<INuisancePathsProvider> paths_provider_;    ///< Provider for default and user nuisance-file paths.

private:
    /**
     * @brief Reads one file and merges all contained entries into a registry.
     *
     * @param path File to parse through @ref DBNodeProviderFactory.
     * @param registry Registry updated in place. Existing entries with the same
     *        @ref ParamId are overwritten.
     *
     * @throws std::runtime_error if the file cannot be converted to a DB node or
     *         if the node content is malformed.
     */
    void merge_file_into_registry(const fs::path& path, NuisanceRegistry& registry) const;

    /**
     * @brief Extracts nuisance entries from a root DB node and merges them.
     *
     * @param root Root node expected to contain a `nuisances` entry.
     * @param registry Registry updated in place.
     *
     * @throws std::runtime_error if a `nuisances` field exists but no valid
     *         entries can be extracted.
     */
    void merge_node_into_registry(const DBNode& root, NuisanceRegistry& registry) const;

    /**
     * @brief Parses a single DB node into a nuisance specification.
     *
     * @param entry Node describing one nuisance parameter.
     * @return Parsed nuisance specification.
     *
     * @throws std::runtime_error if required fields are missing or invalid.
     */
    static NuisanceSpec parse_entry(const DBNode& entry);

    /**
     * @brief Builds the parameter id of a nuisance entry.
     *
     * @param entry Node containing `block` and `code` fields.
     * @return Parameter identifier built from the LHA block and code.
     */
    static ParamId make_param_id(const DBNode& entry);

    /**
     * @brief Returns the first available required value from a list of aliases.
     *
     * @param node Node to inspect.
     * @param candidate_keys Accepted field names, in priority order.
     * @return Value associated with the first matching key.
     *
     * @throws std::runtime_error if none of the candidate keys is present.
     */
    static DBNode::Value get_required_value(const DBNode& node,
                                            std::initializer_list<const char*> candidate_keys);

    /**
     * @brief Converts a DB value to a string representation.
     *
     * @param value Input DB value.
     * @param field_name Name used in error messages.
     * @return String representation with optional surrounding quotes removed.
     *
     * @throws std::runtime_error if @p value is not string-like.
     */
    static std::string value_to_string(const DBNode::Value& value,
                                       const std::string& field_name);

    /**
     * @brief Converts a DB value to a parameter-code string.
     *
     * @param value Input DB value.
     * @param field_name Name used in error messages.
     * @return Code string with optional surrounding quotes removed.
     */
    static std::string value_to_code_string(const DBNode::Value& value,
                                            const std::string& field_name);

    /**
     * @brief Converts a DB value to a floating-point number.
     *
     * @param value Input DB value.
     * @param field_name Name used in error messages.
     * @return Numeric value.
     *
     * @throws std::runtime_error if @p value cannot be interpreted as numeric.
     */
    static double value_to_double(const DBNode::Value& value,
                                  const std::string& field_name);

    /**
     * @brief Parses a textual marginal-distribution name.
     *
     * @param raw Raw distribution label. Case, spaces, hyphens, and optional
     *        quotes are normalized before matching.
     * @return Corresponding marginal type.
     *
     * @throws std::runtime_error if the distribution label is unknown.
     */
    static MarginalType parse_marginal_type(std::string raw);

    /**
     * @brief Normalizes a textual key or distribution label.
     *
     * @param s Input string.
     * @return Lower-case string with spaces and hyphens replaced by underscores.
     */
    static std::string normalise(std::string s);
};

#endif // NUISANCEREADER_H