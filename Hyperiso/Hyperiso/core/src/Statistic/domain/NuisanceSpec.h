#ifndef NUISANCESPEC_H
#define NUISANCESPEC_H

#include <ostream>
#include <utility>
#include <unordered_map>

#include "ParamID.h"
#include "MarginalType.h"

/**
 * @file NuisanceSpec.h
 * @brief Data structures for nuisance-parameter specifications.
 *
 * Nuisance specifications associate an LHA parameter id with allowed bounds and
 * the marginal distribution used when building nuisance constraints.
 */

/**
 * @struct NuisanceSpec
 * @brief Specification of one nuisance parameter.
 *
 * A nuisance specification is the parsed representation of one nuisance entry in
 * a configuration file. It is keyed by @ref ParamId inside a
 * @ref NuisanceRegistry.
 */
struct NuisanceSpec {
    ParamId param_id;                                   ///< Parameter identifier, usually an LHA block/code pair.
    std::pair<double, double> bounds;                   ///< Inclusive lower and upper bounds for the nuisance value.
    MarginalType marginal = MarginalType::GAUSSIAN;     ///< Marginal model used for the nuisance constraint.
};

/**
 * @brief Registry of nuisance specifications indexed by parameter id.
 */
using NuisanceRegistry = std::unordered_map<ParamId, NuisanceSpec>;

/**
 * @brief Converts a marginal type to a stable diagnostic string.
 *
 * @param type Marginal type to stringify.
 *
 * @return String literal representing @p type, or `"UNKNOWN"` for unrecognized
 *         values.
 */
inline const char* to_string(MarginalType type) {
    switch (type) {
        case MarginalType::GAUSSIAN:      return "GAUSSIAN";
        case MarginalType::HALF_GAUSSIAN: return "HALF_GAUSSIAN";
        case MarginalType::FLAT:          return "FLAT";
        case MarginalType::LIKELIHOOD:    return "LIKELIHOOD";
        default:                          return "UNKNOWN";
    }
}

/**
 * @brief Streams a compact representation of one nuisance specification.
 *
 * @param os Output stream.
 * @param spec Nuisance specification to print.
 *
 * @return Reference to @p os.
 */
inline std::ostream& operator<<(std::ostream& os, const NuisanceSpec& spec) {
    os << "NuisanceSpec{"
       << "param_id=" << spec.param_id
       << ", bounds=(" << spec.bounds.first << ", " << spec.bounds.second << ")"
       << ", marginal=" << to_string(spec.marginal)
       << "}";

    return os;
}

/**
 * @struct NuisanceRegistryPrintable
 * @brief Lightweight wrapper enabling pretty-printing of a nuisance registry.
 *
 * Use @ref as_printable to create this wrapper without copying the registry.
 */
struct NuisanceRegistryPrintable {
    const NuisanceRegistry& registry;   ///< Registry referenced for formatted output.
};

/**
 * @brief Wraps a registry for formatted stream output.
 *
 * @param registry Registry to print.
 *
 * @return Non-owning printable wrapper around @p registry.
 */
inline NuisanceRegistryPrintable as_printable(const NuisanceRegistry& registry) {
    return NuisanceRegistryPrintable{registry};
}

/**
 * @brief Streams a multi-line representation of a nuisance registry.
 *
 * @param os Output stream.
 * @param printable Printable registry wrapper.
 *
 * @return Reference to @p os.
 */
inline std::ostream& operator<<(std::ostream& os, const NuisanceRegistryPrintable& printable) {
    const auto& registry = printable.registry;

    os << "NuisanceRegistry{\n";
    for (const auto& [param_id, spec] : registry) {
        os << "  [" << param_id << "] -> " << spec << "\n";
    }
    os << "}";

    return os;
}

#endif // NUISANCESPEC_H