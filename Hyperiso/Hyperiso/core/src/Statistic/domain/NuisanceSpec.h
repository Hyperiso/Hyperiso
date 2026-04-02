#ifndef NUISANCESPEC_H
#define NUISANCESPEC_H

#include <ostream>
#include <utility>
#include <unordered_map>

#include "ParamID.h"
#include "MarginalType.h"


struct NuisanceSpec {
    ParamId param_id;
    std::pair<double, double> bounds;
    MarginalType marginal = MarginalType::GAUSSIAN;
};

using NuisanceRegistry = std::unordered_map<ParamId, NuisanceSpec>;

inline const char* to_string(MarginalType type) {
    switch (type) {
        case MarginalType::GAUSSIAN:      return "GAUSSIAN";
        case MarginalType::HALF_GAUSSIAN: return "HALF_GAUSSIAN";
        case MarginalType::FLAT:          return "FLAT";
        case MarginalType::LIKELIHOOD:    return "LIKELIHOOD";
        default:                          return "UNKNOWN";
    }
}

inline std::ostream& operator<<(std::ostream& os, const NuisanceSpec& spec) {
    os << "NuisanceSpec{"
       << "param_id=" << spec.param_id
       << ", bounds=(" << spec.bounds.first << ", " << spec.bounds.second << ")"
       << ", marginal=" << to_string(spec.marginal)
       << "}";

    return os;
}

struct NuisanceRegistryPrintable {
    const NuisanceRegistry& registry;
};

inline NuisanceRegistryPrintable as_printable(const NuisanceRegistry& registry) {
    return NuisanceRegistryPrintable{registry};
}

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