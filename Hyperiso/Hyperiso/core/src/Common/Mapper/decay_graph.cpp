#include "decay_graph.h"

#include <algorithm>
#include <stdexcept>

namespace {

/**
 * @brief Normalized key used for decay lookup in DecayGraph.
 */
std::string decay_key(const DecayId& decay) {
    return normalize_key(decay.str());
}

/**
 * @brief Normalized key used for observable lookup in DecayGraph.
 */
std::string observable_key(const ObservableId& obs) {
    return normalize_key(obs.str());
}

} // namespace

DecayGraph& DecayGraph::instance() {
    static DecayGraph g;
    return g;
}

void DecayGraph::link(const DecayId& decay, const ObservableId& obs) {
    std::lock_guard<std::mutex> lock(m_);

    const std::string dk = decay_key(decay);
    const std::string ok = observable_key(obs);

    if (auto it = parent_.find(ok); it != parent_.end() && !(it->second == decay)) {
        throw std::runtime_error(
            "Observable '" + obs.str() +
            "' is already linked to decay '" + it->second.str() + "'"
        );
    }

    auto& vec = graph_[dk];
    if (std::find(vec.begin(), vec.end(), obs) == vec.end()) {
        vec.push_back(obs);
    }

    parent_[ok] = decay;
}

std::vector<ObservableId> DecayGraph::observables_of(const DecayId& decay) const {
    std::lock_guard<std::mutex> lock(m_);

    if (auto it = graph_.find(decay_key(decay)); it != graph_.end()) {
        return it->second;
    }

    return {};
}

std::optional<DecayId> DecayGraph::parent_of(const ObservableId& obs) const {
    std::lock_guard<std::mutex> lock(m_);

    if (auto it = parent_.find(observable_key(obs)); it != parent_.end()) {
        return it->second;
    }

    return std::nullopt;
}

bool DecayGraph::has_observables(const DecayId& decay) const {
    std::lock_guard<std::mutex> lock(m_);

    if (auto it = graph_.find(decay_key(decay)); it != graph_.end()) {
        return !it->second.empty();
    }

    return false;
}
