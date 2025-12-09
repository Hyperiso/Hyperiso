#include "decay_graph.h"

DecayGraph& DecayGraph::instance() { static DecayGraph g; return g; }

void DecayGraph::link(const DecayId& decay, const ObservableId& obs) {
    std::lock_guard<std::mutex> lock(m_);
    auto& vec = graph_[decay.str()];
    if (std::find(vec.begin(), vec.end(), obs) == vec.end())
    vec.push_back(obs);
    parent_[obs.str()] = decay;
}

std::vector<ObservableId> DecayGraph::observables_of(const DecayId& decay) const {
    std::lock_guard<std::mutex> lock(m_);
    if (auto it = graph_.find(decay.str()); it != graph_.end()) return it->second;
    return {};
}

std::optional<DecayId> DecayGraph::parent_of(const ObservableId& obs) const {
    std::lock_guard<std::mutex> lock(m_);
    if (auto it = parent_.find(obs.str()); it != parent_.end()) return it->second;
    return std::nullopt;
}