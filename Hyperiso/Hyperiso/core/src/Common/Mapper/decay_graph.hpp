// decay_graph.hpp
#pragma once
#include "dynamic_registry.hpp"   // for SymbolId
#include <unordered_map>
#include <vector>
#include <string>
#include <mutex>

// Forward declarations to break include cycle
struct DecayTag;
struct ObservableTag;

using DecayId      = SymbolId<DecayTag>;
using ObservableId = SymbolId<ObservableTag>;

class DecayGraph {
public:
    static DecayGraph& instance() { static DecayGraph g; return g; }

    void link(const DecayId& decay, const ObservableId& obs) {
        std::lock_guard<std::mutex> lock(m_);
        auto& vec = graph_[decay.str()];
        if (std::find(vec.begin(), vec.end(), obs) == vec.end())
        vec.push_back(obs);
        parent_[obs.str()] = decay;
    }

    std::vector<ObservableId> observables_of(const DecayId& decay) const {
        std::lock_guard<std::mutex> lock(m_);
        if (auto it = graph_.find(decay.str()); it != graph_.end()) return it->second;
        return {};
    }

    std::optional<DecayId> parent_of(const ObservableId& obs) const {
        std::lock_guard<std::mutex> lock(m_);
        if (auto it = parent_.find(obs.str()); it != parent_.end()) return it->second;
        return std::nullopt;
    }

private:
    mutable std::mutex m_;
    std::unordered_map<std::string, std::vector<ObservableId>> graph_;
    std::unordered_map<std::string, DecayId> parent_;
};
