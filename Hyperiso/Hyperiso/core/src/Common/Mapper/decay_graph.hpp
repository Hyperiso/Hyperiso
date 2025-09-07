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
        graph_[decay.str()].push_back(obs);
    }

    std::vector<ObservableId> observables_of(const DecayId& decay) const {
        std::lock_guard<std::mutex> lock(m_);
        auto it = graph_.find(decay.str());
        if (it == graph_.end()) return {};
        return it->second;
    }

private:
    mutable std::mutex m_;
    std::unordered_map<std::string, std::vector<ObservableId>> graph_;
};
