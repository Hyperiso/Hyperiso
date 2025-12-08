#ifndef DECAY_GRAPH_H
#define DECAY_GRAPH_H

#include <vector>
#include <string>
#include <mutex>

#include "dynamic_registry.h"
#include <unordered_map>

/**
 * @file decay_graph.h
 * @brief Lightweight graph linking decay channels to observables.
 *
 * This header defines:
 *   - DecayId      : strong-typed string identifier for decays,
 *   - ObservableId : strong-typed string identifier for observables,
 *   - DecayGraph   : a singleton storing the relation
 *                    "decay channel → list of observables"
 *                    and the inverse "observable → parent decay".
 *
 * The graph is stored in memory and is thread-safe via an internal mutex.
 */

struct DecayTag;
struct ObservableTag;

/** @brief Strongly typed identifier for a decay channel. */
using DecayId      = SymbolId<DecayTag>;

/** @brief Strongly typed identifier for an observable. */
using ObservableId = SymbolId<ObservableTag>;

/**
 * @class DecayGraph
 * @brief Singleton managing decay–observable relationships.
 *
 * DecayGraph keeps two internal maps:
 *   - graph_  : decay name (string)  -> list of ObservableId
 *   - parent_ : observable name      -> DecayId (its parent decay)
 *
 * It is mainly used by ObservableMapper::register_custom to
 * associate newly registered observables to a given decay.
 *
 * All operations are thread-safe and protected by an internal mutex.
 */
class DecayGraph {
public:
    /**
     * @brief Returns the unique DecayGraph instance (singleton).
     */
    static DecayGraph& instance();

    /**
     * @brief Registers a link between a decay and an observable.
     *
     * If the observable is not yet present in the list for @p decay,
     * it is appended. The reverse mapping (parent decay of the observable)
     * is also updated.
     *
     * This function is thread-safe.
     *
     * @param decay Decay identifier.
     * @param obs   Observable identifier.
     */
    void link(const DecayId& decay, const ObservableId& obs);

    /**
     * @brief Returns all observables associated to a given decay.
     *
     * If the decay has no registered observables, an empty vector is returned.
     * The order corresponds to the insertion order.
     *
     * This function is thread-safe.
     *
     * @param decay Decay identifier.
     * @return Vector of ObservableId linked to @p decay.
     */
    std::vector<ObservableId> observables_of(const DecayId& decay) const;

    /**
     * @brief Returns the parent decay of a given observable, if any.
     *
     * If the observable is not known to the graph, std::nullopt is returned.
     *
     * This function is thread-safe.
     *
     * @param obs Observable identifier.
     * @return Optional DecayId corresponding to the parent decay.
     */
    std::optional<DecayId> parent_of(const ObservableId& obs) const;

private:
    mutable std::mutex m_;  ///< Protects concurrent access to the maps.
    std::unordered_map<std::string, std::vector<ObservableId>> graph_;
    std::unordered_map<std::string, DecayId> parent_;
};

#endif