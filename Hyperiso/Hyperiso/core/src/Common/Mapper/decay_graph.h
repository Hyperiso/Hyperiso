#ifndef DECAY_GRAPH_H
#define DECAY_GRAPH_H

#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "dynamic_registry.h"

/**
 * @file decay_graph.h
 * @brief Runtime graph linking decay identifiers to observable identifiers.
 *
 * The static enum mappers (`Decays`, `Observables`) are not sufficient once the
 * user can register new decays and observables at runtime.  `DecayGraph` stores
 * the dynamic relation:
 *
 *   - one decay channel (`DecayId`) owns zero or more observables
 *     (`ObservableId`),
 *   - one observable has at most one parent decay.
 *
 * The graph is intentionally small and in-memory only.  It complements the
 * static `decay_observable_mapping()` table used for builtin decays.
 *
 * @note Keys are normalized with `normalize_key()` so graph lookups follow the
 * same case-insensitive semantics as `DynamicRegistry`.
 */

struct DecayTag;
struct ObservableTag;

/** @brief Strongly typed string identifier for a decay channel. */
using DecayId = SymbolId<DecayTag>;

/** @brief Strongly typed string identifier for an observable. */
using ObservableId = SymbolId<ObservableTag>;

/**
 * @class DecayGraph
 * @brief Thread-safe singleton storing dynamic decay-observable links.
 *
 * `DecayGraph` is used when custom symbols are registered.  Builtin relations
 * still live in `decay_observable_mapping()`, while this graph stores relations
 * that cannot be represented by the static enum maps.
 *
 * Invariants enforced by this class:
 *   - an observable can only have one parent decay,
 *   - duplicate links are ignored,
 *   - lookup keys are case-insensitive.
 */
class DecayGraph {
public:
    /**
     * @brief Return the unique global graph instance.
     *
     * @return Mutable singleton instance.
     */
    static DecayGraph& instance();

    /**
     * @brief Link an observable to a parent decay.
     *
     * If the link already exists, this function is a no-op.  If the observable
     * is already linked to a different decay, an exception is thrown because an
     * observable must not belong to several decay channels.
     *
     * @param decay Parent decay identifier.
     * @param obs Observable identifier to attach to @p decay.
     *
     * @throws std::runtime_error if @p obs is already linked to another decay.
     */
    void link(const DecayId& decay, const ObservableId& obs);

    /**
     * @brief Return all dynamically linked observables for a decay.
     *
     * This only returns runtime links stored in the graph.  Builtin static
     * observables are merged at `DecayMapper` level.
     *
     * @param decay Decay identifier.
     * @return Observables dynamically linked to @p decay, in insertion order.
     */
    std::vector<ObservableId> observables_of(const DecayId& decay) const;

    /**
     * @brief Return the dynamic parent decay of an observable, if known.
     *
     * This only checks runtime links.  Builtin static fallback is handled by
     * `DecayMapper::get_decay_id()`.
     *
     * @param obs Observable identifier.
     * @return Parent decay if @p obs is known by the dynamic graph.
     */
    std::optional<DecayId> parent_of(const ObservableId& obs) const;

    /**
     * @brief Check whether a decay has at least one dynamic observable.
     *
     * @param decay Decay identifier.
     * @return true if @p decay has at least one runtime-linked observable.
     */
    bool has_observables(const DecayId& decay) const;

private:
    mutable std::mutex m_;  ///< Protects all maps below.

    /** @brief Normalized decay name -> dynamic observable list. */
    std::unordered_map<std::string, std::vector<ObservableId>> graph_;

    /** @brief Normalized observable name -> dynamic parent decay. */
    std::unordered_map<std::string, DecayId> parent_;
};

#endif // DECAY_GRAPH_H
