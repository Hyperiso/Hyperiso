#ifndef DECAY_IDS_H
#define DECAY_IDS_H

#include <algorithm>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "generic_mapper.h"
#include "Map.h"
#include "LhaID.h"
#include "observable_ids.hpp"

/**
 * @file decay_ids.hpp
 * @brief Mapper and dynamic identifiers for decay channels.
 *
 * The decay layer mirrors the observable layer:
 *
 *   - `Decays` is the historical static enum for builtin decay channels;
 *   - `DecayId` is the string-based dynamic identifier used for both builtin
 *     and custom decays.
 *
 * Static relations (`Decays -> Observables`) are still read from
 * `decay_observable_mapping()`. Runtime relations (`DecayId -> ObservableId`)
 * are stored in `DecayGraph`.
 */

/** @brief Tag type for decay identifiers. */
struct DecayTag {};

/** @brief Strongly typed dynamic identifier for decays. */
using DecayId = IdOf<DecayTag>;

/**
 * @brief Placeholder external mapping for decays.
 *
 * Decays currently have no builtin FLHA-like external identifier.  The empty
 * map keeps `DecayMapper` compatible with `GenericMapperWithExt`, so a real
 * external mapping can be added later without changing the public mapper shape.
 *
 * @return Empty static map.
 */
inline const std::map<Decays, LhaID>& decay_external_mapping() {
    static const std::map<Decays, LhaID> empty{};
    return empty;
}

/**
 * @struct CustomObservableSpec
 * @brief Small value object used to register a custom decay with observables.
 *
 * This helper avoids temporarily creating an empty custom decay when the caller
 * already knows which custom observables belong to it.
 */
struct CustomObservableSpec {
    std::string canonical;                  ///< Canonical observable name.
    std::vector<std::string> aliases{};     ///< Optional aliases.
    std::optional<LhaID> ext = std::nullopt; ///< Optional FLHA id.
};

/**
 * @class DecayMapper
 * @brief Mapper for decay names, dynamic ids and observable membership.
 *
 * `DecayMapper` resolves names to `DecayId` through the dynamic registry and
 * provides both legacy static accessors and dynamic `DecayId`-based accessors.
 *
 * @important For dynamic/user-facing code, prefer `DecayId` and
 * `ObservableId`.  Use `Decays`/`Observables` only when a builtin enum value is
 * explicitly required.
 */
class DecayMapper
: public GenericMapperWithExt<
      DecayTag,
      Decays,
      LhaID,
      std::hash<LhaID>,
      decays_mapping,
      decay_external_mapping
  >
{
public:
    using Base = GenericMapperWithExt<
        DecayTag,
        Decays,
        LhaID,
        std::hash<LhaID>,
        decays_mapping,
        decay_external_mapping
    >;

    using Base::init_builtins;
    using Base::id_of;           ///< Dynamic string/alias -> DecayId.
    using Base::enum_elt;        ///< Alias kept for generic code: string -> DecayId.
    using Base::str;             ///< str(DecayId) or str(Decays).
    using Base::to_id;           ///< Builtin Decays -> DecayId.
    using Base::enum_of;         ///< DecayId -> optional builtin Decays.
    using Base::list_all;
    using Base::from_external;
    using Base::external_of;
    using Base::set_external;

    /**
     * @brief Disabled low-level decay registration.
     *
     * A decay without observables violates the public invariant of the dynamic
     * mapper layer.  Use `register_custom_with_observables()` to create a
     * custom decay together with at least one custom observable.
     */
    static bool register_custom(std::string canonical,
                                std::vector<std::string> aliases = {},
                                std::optional<LhaID> ext = std::nullopt) = delete;

    /**
     * @brief Legacy static lookup of builtin observables for a builtin decay.
     *
     * @param d Builtin decay enum value.
     * @return Builtin observable enum values attached to @p d.
     *
     * @throws std::out_of_range if @p d is not present in the static map.
     */
    static std::vector<Observables> get_observables(Decays d) {
        return decay_observable_mapping().at(d);
    }

    /**
     * @brief Convert builtin observables of a builtin decay to dynamic ids.
     *
     * @param d Builtin decay enum value.
     * @return Observable ids corresponding to the builtin observables.
     */
    static std::vector<ObservableId> get_observable_ids(Decays d) {
        std::vector<ObservableId> out;
        const auto& obs_enums = decay_observable_mapping().at(d);
        out.reserve(obs_enums.size());

        for (auto obs : obs_enums) {
            out.push_back(ObservableMapper::to_id(obs));
        }

        return out;
    }

    /**
     * @brief Dynamic lookup of observables attached to a decay id.
     *
     * This merges:
     *   - runtime links stored in `DecayGraph`,
     *   - static builtin links if @p d corresponds to a builtin `Decays` enum.
     *
     * Duplicates are removed while preserving insertion order.
     *
     * @param d Dynamic decay identifier.
     * @return Observable ids attached to @p d.
     *
     * @throws std::runtime_error if @p d has no observable.
     */
    static std::vector<ObservableId> get_observables(const DecayId& d) {
        std::vector<ObservableId> out;
        std::unordered_set<ObservableId> seen;

        auto add = [&](const ObservableId& obs) {
            if (seen.insert(obs).second) {
                out.push_back(obs);
            }
        };

        for (const auto& obs : DecayGraph::instance().observables_of(d)) {
            add(obs);
        }

        if (auto builtin = enum_of(d); builtin) {
            for (const auto& obs : get_observable_ids(*builtin)) {
                add(obs);
            }
        }

        if (out.empty()) {
            throw std::runtime_error("Decay has no observable: " + d.str());
        }

        return out;
    }

    /**
     * @brief Dynamic lookup of observables attached to a decay resolved by name.
     *
     * @param d Decay name or alias.
     * @return Observable ids attached to the resolved decay.
     */
    static std::vector<ObservableId> get_observables(std::string_view d) {
        return get_observables(id_of(d));
    }

    /**
     * @brief Test whether a decay has at least one observable.
     *
     * This checks both the runtime graph and the builtin static mapping.
     *
     * @param d Dynamic decay identifier.
     * @return true if @p d has at least one observable.
     */
    static bool has_observables(const DecayId& d) {
        if (DecayGraph::instance().has_observables(d)) {
            return true;
        }

        if (auto builtin = enum_of(d); builtin) {
            auto it = decay_observable_mapping().find(*builtin);
            return it != decay_observable_mapping().end() && !it->second.empty();
        }

        return false;
    }

    /**
     * @brief Throw if a decay has no observable.
     *
     * @param d Dynamic decay identifier.
     *
     * @throws std::runtime_error if @p d has no observable.
     */
    static void require_observables(const DecayId& d) {
        if (!has_observables(d)) {
            throw std::runtime_error("Decay must have at least one observable: " + d.str());
        }
    }

    /**
     * @brief Register a custom decay and its custom observables in one call.
     *
     * This helper enforces the invariant that a custom decay must have at least
     * one observable.  It first registers the decay, then registers every custom
     * observable under that decay through `ObservableMapper::register_custom()`.
     *
     * @param canonical Canonical decay name.
     * @param aliases Optional decay aliases.
     * @param observables Custom observable definitions to attach.
     * @return true if the decay was registered and observable registrations were
     *         attempted successfully; false if the decay registration is blocked
     *         by a builtin symbol.
     *
     * @throws std::invalid_argument if @p observables is empty.
     * @throws std::runtime_error if one observable cannot be linked.
     *
     * @note Registration is not transactional because the underlying registry has
     * no unregister operation.  Avoid reusing canonical names on failure.
     */
    static bool register_custom_with_observables(std::string canonical,
                                                 std::vector<std::string> aliases,
                                                 std::vector<CustomObservableSpec> observables)
    {
        if (observables.empty()) {
            throw std::invalid_argument(
                "Cannot register custom decay '" + canonical + "' without observables"
            );
        }

        const bool ok = Base::register_custom(canonical, std::move(aliases), std::nullopt);
        if (!ok) {
            return false;
        }

        const DecayId decay_id = Base::id_of(canonical);

        for (auto& obs : observables) {
            ObservableMapper::register_custom(
                obs.canonical,
                std::move(obs.aliases),
                std::move(obs.ext),
                decay_id
            );
        }

        return true;
    }

    /**
     * @brief Legacy static lookup of the builtin decay owning a builtin observable.
     *
     * @param obs Builtin observable enum.
     * @return Builtin parent decay enum.
     *
     * @throws std::runtime_error if @p obs is not mapped statically.
     */
    static Decays get_decay(Observables obs) {
        for (const auto& [d, vec] : decay_observable_mapping()) {
            if (std::find(vec.begin(), vec.end(), obs) != vec.end()) {
                return d;
            }
        }
        throw std::runtime_error("Observable not mapped to any decay");
    }

    /**
     * @brief Static cache mapping builtin observable ids to builtin decay ids.
     *
     * @return Map built once from `decay_observable_mapping()`.
     */
    inline static const std::unordered_map<ObservableId, DecayId>& observable_to_decay() {
        static std::unordered_map<ObservableId, DecayId> cache;
        static std::once_flag init;

        std::call_once(init, [] {
            for (const auto& [decayEnum, vec] : decay_observable_mapping()) {
                DecayId did = DecayMapper::to_id(decayEnum);
                for (auto o : vec) {
                    cache.emplace(ObservableMapper::to_id(o), did);
                }
            }
        });

        return cache;
    }

    /**
     * @brief Return the parent decay id of an observable id.
     *
     * Resolution order:
     *   1. runtime graph (`DecayGraph`) for custom links,
     *   2. static cache built from builtin mappings.
     *
     * @param obs Observable identifier.
     * @return Parent decay id, if known.
     */
    inline static std::optional<DecayId> get_decay_id(const ObservableId& obs) {
        if (auto d = DecayGraph::instance().parent_of(obs); d) {
            return d;
        }

        const auto& m = observable_to_decay();
        if (auto it = m.find(obs); it != m.end()) {
            return it->second;
        }

        return std::nullopt;
    }

    /**
     * @brief Return the parent decay id of a builtin observable enum.
     */
    inline static std::optional<DecayId> get_decay_id(Observables obs) {
        return get_decay_id(ObservableMapper::to_id(obs));
    }

    /**
     * @brief Throwing variant of `get_decay_id()`.
     *
     * @param obs Observable identifier.
     * @return Parent decay id.
     *
     * @throws std::runtime_error if @p obs has no parent decay.
     */
    inline static DecayId get_decay_id_or_throw(const ObservableId& obs) {
        auto d = get_decay_id(obs);
        if (!d) {
            throw std::runtime_error("Observable not mapped to any decay: " + obs.str());
        }
        return *d;
    }
};

#endif // DECAY_IDS_H
