#ifndef DECAY_IDS_H
#define DECAY_IDS_H

#include "generic_mapper.h"
#include "Map.h"
#include "LhaID.h" 
#include "observable_ids.hpp"

/**
 * @file decay_ids.h
 * @brief Mapper and identifiers for decay channels.
 *
 * This header defines:
 *   - DecayTag / DecayId: strongly typed identifiers for decays,
 *   - a trivial external mapping placeholder for decays (FLHA-like),
 *   - DecayMapper: high-level mapping Decays <-> text <-> LhaID,
 *                  plus helpers to relate decays and observables.
 */

/** @brief Tag type for decay identifiers. */

struct DecayTag {};

/** @brief Strongly typed identifier for decays (string-based). */
using DecayId = IdOf<DecayTag>;

/**
 * @brief Placeholder for a Decays -> LhaID mapping.
 *
 * For now, decays do not have a builtin external (FLHA) mapping, so this
 * function returns an empty map. It exists to satisfy the GenericMapperWithExt
 * interface and can be extended later when a proper external mapping is added.
 */
inline const std::map<Decays, LhaID>& decay_external_mapping(){
    static const std::map<Decays, LhaID> empty{};
    return empty;
}

/**
 * @class DecayMapper
 * @brief Mapper for Decays <-> DecayId <-> external LhaID (optional).
 *
 * This class specializes GenericMapperWithExt with:
 *   - Tag         = DecayTag
 *   - EnumT       = Decays
 *   - ExternalKey = LhaID
 *   - Hash        = std::hash<LhaID>
 *   - MapFn       = decays_mapping
 *   - ExtFn       = decay_external_mapping
 *
 * In addition to the generic mapping features, it provides:
 *   - helpers to retrieve observables associated to a decay enum,
 *   - helpers to find the decay associated to a given observable,
 *     combining dynamic information from DecayGraph and a legacy
 *     static cache built from decay_observable_mapping().
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
        DecayTag, Decays, LhaID, std::hash<LhaID>,
        decays_mapping, decay_external_mapping
    >;

    using Base::init_builtins;
    using Base::id_of;
    using Base::str;         // str(DecayId), str(Decays)
    using Base::to_id;       // Decays -> DecayId
    using Base::enum_of;
    using Base::list_all;
    using Base::from_external;   // LhaID -> optional<DecayId>
    using Base::external_of;     // DecayId -> optional<LhaID>
    using Base::set_external; 

    /**
     * @brief Returns the list of observables associated with a decay enum.
     *
     * Uses the static mapping decay_observable_mapping() and throws if the
     * decay is not present.
     *
     * @param d Decay enum value.
     * @return Vector of Observables belonging to this decay channel.
     */
    static std::vector<Observables> get_observables(Decays d) {
        return decay_observable_mapping().at(d);
    }

    /**
     * @brief Returns the decay enum associated with a given observable enum.
     *
     * Scans decay_observable_mapping() and returns the first decay for which
     * the observable appears in the associated vector.
     *
     * @param obs Observable enum value.
     * @return Decays enum value corresponding to the parent decay.
     *
     * @throws std::runtime_error if the observable is not mapped to any decay.
     */
    static Decays get_decay(Observables obs) {
        for (const auto& [d, vec] : decay_observable_mapping())
            if (std::find(vec.begin(), vec.end(), obs) != vec.end()) return d;
        throw std::runtime_error("Observable not mapped to any decay");
    }

    /**
     * @brief Lazily built cache ObservableId -> DecayId from static mapping.
     *
     * The first time this is called, it constructs a map by iterating over
     * decay_observable_mapping(), converting each Decays and Observables
     * entry into DecayId and ObservableId via DecayMapper::to_id and
     * ObservableMapper::to_id.
     *
     * @return Reference to the static cache map.
     */
    inline static const std::unordered_map<ObservableId, DecayId>& observable_to_decay() {
        static std::unordered_map<ObservableId, DecayId> cache;
        static std::once_flag init;
        std::call_once(init, []{
            for (const auto& [decayEnum, vec] : decay_observable_mapping()) {
                DecayId did = DecayMapper::to_id(decayEnum);
                for (auto o : vec)
                    cache.emplace(ObservableMapper::to_id(o), did);
            }
        });
        return cache;
    }

    /**
     * @brief Returns the DecayId associated with a given ObservableId, if any.
     *
     * Resolution is done in two steps:
     *   1. Query DecayGraph (dynamic links created at runtime via
     *      ObservableMapper::register_custom). If a parent decay is found,
     *      it is returned.
     *   2. Fallback to the legacy static cache built from
     *      decay_observable_mapping() via observable_to_decay().
     *
     * If neither source provides a match, std::nullopt is returned.
     *
     * @param obs Observable identifier.
     * @return Optional DecayId of the parent decay.
     */
    inline static std::optional<DecayId> get_decay_id(ObservableId obs) {
        if (auto d = DecayGraph::instance().parent_of(obs); d)  // dynamique
            return d;

        // Fallback legacy (cache construit depuis decay_observable_mapping()) :
        const auto& m = observable_to_decay();
        if (auto it = m.find(obs); it != m.end()) return it->second;
        return std::nullopt;
    }

};

#endif