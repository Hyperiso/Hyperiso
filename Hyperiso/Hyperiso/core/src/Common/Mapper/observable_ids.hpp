#ifndef OBSERVABLE_IDS_H
#define OBSERVABLE_IDS_H

#include "generic_mapper.h"
#include "Map.h"
#include "LhaID.h"
#include "decay_graph.h"

/**
 * @file ObservableMapper.h
 * @brief Mapper for Observables <-> text ids <-> FLHA ids, with decay links.
 *
 * This mapper combines:
 *   - the builtin mapping Observables <-> string names,
 *   - the FLHA mapping Observables <-> LhaID,
 *   - a dynamic registry for custom observables and aliases,
 *   - a link to DecayGraph to record which observables belong to which decay.
 */

struct ObservableTag {};
/** @brief Strongly typed identifier for observables. */
using ObservableId = IdOf<ObservableTag>;

/**
 * @class ObservableMapper
 * @brief High-level mapping for observables, including FLHA and decay associations.
 *
 * Inherits from GenericMapperWithExt with:
 *   - Tag         = ObservableTag
 *   - EnumT       = Observables
 *   - ExternalKey = LhaID (FLHA identifier)
 *   - Hash        = std::hash<LhaID>
 *   - MapFn       = observable_mapping
 *   - ExtFn       = observable_flha_mapping
 *
 * Provided functionality includes:
 *   - string (case-insensitive) <-> ObservableId,
 *   - Observables <-> ObservableId,
 *   - LhaID (FLHA) <-> ObservableId,
 *   - helper overloads enum_elt() from string or from FLHA,
 *   - registration of custom observables optionally attached to a decay.
 */
class ObservableMapper
: public GenericMapperWithExt<
      ObservableTag,
      Observables,
      LhaID,
      std::hash<LhaID>,
      observable_mapping,
      observable_flha_mapping
  >
{
public:
    using Base = GenericMapperWithExt<
        ObservableTag, Observables, LhaID, std::hash<LhaID>,
        observable_mapping, observable_flha_mapping
    >;

    using Base::init_builtins;
    using Base::id_of;           // string -> ObservableId
    using Base::str;             // str(ObservableId), str(Observables)
    using Base::to_id;           // Observables -> ObservableId
    using Base::enum_of;         // ObservableId -> optional<Observables>
    using Base::list_all;
    using Base::from_external;   // LhaID -> optional<ObservableId>
    using Base::external_of;     // ObservableId -> optional<LhaID>
    using Base::set_external;

    /**
     * @brief Legacy lookup: name -> Observables (builtin only).
     *
     * Performs a case-insensitive search over observable_mapping()
     * and returns the corresponding enum value.
     *
     * @param s Observable name.
     * @return Matching Observables enum value.
     *
     * @throws std::out_of_range if the name is unknown.
     */
    static Observables enum_elt(std::string_view s);

    /**
     * @brief Lookup from FLHA identifier to Observables enum (builtin only).
     *
     * This first resolves the FLHA id to an ObservableId via from_external(),
     * then maps that id to its builtin enum value via enum_of().
     *
     * @param ext FLHA identifier (LhaID).
     * @return Matching Observables enum value.
     *
     * @throws std::out_of_range if the FLHA id is unknown or not mapped
     *         to a builtin enum.
     */
    static Observables enum_elt(const LhaID& ext);

    /**
     * @brief Shortcut: FLHA -> ObservableId (dynamic + builtin).
     */
    static std::optional<ObservableId> from_flha(const LhaID& ext);

    /**
     * @brief Shortcut: ObservableId -> FLHA (dynamic + builtin).
     */
    static std::optional<LhaID>        flha_of(const ObservableId& id);

    /**
     * @brief Returns the FLHA id associated with a builtin Observables enum.
     *
     * This uses the static FLHA mapping observable_flha_mapping().
     *
     * @param id Observables enum value.
     * @return Corresponding LhaID.
     *
     * @throws std::out_of_range if @p id is not present in the mapping.
     */
    static LhaID        flha(const Observables& id);

    /**
     * @brief Registers a custom observable, optionally attached to a decay.
     *
     * The observable is inserted in the dynamic registry with:
     *   - canonical name,
     *   - optional aliases,
     *   - optional FLHA identifier.
     *
     * If @p parent_decay is provided and registration succeeds, a link
     * is also recorded in DecayGraph so that:
     *   - DecayGraph::observables_of(parent_decay) includes the new observable.
     *
     * @param canonical    Canonical name for the observable.
     * @param aliases      Optional aliases for name lookup.
     * @param ext          Optional FLHA identifier.
     * @param parent_decay Optional parent decay identifier.
     * @return true if registration succeeded, false if blocked by a builtin.
     */
    static bool register_custom(const std::string& canonical,
                            std::vector<std::string> aliases = {},
                            std::optional<LhaID> ext = std::nullopt,
                            std::optional<DecayId> parent_decay = std::nullopt);

    /**
     * @brief Returns the FLHA id associated with a given ObservableId.
     *
     * This queries the dynamic registry. If the observable has no FLHA
     * attached, an exception is thrown.
     *
     * @param id Observable identifier.
     * @return LhaID corresponding to @p id.
     *
     * @throws std::runtime_error if no FLHA id is defined for @p id.
     */
    static LhaID flha(const ObservableId& id);

    /**
     * @brief Convenience overload: register custom observable using a Decays enum.
     *
     * Same semantics as the generic register_custom(), but the parent decay
     * is specified via its enum value instead of a DecayId string.
     */
    static bool register_custom(const std::string& canonical,
                                std::vector<std::string> aliases,
                                std::optional<LhaID> ext,
                                Decays parent_decay_enum);

    /**
     * @brief Convenience overload: register custom observable using a decay name.
     *
     * Same semantics as the generic register_custom(), but the parent decay
     * is specified as a string (to be resolved to a DecayId) instead of a
     * DecayId or Decays enum.
     */
    static bool register_custom(const std::string& canonical,
                                std::vector<std::string> aliases,
                                std::optional<LhaID> ext,
                                std::string_view parent_decay_str);
};

#endif