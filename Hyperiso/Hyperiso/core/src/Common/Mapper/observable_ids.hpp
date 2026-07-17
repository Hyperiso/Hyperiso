#ifndef OBSERVABLE_IDS_H
#define OBSERVABLE_IDS_H

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "generic_mapper.h"
#include "Map.h"
#include "LhaID.h"
#include "decay_graph.h"

/**
 * @file observable_ids.hpp
 * @brief Mapper for builtin and runtime observable identifiers.
 *
 * The observable layer has two representations:
 *
 *   - `Observables`: the historical static enum, valid only for builtin
 *     observables present in `observable_mapping()`;
 *   - `ObservableId`: a string-based, case-insensitive identifier that can
 *     represent both builtin and custom observables.
 *
 * New code should use `ObservableId` for user-facing and dynamic paths.  The
 * static enum remains useful for legacy code and for places where a builtin
 * observable is explicitly required.
 */

struct ObservableTag {};

/** @brief Strongly typed dynamic identifier for observables. */
using ObservableId = IdOf<ObservableTag>;

struct BinnedObservableId;

/**
 * @class ObservableMapper
 * @brief High-level mapper for observable names, FLHA ids and parent decays.
 *
 * `ObservableMapper` combines:
 *   - builtin enum <-> string mappings from `observable_mapping()`,
 *   - optional FLHA ids from `observable_flha_mapping()`,
 *   - runtime registration through `DynamicRegistry`,
 *   - parent decay links through `DecayGraph`.
 *
 * @important A custom observable must always be registered with a parent decay.
 * This is enforced by the public `register_custom()` overloads below.
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
        ObservableTag,
        Observables,
        LhaID,
        std::hash<LhaID>,
        observable_mapping,
        observable_flha_mapping
    >;

    using Base::init_builtins;
    using Base::id_of;           ///< Dynamic string/alias -> ObservableId.
    using Base::enum_elt;        ///< Alias kept for generic code: string -> ObservableId.
    using Base::str;             ///< str(ObservableId) or str(Observables).
    using Base::to_id;           ///< Builtin Observables -> ObservableId.
    using Base::enum_of;         ///< ObservableId -> optional builtin Observables.
    using Base::list_all;
    using Base::from_external;   ///< LhaID -> optional ObservableId.
    using Base::external_of;     ///< ObservableId -> optional LhaID.
    using Base::set_external;

    /**
     * @brief Legacy builtin-only lookup from string to enum.
     *
     * This function deliberately ignores custom observables because a custom
     * observable cannot be represented by the static `Observables` enum.  Use
     * `id_of()` for all dynamic/user-facing lookup.
     *
     * @param s Builtin observable name or alias present in `observable_mapping()`.
     * @return Matching builtin enum value.
     *
     * @throws std::out_of_range if @p s is not a builtin observable.
     */
    static Observables enum_elt_legacy(std::string_view s);

    /**
     * @brief Legacy builtin-only lookup from FLHA id to enum.
     *
     * This first resolves the FLHA id to an `ObservableId`, then attempts to
     * recover the corresponding builtin enum.  Custom observables with a FLHA id
     * can be resolved through `from_flha()` but cannot necessarily be converted
     * to `Observables`.
     *
     * @param ext FLHA identifier.
     * @return Matching builtin enum value.
     *
     * @throws std::out_of_range if @p ext is unknown or points to a custom id.
     */
    static Observables enum_elt(const LhaID& ext);

    /**
     * @brief Resolve a FLHA id to a dynamic observable id.
     *
     * The canonical v1.0.2 polarization pattern is @c 92ij. For input
     * compatibility, the unambiguous v1.0.0--v1.0.1 polarization identifiers
     * @c 92015 and @c 921423 are accepted and translated to @c 9212 and
     * @c 9221, respectively. Output always uses the canonical identifiers.
     *
     * @param ext FLHA identifier.
     * @return Observable id if the FLHA id is known.
     */
    static std::optional<ObservableId> from_flha(const LhaID& ext);

    /**
     * @brief Return the FLHA id attached to an observable id, if any.
     *
     * @param id Observable identifier.
     * @return FLHA id if @p id has one.
     */
    static std::optional<LhaID> flha_of(const ObservableId& id);

    /**
     * @brief Decode a binned FLHA id into a binned observable id.
     *
     * Binned FLHA ids are normal observable FLHA ids where six integer fields
     * encoding the two bin boundaries are inserted after the first two fields.
     *
     * @param ext Binned FLHA identifier.
     * @return Decoded binned observable id, or nullopt if decoding fails.
     */
    static std::optional<BinnedObservableId> from_binned_flha(const LhaID& ext);

    /**
     * @brief Encode a binned observable id into its FLHA representation.
     *
     * @param id Binned observable identifier.
     * @return Encoded FLHA id, or nullopt if the observable has no FLHA id.
     */
    static std::optional<LhaID> binned_flha_of(const BinnedObservableId& id);

    /**
     * @brief Throwing variant of `from_binned_flha()`.
     *
     * @param ext Binned FLHA identifier.
     * @return Decoded binned observable id.
     *
     * @throws std::runtime_error if @p ext cannot be decoded.
     */
    static BinnedObservableId binned_from_flha(const LhaID& ext);

    /**
     * @brief Return the FLHA id associated with a builtin enum value.
     *
     * @param id Builtin observable enum.
     * @return FLHA id from `observable_flha_mapping()`.
     *
     * @throws std::out_of_range if @p id has no static FLHA entry.
     */
    static LhaID flha(const Observables& id);

    /**
     * @brief Register a custom observable and attach it to a parent decay.
     *
     * @param canonical Canonical observable name.
     * @param aliases Case-insensitive aliases for lookup.
     * @param ext Optional FLHA id.
     * @param parent_decay Parent decay id. It must already be registered, either
     *                     as builtin or custom.
     * @return true on success, false if registration is blocked by a builtin.
     *
     * @throws std::runtime_error if @p parent_decay is unknown.
     * @throws std::runtime_error if the observable is already linked to another
     *                            decay.
     */
    static bool register_custom(const std::string& canonical,
                                std::vector<std::string> aliases,
                                std::optional<LhaID> ext,
                                const DecayId& parent_decay);

    /**
     * @brief Register a custom observable under a builtin parent decay.
     */
    static bool register_custom(const std::string& canonical,
                                std::vector<std::string> aliases,
                                std::optional<LhaID> ext,
                                Decays parent_decay_enum);

    /**
     * @brief Register a custom observable under a parent decay resolved by name.
     */
    static bool register_custom(const std::string& canonical,
                                std::vector<std::string> aliases,
                                std::optional<LhaID> ext,
                                std::string_view parent_decay_str);

    /**
     * @brief Return the FLHA id associated with an observable id.
     *
     * @param id Observable id.
     * @return FLHA id.
     *
     * @throws std::runtime_error if @p id has no FLHA id.
     */
    static LhaID flha(const ObservableId& id);

    /**
     * @brief Return the binned FLHA id associated with a binned observable id.
     *
     * @param id Binned observable id.
     * @return Binned FLHA id.
     *
     * @throws std::runtime_error if encoding fails.
     */
    static LhaID flha(const BinnedObservableId& id);
};

#endif // OBSERVABLE_IDS_H
