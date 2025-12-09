#ifndef WGROUP_IDS_H
#define WGROUP_IDS_H

#include "generic_mapper.h"
#include "Map.h"
#include "scaletype_ids.hpp"
#include "wilsonbasis_ids.hpp"

/**
 * @file wgroup_ids.h
 * @brief Mapping and helpers for Wilson operator groups.
 *
 * This header defines:
 *   - WGroupTag / WGroupId: strongly typed identifiers for Wilson groups,
 *   - a placeholder external mapping for WGroup (currently empty),
 *   - GroupMapper: high-level mapper WGroup <-> text <-> optional external string,
 *                  with helpers to build scale/basis-dependent block names.
 */

/**
 * @brief Placeholder for a WGroup -> external string mapping.
 *
 * At the moment no external string keys are defined for WGroup, so this
 * returns an empty map. It is provided to satisfy GenericMapperWithExt
 * and can be filled in later if a convention is adopted.
 */
inline const std::map<WGroup, std::string>& group_external_mapping(){
    static const std::map<WGroup, std::string> empty{};
    return empty;
}

/** @brief Tag type for Wilson-group identifiers. */
struct WGroupTag {};

/** @brief Strongly typed identifier for Wilson operator groups. */
using WGroupId = IdOf<WGroupTag>;

/**
 * @class GroupMapper
 * @brief Mapper for WGroup <-> WGroupId <-> optional external string.
 *
 * This class specializes GenericMapperWithExt with:
 *   - Tag         = WGroupTag
 *   - EnumT       = WGroup
 *   - ExternalKey = std::string (optional external label)
 *   - Hash        = std::hash<std::string>
 *   - MapFn       = group_mapping
 *   - ExtFn       = group_external_mapping
 *
 * It exposes:
 *   - the usual enum/string/id mapping utilities from GenericMapperWithExt,
 *   - a convenience register_custom() with an optional external string,
 *   - domain-specific helpers to build group block names including
 *     scale and Wilson basis information.
 */
class GroupMapper
: public GenericMapperWithExt<
      WGroupTag,
      WGroup,
      std::string,
      std::hash<std::string>,
      group_mapping,
      group_external_mapping
  >
{
public:
    using Base = GenericMapperWithExt<
        WGroupTag, WGroup, std::string, std::hash<std::string>,
        group_mapping, group_external_mapping
    >;

    using Base::init_builtins;
    using Base::id_of;
    using Base::list_all;
    using Base::from_external;
    using Base::external_of;
    using Base::to_id;          // WGroup   -> WGroupId
    using Base::enum_of;        // WGroupId -> optional<WGroup>
    using Base::str;            // overloads: str(WGroupId) and str(WGroup)
    using Base::enum_elt;       // string -> WGroupId (keeps new runtime id API)
    using Base::enum_elt_legacy;// string -> WGroup  (legacy enum API)

    /**
     * @brief Registers a custom Wilson group with optional external label.
     *
     * The external string, if provided, is stored as an additional key
     * that can be used with from_external() / external_of().
     *
     * @param canonical Canonical name for the group.
     * @param aliases   Optional aliases used for case-insensitive lookup.
     * @param ext       Optional external string key.
     * @return true if registration succeeded, false if blocked by a builtin.
     */
    static bool register_custom(std::string canonical,
                                std::vector<std::string> aliases = {},
                                std::optional<std::string> ext = std::nullopt)
    {
        return Base::register_custom(std::move(canonical), std::move(aliases), std::move(ext));
    }

    /**
     * @brief Builds a composite block name using a WGroupId, scale and basis.
     *
     * The format is:
     *   "<group>_<scale>[_<basis>]"
     *
     * where the Wilson basis is appended only for hadronic scale.
     *
     * @param gid   Group identifier.
     * @param s     Scale type (matching / hadronic).
     * @param b     Wilson basis (used only if s == HADRONIC).
     * @return Composite block name.
     */
    static std::string str(const WGroupId& gid, ScaleType s, WilsonBasis b = WilsonBasis::B_STANDARD){
        return str_impl(gid.str(), s, b);
    }

    /**
     * @brief Builds a composite block name using a WGroup enum, scale and basis.
     *
     * Same format as the overload taking WGroupId:
     *   "<group>_<scale>[_<basis>]"
     *
     * @param g     Wilson group enum.
     * @param s     Scale type (matching / hadronic).
     * @param b     Wilson basis (used only if s == HADRONIC).
     * @return Composite block name.
     */
    static std::string str(WGroup g, ScaleType s, WilsonBasis b = WilsonBasis::B_STANDARD){
        return str_impl(Base::str(g), s, b);
    }

private:
    /**
     * @brief Internal helper that composes "<base>_<scale>[_<basis>]".
     */
    static std::string str_impl(const std::string& base, ScaleType s, WilsonBasis b){
        std::string out = base + "_" + ScaleTypeMapper::str(s);
        if (s == ScaleType::HADRONIC) out += "_" + WilsonBasisMapper::str(b);
        return out;
    }
};

#endif